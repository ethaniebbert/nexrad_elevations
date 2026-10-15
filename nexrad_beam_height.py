#!/usr/bin/env python
"""
Lowest-tilt beam height above ground level (AGL) for every U.S. NEXRAD (WSR-88D)
radar, computed against the Copernicus GLO-90 digital surface model.

For each radar the beam is traced along a polar grid using the standard
effective-earth-radius (4/3) refraction model (Doviak & Zrnic 1993, eq. 2.28):

    h(r) = sqrt(r^2 + (k a)^2 + 2 r (k a) sin(theta)) - k a + h_antenna
    s(r) = k a * asin(r cos(theta) / (k a + h(r) - h_antenna))

Terrain is sampled at every range gate, partial beam blockage (PBB) is computed
for a circular half-power beam cross-section (Bech et al. 2003), and the
cumulative beam blockage (CBB) is the running maximum of PBB along each radial.
Points downrange of >= --block-threshold blockage are flagged as terrain-blocked.

Results are resampled to a regular lat/lon grid, composited across radars
(lowest unblocked beam wins), and written as:

  * nexrad_beam_height_agl.kmz  - Google Earth overlays + radar placemarks
  * <region>_beam_height_agl.tif - GeoTIFFs (meters) for GIS work
  * radar_summary.csv           - per-radar statistics

Usage:
    python nexrad_beam_height.py                     # every U.S. radar
    python nexrad_beam_height.py --sites KABX,KFSX   # a subset
    python nexrad_beam_height.py --help
"""

from __future__ import annotations

import os
import sys

# Point PROJ at this environment's database. A PROJ_LIB inherited from another
# conda install otherwise breaks pyproj/rasterio CRS lookups.
_PROJ_DIR = os.path.join(sys.prefix, "share", "proj")
if os.path.exists(os.path.join(_PROJ_DIR, "proj.db")):
    os.environ["PROJ_DATA"] = _PROJ_DIR
    os.environ.pop("PROJ_LIB", None)

import argparse
import csv
import io
import math
import re
import time
import zipfile
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import escape

import numpy as np
import requests
from pyproj import Geod
from scipy.ndimage import map_coordinates
from tqdm import tqdm

GEOD = Geod(ellps="WGS84")
EARTH_RADIUS_M = 6_371_000.0
FT_PER_M = 1.0 / 0.3048

STATION_URL = "https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt"
DEM_BUCKET = "https://copernicus-dem-90m.s3.amazonaws.com"
REPO_DIR = Path(__file__).resolve().parent

# ROC redundant test RDA and the NSSL research radar, both in Norman, OK and
# co-located with KTLX. Pass them explicitly with --sites to include them.
NON_OPERATIONAL = {"KCRI", "KOUN"}


# --------------------------------------------------------------------------- #
# Radar sites
# --------------------------------------------------------------------------- #
@dataclass
class Site:
    icao: str
    name: str
    state: str
    lat: float
    lon: float
    antenna_m: float  # antenna (feedhorn) elevation, m MSL
    elev_angle: float = 0.5
    region: str = ""


def region_of(lat: float, lon: float) -> str:
    if lon > 0:
        return "Guam"
    if lat > 50 and lon < -129:
        return "Alaska"
    if lat < 25 and lon < -150:
        return "Hawaii"
    if lat < 20 and lon > -70:
        return "Puerto Rico"
    return "CONUS"


def load_sites(path: Path) -> list[Site]:
    """Parse NCEI's fixed-width nexrad-stations.txt (ELEV = antenna elevation, ft MSL)."""
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        resp = requests.get(STATION_URL, timeout=60)
        resp.raise_for_status()
        path.write_bytes(resp.content)

    lines = path.read_text(encoding="latin-1").splitlines()
    spans = [m.span() for m in re.finditer(r"-+", lines[1])]
    names = [lines[0][a:b].strip() for a, b in spans]

    sites = []
    for line in lines[2:]:
        if not line.strip():
            continue
        rec = {n: line[a:b].strip() for n, (a, b) in zip(names, spans)}
        if "NEXRAD" not in rec["STNTYPE"]:
            continue
        if rec["COUNTRY"] not in ("UNITED STATES", "GUAM"):
            continue
        lat, lon = float(rec["LAT"]), float(rec["LON"])
        sites.append(Site(
            icao=rec["ICAO"],
            name=rec["NAME"].title(),
            state=rec["ST"],
            lat=lat,
            lon=lon,
            antenna_m=float(rec["ELEV"]) / FT_PER_M,
            region=region_of(lat, lon),
        ))
    return sites


# --------------------------------------------------------------------------- #
# Beam geometry
# --------------------------------------------------------------------------- #
def beam_height(r, elev_deg, antenna_m, ae):
    """Beam-center height (m MSL) and ground distance (m) at slant range r (m)."""
    th = np.deg2rad(elev_deg)
    h = np.sqrt(r**2 + ae**2 + 2.0 * r * ae * np.sin(th)) - ae + antenna_m
    s = ae * np.arcsin(r * np.cos(th) / (ae + h - antenna_m))
    return h, s


def partial_blockage(y, a):
    """Fraction of a circular beam of radius a blocked by terrain y above beam center."""
    yc = np.clip(y, -a, a)
    return (yc * np.sqrt(a * a - yc * yc) + a * a * np.arcsin(yc / a)
            + 0.5 * np.pi * a * a) / (np.pi * a * a)


# --------------------------------------------------------------------------- #
# Copernicus GLO-90 DEM
# --------------------------------------------------------------------------- #
def tile_name(lat0: int, lon0: int) -> str:
    ns = "N" if lat0 >= 0 else "S"
    ew = "E" if lon0 >= 0 else "W"
    return f"Copernicus_DSM_COG_30_{ns}{abs(lat0):02d}_00_{ew}{abs(lon0):03d}_00_DEM"


def glo90_cols_per_degree(lat: float) -> int:
    """GLO-90 tiles thin out in longitude toward the poles."""
    a = abs(lat)
    for limit, cols in ((50, 1200), (60, 800), (70, 600), (80, 400), (85, 240)):
        if a < limit:
            return cols
    return 120


def tiles_for_bbox(bbox) -> list[tuple[int, int]]:
    """1-degree tiles whose pixel centers cover bbox. Tile (lat0, lon0) holds
    pixel centers lat in (lat0, lat0+1], lon in [lon0, lon0+1)."""
    west, south, east, north = bbox
    lats = range(math.ceil(south) - 1, math.ceil(north))
    lons = range(math.floor(west), math.floor(east) + 1)
    return [(la, lo) for la in lats for lo in lons]


def load_tile_list(cache: Path) -> set[str]:
    path = cache / "tileList.txt"
    if not path.exists():
        cache.mkdir(parents=True, exist_ok=True)
        resp = requests.get(f"{DEM_BUCKET}/tileList.txt", timeout=60)
        resp.raise_for_status()
        path.write_bytes(resp.content)
    return set(path.read_text().split())


def _download(name: str, cache: Path) -> None:
    dest = cache / f"{name}.tif"
    if dest.exists():
        return
    tmp = dest.with_suffix(".part")
    for attempt in range(5):
        try:
            with requests.get(f"{DEM_BUCKET}/{name}/{name}.tif", stream=True, timeout=120) as resp:
                resp.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in resp.iter_content(1 << 20):
                        fh.write(chunk)
            tmp.rename(dest)
            return
        except requests.RequestException:
            if attempt == 4:
                raise
            time.sleep(2 ** attempt)


def ensure_tiles(names: set[str], cache: Path, workers: int) -> None:
    missing = sorted(n for n in names if not (cache / f"{n}.tif").exists())
    if not missing:
        return
    with ThreadPoolExecutor(workers) as ex:
        list(tqdm(ex.map(lambda n: _download(n, cache), missing), total=len(missing),
                  desc="Downloading GLO-90 tiles", unit="tile"))


@dataclass
class DemMosaic:
    data: np.ndarray  # float32, rows north->south
    lat_top: float    # latitude of first row's pixel centers
    lon_left: float   # longitude of first column's pixel centers
    dlat: float
    dlon: float

    def sample(self, lat, lon):
        """Bilinear terrain height (m) at arbitrary lat/lon arrays."""
        rows = (self.lat_top - lat) / self.dlat
        cols = (lon - self.lon_left) / self.dlon
        coords = np.vstack([rows.ravel(), cols.ravel()])
        return map_coordinates(self.data, coords, order=1, mode="nearest").reshape(np.shape(lat))


def build_mosaic(bbox, cache: Path, available: set[str], ref_lat: float) -> DemMosaic:
    """Mosaic GLO-90 tiles over bbox on the native grid at ref_lat. Ocean tiles
    (absent from the bucket) are left at 0 m."""
    import rasterio

    rows_per_deg = 1200
    cols_per_deg = glo90_cols_per_degree(ref_lat)
    west, south, east, north = bbox
    lat_top = math.ceil(north * rows_per_deg) / rows_per_deg
    lat_bot = math.floor(south * rows_per_deg) / rows_per_deg
    lon_left = math.floor(west * cols_per_deg) / cols_per_deg
    lon_right = math.ceil(east * cols_per_deg) / cols_per_deg
    nrows = round((lat_top - lat_bot) * rows_per_deg) + 1
    ncols = round((lon_right - lon_left) * cols_per_deg) + 1
    data = np.zeros((nrows, ncols), dtype=np.float32)

    for lat0, lon0 in tiles_for_bbox(bbox):
        name = tile_name(lat0, lon0)
        if name not in available:
            continue
        with rasterio.open(cache / f"{name}.tif") as src:
            tile = src.read(1)
        tile_cols = tile.shape[1]
        if tile_cols != cols_per_deg:
            # Resample columns onto this mosaic's longitude spacing.
            tgt_cols = np.arange(cols_per_deg) * tile_cols / cols_per_deg
            rr, cc = np.meshgrid(np.arange(tile.shape[0]), tgt_cols, indexing="ij")
            tile = map_coordinates(tile, [rr, cc], order=1, mode="nearest").astype(np.float32)
        r0 = round((lat_top - (lat0 + 1)) * rows_per_deg)
        c0 = round((lon0 - lon_left) * cols_per_deg)
        _paste(data, tile, r0, c0)

    return DemMosaic(data, lat_top, lon_left, 1.0 / rows_per_deg, 1.0 / cols_per_deg)


def _paste(dst: np.ndarray, src: np.ndarray, r0: int, c0: int) -> None:
    r1, c1 = r0 + src.shape[0], c0 + src.shape[1]
    dr0, dc0 = max(r0, 0), max(c0, 0)
    dr1, dc1 = min(r1, dst.shape[0]), min(c1, dst.shape[1])
    if dr0 >= dr1 or dc0 >= dc1:
        return
    dst[dr0:dr1, dc0:dc1] = src[dr0 - r0:dr1 - r0, dc0 - c0:dc1 - c0]


# --------------------------------------------------------------------------- #
# Per-radar computation
# --------------------------------------------------------------------------- #
@dataclass
class Config:
    max_range_km: float = 230.0
    gate_m: float = 250.0
    az_step: float = 0.5
    beamwidth: float = 0.95
    k_factor: float = 4.0 / 3.0
    block_threshold: float = 0.5
    min_block_range_km: float = 1.5
    res_deg: float = 0.01
    cache: Path = Path("dem_cache")

    @property
    def ae(self) -> float:
        return self.k_factor * EARTH_RADIUS_M

    def ground_range(self, site: Site) -> float:
        return float(beam_height(self.max_range_km * 1e3, site.elev_angle, site.antenna_m, self.ae)[1])


def coverage_bbox(site: Site, ground_m: float, pad_deg: float = 0.02):
    az = np.arange(0, 360, 1.0)
    lon, lat, _ = GEOD.fwd(np.full_like(az, site.lon), np.full_like(az, site.lat), az,
                           np.full_like(az, ground_m))
    return (lon.min() - pad_deg, lat.min() - pad_deg, lon.max() + pad_deg, lat.max() + pad_deg)


@dataclass
class Grid:
    """Regular lat/lon grid addressed by pixel centers."""
    lon_left: float
    lat_top: float
    res: float
    nrows: int
    ncols: int

    @classmethod
    def covering(cls, bbox, res):
        west, south, east, north = bbox
        lon_left = math.floor(west / res) * res
        lat_top = math.ceil(north / res) * res
        return cls(lon_left, lat_top, res,
                   round((lat_top - math.floor(south / res) * res) / res) + 1,
                   round((math.ceil(east / res) * res - lon_left) / res) + 1)

    def offset_in(self, parent: "Grid") -> tuple[int, int]:
        return (round((parent.lat_top - self.lat_top) / self.res),
                round((self.lon_left - parent.lon_left) / self.res))

    def centers(self):
        lats = self.lat_top - np.arange(self.nrows) * self.res
        lons = self.lon_left + np.arange(self.ncols) * self.res
        return np.meshgrid(lons, lats)

    @property
    def bounds(self):  # west, south, east, north outer pixel edges
        h = self.res / 2
        return (self.lon_left - h, self.lat_top - (self.nrows - 1) * self.res - h,
                self.lon_left + (self.ncols - 1) * self.res + h, self.lat_top + h)


@dataclass
class RadarResult:
    site: Site
    grid: Grid
    ground_range_m: float
    center_agl: np.ndarray  # m, NaN outside coverage
    bottom_agl: np.ndarray  # m
    cbb: np.ndarray         # cumulative beam blockage fraction
    stats: dict = field(default_factory=dict)


def process_radar(site: Site, cfg: Config, available: set[str]) -> RadarResult:
    ae = cfg.ae
    half_bw = np.deg2rad(cfg.beamwidth / 2.0)

    # --- polar grid: trace beam & blockage along every radial ---
    r = (np.arange(int(cfg.max_range_km * 1e3 / cfg.gate_m)) + 0.5) * cfg.gate_m
    h_msl, s = beam_height(r, site.elev_angle, site.antenna_m, ae)
    s_max = float(s[-1])
    bbox = coverage_bbox(site, s_max)
    dem = build_mosaic(bbox, cfg.cache, available, site.lat)

    az = (np.arange(int(round(360 / cfg.az_step))) + 0.5) * cfg.az_step
    AZ, S = np.meshgrid(az, s, indexing="ij")
    lon_p, lat_p, _ = GEOD.fwd(np.full_like(AZ, site.lon), np.full_like(AZ, site.lat), AZ, S)
    terrain = dem.sample(lat_p, lon_p)

    radius = r * np.tan(half_bw)  # half-power beam radius (m)
    pbb = partial_blockage(terrain - h_msl[None, :], radius[None, :])
    # Inside the antenna's near field (~2D^2/lambda = 1.4 km) the beam isn't
    # formed, and the DSM there often contains the radome/tower itself.
    pbb[:, r < cfg.min_block_range_km * 1e3] = 0.0
    cbb = np.maximum.accumulate(pbb, axis=1)

    # --- resample onto the lat/lon output grid ---
    grid = Grid.covering(bbox, cfg.res_deg)
    lon_g, lat_g = grid.centers()
    az_g, _, dist_g = GEOD.inv(np.full_like(lon_g, site.lon), np.full_like(lat_g, site.lat),
                               lon_g, lat_g)
    inside = dist_g <= s_max
    r_g = np.interp(dist_g, s, r)
    h_g, _ = beam_height(r_g, site.elev_angle, site.antenna_m, ae)
    z_g = dem.sample(lat_g, lon_g)

    center = (h_g - z_g).astype(np.float32)
    bottom = (center - r_g * np.tan(half_bw)).astype(np.float32)
    ia = (np.floor(np.mod(az_g, 360.0) / cfg.az_step).astype(int)) % len(az)
    ir = np.clip(np.floor(r_g / cfg.gate_m).astype(int), 0, len(r) - 1)
    cbb_g = cbb[ia, ir].astype(np.float32)
    # The polar gate and the grid pixel see slightly different terrain; a beam
    # center below this pixel's ground is blocked here regardless.
    cbb_g[center <= 0] = 1.0
    for arr in (center, bottom, cbb_g):
        arr[~inside] = np.nan

    # --- area-weighted stats over the coverage disk ---
    w = np.cos(np.deg2rad(lat_g))[inside]
    blocked = cbb_g[inside] >= cfg.block_threshold
    clear_center = center[inside]
    stats = {
        "blocked_pct": 100.0 * w[blocked].sum() / w.sum(),
        **{f"below_{t}ft_pct": 100.0 * w[~blocked & (clear_center < t / FT_PER_M)].sum() / w.sum()
           for t in (1000, 3000, 6000, 10000)},
    }
    return RadarResult(site, grid, s_max, center, bottom, cbb_g, stats)


def _worker(args):
    return process_radar(*args)


# --------------------------------------------------------------------------- #
# Compositing
# --------------------------------------------------------------------------- #
@dataclass
class Composite:
    grid: Grid
    center: np.ndarray
    bottom: np.ndarray
    owner: np.ndarray     # index into site list of the radar providing the lowest beam
    in_range: np.ndarray

    @classmethod
    def empty(cls, grid: Grid):
        shape = (grid.nrows, grid.ncols)
        return cls(grid, np.full(shape, np.inf, np.float32), np.full(shape, np.nan, np.float32),
                   np.full(shape, -1, np.int16), np.zeros(shape, bool))

    def add(self, res: RadarResult, idx: int, threshold: float):
        r0, c0 = res.grid.offset_in(self.grid)
        sl = (slice(r0, r0 + res.grid.nrows), slice(c0, c0 + res.grid.ncols))
        valid = np.isfinite(res.center_agl) & (res.cbb < threshold)
        better = valid & (res.center_agl < self.center[sl])
        self.center[sl][better] = res.center_agl[better]
        self.bottom[sl][better] = res.bottom_agl[better]
        self.owner[sl][better] = idx
        self.in_range[sl] |= np.isfinite(res.center_agl)

    def finish(self):
        self.center[~np.isfinite(self.center)] = np.nan

    @property
    def blocked(self):
        return self.in_range & np.isnan(self.center)


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #
LEVELS = {
    "ft": [0, 1000, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000, 15000],
    "m": [0, 250, 500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000],
}
BLOCKED_RGB = (90, 90, 90)


def level_colors(n: int):
    """Lightness-ordered ramp: bright = low beam (good coverage), dark = high beam."""
    from matplotlib import colormaps
    cmap = colormaps["viridis"]
    return [tuple(int(255 * c) for c in cmap(1.0 - i / (n - 1))[:3]) for i in range(n)]


def classify(values_m: np.ndarray, units: str) -> np.ndarray:
    """Palette index per pixel: 0 = transparent, 1..n = height bins."""
    vals = values_m * (FT_PER_M if units == "ft" else 1.0)
    idx = np.digitize(vals, LEVELS[units][1:]) + 1
    idx[np.isnan(vals)] = 0
    return idx.astype(np.uint8)


def png_bytes(index: np.ndarray, palette: list[tuple[int, int, int]]) -> bytes:
    from PIL import Image
    img = Image.fromarray(index, mode="P")
    flat = [0, 0, 0] + [c for rgb in palette for c in rgb]
    img.putpalette(flat + [0] * (768 - len(flat)))
    buf = io.BytesIO()
    img.save(buf, format="PNG", transparency=0, optimize=True)
    return buf.getvalue()


def legend_png(units: str) -> bytes:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    levels = LEVELS[units]
    colors = level_colors(len(levels))
    fmt = (lambda v: f"{v / 1000:g}k") if units == "ft" else (lambda v: f"{v:g}")
    labels = [f"{fmt(a)}–{fmt(b)}" for a, b in zip(levels[:-1], levels[1:])] + [f"≥ {fmt(levels[-1])}"]

    fig, ax = plt.subplots(figsize=(2.3, 0.28 * (len(labels) + 2) + 0.5), dpi=150)
    fig.patch.set_facecolor((1, 1, 1, 0.88))
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, len(labels) + 2)
    ax.text(0.04, len(labels) + 1.4, f"Lowest beam height AGL ({units})",
            fontsize=8.5, fontweight="bold", va="center", color="#1a1a1a")
    rows = list(zip(labels, colors)) + [("Terrain-blocked", BLOCKED_RGB)]
    for i, (lab, rgb) in enumerate(rows):
        y = len(labels) - i
        ax.add_patch(Rectangle((0.04, y - 0.35), 0.18, 0.7, color=[c / 255 for c in rgb], lw=0))
        ax.text(0.28, y, lab, fontsize=8, va="center", color="#1a1a1a")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", facecolor=fig.get_facecolor(), bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return buf.getvalue()


def chunk_overlays(index: np.ndarray, grid: Grid, palette, prefix: str, max_px: int = 2048):
    """Split a big image into Google-Earth-friendly tiles: [(href, png, (w,s,e,n))]."""
    out = []
    for r0 in range(0, grid.nrows, max_px):
        for c0 in range(0, grid.ncols, max_px):
            block = index[r0:r0 + max_px, c0:c0 + max_px]
            if not block.any():
                continue
            sub = Grid(grid.lon_left + c0 * grid.res, grid.lat_top - r0 * grid.res, grid.res,
                       block.shape[0], block.shape[1])
            out.append((f"{prefix}_{r0}_{c0}.png", png_bytes(np.ascontiguousarray(block), palette),
                        sub.bounds))
    return out


# --------------------------------------------------------------------------- #
# KML
# --------------------------------------------------------------------------- #
def kml_ground_overlay(name, href, bounds, visible=True, opacity=0.8):
    w, s, e, n = bounds
    alpha = f"{int(round(opacity * 255)):02x}"
    return (f"<GroundOverlay><name>{escape(name)}</name><visibility>{int(visible)}</visibility>"
            f"<color>{alpha}ffffff</color><Icon><href>{href}</href></Icon>"
            f"<LatLonBox><north>{n:.6f}</north><south>{s:.6f}</south>"
            f"<east>{e:.6f}</east><west>{w:.6f}</west></LatLonBox></GroundOverlay>")


def kml_folder(name, children, visible=True, open_=False, description=""):
    desc = f"<description>{escape(description)}</description>" if description else ""
    return (f"<Folder><name>{escape(name)}</name><visibility>{int(visible)}</visibility>"
            f"<open>{int(open_)}</open>{desc}{''.join(children)}</Folder>")


def range_ring(site: Site, ground_m: float) -> str:
    az = np.linspace(0, 360, 181)
    lon, lat, _ = GEOD.fwd(np.full_like(az, site.lon), np.full_like(az, site.lat), az,
                           np.full_like(az, ground_m))
    coords = " ".join(f"{x:.5f},{y:.5f},0" for x, y in zip(lon, lat))
    return (f"<Placemark><name>{site.icao} range ring</name><visibility>0</visibility>"
            f"<styleUrl>#ring</styleUrl><LineString><tessellate>1</tessellate>"
            f"<coordinates>{coords}</coordinates></LineString></Placemark>")


def site_placemark(res: RadarResult, cfg: Config, units: str) -> str:
    s = res.site
    st = res.stats
    to_u = FT_PER_M if units == "ft" else 1.0
    rows = [
        ("Name", f"{s.name}, {s.state}"),
        ("Antenna elevation", f"{s.antenna_m * FT_PER_M:,.0f} ft / {s.antenna_m:,.0f} m MSL"),
        ("Elevation angle", f"{s.elev_angle:g}°"),
        ("Max slant range", f"{cfg.max_range_km:g} km"),
        ("Beam height at max range", f"{(beam_height(cfg.max_range_km * 1e3, s.elev_angle, 0, cfg.ae)[0]) * to_u:,.0f} {units} above antenna"),
        ("Coverage area terrain-blocked", f"{st['blocked_pct']:.1f}%"),
        ("Unblocked & beam below 1k ft AGL", f"{st['below_1000ft_pct']:.1f}%"),
        ("Unblocked & beam below 3k ft AGL", f"{st['below_3000ft_pct']:.1f}%"),
        ("Unblocked & beam below 6k ft AGL", f"{st['below_6000ft_pct']:.1f}%"),
        ("Unblocked & beam below 10k ft AGL", f"{st['below_10000ft_pct']:.1f}%"),
    ]
    table = "".join(f"<tr><td><b>{k}</b></td><td>{v}</td></tr>" for k, v in rows)
    return (f"<Placemark><name>{s.icao}</name><styleUrl>#radar</styleUrl>"
            f"<description><![CDATA[<table>{table}</table>]]></description>"
            f"<Point><coordinates>{s.lon:.6f},{s.lat:.6f},0</coordinates></Point></Placemark>")


def write_google_earth(path: Path, results: list[RadarResult], composites: dict[str, Composite],
                       cfg: Config, units: str, individual: bool, image_base_url: str | None = None):
    """Write a .kmz (images zipped inside) or a .kml with images in files/ beside it.
    image_base_url makes the KML's image links absolute, e.g. to a web host;
    Google Earth Web can only load KML overlays whose images are online."""
    levels = LEVELS[units]
    palette = level_colors(len(levels))
    files: dict[str, bytes] = {"files/legend.png": legend_png(units)}
    prefix = image_base_url.rstrip("/") + "/" if image_base_url else ""

    def overlays(name, index, grid, pal, prefix, visible=True, opacity=0.8):
        items = []
        for href, png, bounds in chunk_overlays(index, grid, pal, prefix):
            files[f"files/{href}"] = png
            items.append(kml_ground_overlay(name, f"{prefix}files/{href}", bounds, visible, opacity))
        return items

    center_folders, bottom_folders, blocked_folders = [], [], []
    for region, comp in composites.items():
        slug = region.lower().replace(" ", "_")
        center_folders.append(kml_folder(region, overlays(
            f"{region} beam center", classify(comp.center, units), comp.grid, palette, f"center_{slug}")))
        bottom_folders.append(kml_folder(region, overlays(
            f"{region} beam bottom", classify(comp.bottom, units), comp.grid, palette, f"bottom_{slug}",
            visible=False), visible=False))
        blocked_folders.append(kml_folder(region, overlays(
            f"{region} blocked", comp.blocked.astype(np.uint8), comp.grid, [BLOCKED_RGB], f"blocked_{slug}",
            opacity=0.6)))

    per_radar = []
    if individual:
        for res in tqdm(results, desc="Rendering per-radar overlays", unit="radar"):
            idx = classify(res.center_agl, units)
            idx[res.cbb >= cfg.block_threshold] = len(palette) + 1
            items = overlays(f"{res.site.icao} beam center", idx, res.grid, palette + [BLOCKED_RGB],
                             f"radar_{res.site.icao}", visible=False)
            items.append(range_ring(res.site, res.ground_range_m))
            per_radar.append(kml_folder(f"{res.site.icao} – {res.site.name}, {res.site.state}",
                                        items, visible=False))

    placemarks = [site_placemark(r, cfg, units) for r in results]
    angles = sorted({r.site.elev_angle for r in results})
    about = (f"Lowest-tilt ({'/'.join(f'{a:g}' for a in angles)}°) beam-center height above the "
             f"Copernicus GLO-90 surface, {cfg.k_factor:.3g}·Re refraction, {cfg.beamwidth}° beam, "
             f"{cfg.max_range_km:g} km max range. Composite shows the lowest beam among radars with "
             f"< {cfg.block_threshold:.0%} cumulative terrain blockage. "
             f"Terrain: Copernicus DEM GLO-90 © DLR e.V. 2010-2014 and © Airbus Defence and Space "
             f"GmbH 2014-2018 provided under COPERNICUS by the European Union and ESA; all rights reserved.")

    kml = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2"><Document>
<name>NEXRAD lowest-beam height AGL</name><open>1</open>
<description>{escape(about)}</description>
<Style id="radar"><IconStyle><scale>0.8</scale><Icon><href>https://maps.google.com/mapfiles/kml/shapes/target.png</href></Icon></IconStyle><LabelStyle><scale>0.7</scale></LabelStyle></Style>
<Style id="ring"><LineStyle><color>ccffffff</color><width>1.5</width></LineStyle></Style>
<ScreenOverlay><name>Legend</name><Icon><href>{prefix}files/legend.png</href></Icon>
<overlayXY x="0" y="0" xunits="fraction" yunits="fraction"/><screenXY x="10" y="30" xunits="pixels" yunits="pixels"/>
<size x="0" y="0" xunits="pixels" yunits="pixels"/></ScreenOverlay>
{kml_folder("Radar sites", placemarks)}
{kml_folder("Composite: lowest beam-center height AGL", center_folders, open_=True)}
{kml_folder("Composite: lowest beam-bottom height AGL (half-power)", bottom_folders, visible=False)}
{kml_folder(f"Terrain-blocked (≥{cfg.block_threshold:.0%} beam blockage, every radar)", blocked_folders)}
{kml_folder("Individual radars", per_radar, visible=False) if per_radar else ""}
</Document></kml>
"""
    if path.suffix == ".kmz":
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as zf:
            zf.writestr("doc.kml", kml)
            for name, data in files.items():
                zf.writestr(name, data, compress_type=zipfile.ZIP_STORED)
    else:
        path.write_text(kml, encoding="utf-8")
        for name, data in files.items():
            (path.parent / name).parent.mkdir(parents=True, exist_ok=True)
            (path.parent / name).write_bytes(data)


# --------------------------------------------------------------------------- #
# GeoTIFF / CSV
# --------------------------------------------------------------------------- #
def write_geotiff(path: Path, comp: Composite):
    import rasterio
    from rasterio.transform import from_origin

    w, s, e, n = comp.grid.bounds
    profile = dict(driver="GTiff", width=comp.grid.ncols, height=comp.grid.nrows, count=3,
                   dtype="float32", crs="EPSG:4326", nodata=np.nan, compress="deflate",
                   predictor=3, tiled=True,
                   transform=from_origin(w, n, comp.grid.res, comp.grid.res))
    owner = comp.owner.astype(np.float32)
    owner[comp.owner < 0] = np.nan
    with rasterio.open(path, "w", **profile) as dst:
        for band, (arr, desc) in enumerate([
            (comp.center, "lowest beam-center height AGL (m)"),
            (comp.bottom, "beam-bottom height AGL (m), same radar"),
            (owner, "radar index (see radar_summary.csv)"),
        ], start=1):
            dst.write(arr, band)
            dst.set_band_description(band, desc)


def write_summary(path: Path, results: list[RadarResult]):
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["index", "icao", "name", "state", "region", "lat", "lon", "antenna_m_msl",
                     "elev_angle_deg", "blocked_pct", "below_1000ft_pct", "below_3000ft_pct",
                     "below_6000ft_pct", "below_10000ft_pct"])
        for i, r in enumerate(results):
            s = r.site
            wr.writerow([i, s.icao, s.name, s.state, s.region, s.lat, s.lon, f"{s.antenna_m:.1f}",
                         s.elev_angle, *(f"{r.stats[k]:.2f}" for k in (
                             "blocked_pct", "below_1000ft_pct", "below_3000ft_pct",
                             "below_6000ft_pct", "below_10000ft_pct"))])


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--sites", help="comma-separated ICAO list (default: every U.S. NEXRAD)")
    p.add_argument("--stations-file", type=Path, default=REPO_DIR / "data" / "nexrad-stations.txt",
                   help="NCEI nexrad-stations.txt (downloaded if missing)")
    p.add_argument("--elev-angle", type=float, default=0.5, help="lowest tilt, degrees (default 0.5)")
    p.add_argument("--elev-override", action="append", default=[], metavar="ICAO=DEG",
                   help="per-site lowest tilt, e.g. --elev-override KMSX=0.2 (repeatable)")
    p.add_argument("--max-range-km", type=float, default=230.0, help="max slant range (default 230)")
    p.add_argument("--gate-m", type=float, default=250.0, help="polar range spacing (default 250 m)")
    p.add_argument("--az-step", type=float, default=0.5, help="polar azimuth spacing (default 0.5°)")
    p.add_argument("--beamwidth", type=float, default=0.95, help="half-power beamwidth (default 0.95°)")
    p.add_argument("--k-factor", type=float, default=4.0 / 3.0, help="effective earth radius factor")
    p.add_argument("--block-threshold", type=float, default=0.5,
                   help="cumulative blockage fraction treated as blocked (default 0.5)")
    p.add_argument("--min-block-range-km", type=float, default=1.5,
                   help="ignore terrain closer than this (antenna near field; default 1.5 km)")
    p.add_argument("--res-deg", type=float, default=0.01, help="output grid spacing (default 0.01°)")
    p.add_argument("--units", choices=["ft", "m"], default="ft", help="KMZ color-scale units")
    p.add_argument("--no-individual", action="store_true", help="skip per-radar overlays")
    p.add_argument("--format", choices=["kmz", "kml", "both"], default="kmz",
                   help="kmz = one self-contained file (Google Earth Pro); kml = .kml plus files/ "
                        "folder of images (default kmz)")
    p.add_argument("--image-base-url", help="URL where the files/ folder will be hosted; makes the "
                   ".kml's image links absolute so it loads in Google Earth Web")
    p.add_argument("--cache-dir", type=Path, default=REPO_DIR / "dem_cache")
    p.add_argument("--out-dir", type=Path, default=REPO_DIR / "output")
    p.add_argument("--workers", type=int, default=max(1, min(4, (os.cpu_count() or 2) - 1)))
    p.add_argument("--download-workers", type=int, default=8)
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    cfg = Config(max_range_km=args.max_range_km, gate_m=args.gate_m, az_step=args.az_step,
                 beamwidth=args.beamwidth, k_factor=args.k_factor,
                 block_threshold=args.block_threshold, min_block_range_km=args.min_block_range_km,
                 res_deg=args.res_deg, cache=args.cache_dir)

    sites = load_sites(args.stations_file)
    if args.sites:
        wanted = [s.strip().upper() for s in args.sites.split(",") if s.strip()]
        unknown = set(wanted) - {s.icao for s in sites}
        if unknown:
            sys.exit(f"Unknown site(s): {', '.join(sorted(unknown))}")
        sites = [s for s in sites if s.icao in wanted]
    else:
        sites = [s for s in sites if s.icao not in NON_OPERATIONAL]

    overrides = {}
    for item in args.elev_override:
        icao, _, deg = item.partition("=")
        overrides[icao.strip().upper()] = float(deg)
    for s in sites:
        s.elev_angle = overrides.get(s.icao, args.elev_angle)
    sites.sort(key=lambda s: s.icao)
    print(f"{len(sites)} radars, {cfg.max_range_km:g} km range, lowest tilt {args.elev_angle:g}° "
          f"({len(overrides)} override(s))")

    # Download every DEM tile up front so workers only read from disk.
    available = load_tile_list(cfg.cache)
    needed = set()
    for s in sites:
        bbox = coverage_bbox(s, cfg.ground_range(s))
        needed |= {tile_name(*t) for t in tiles_for_bbox(bbox)} & available
    ensure_tiles(needed, cfg.cache, args.download_workers)

    # Region composites sized to cover every radar in the region.
    regions: dict[str, list] = {}
    for s in sites:
        regions.setdefault(s.region, []).append(coverage_bbox(s, cfg.ground_range(s)))
    composites = {}
    for region, boxes in regions.items():
        b = np.array(boxes)
        composites[region] = Composite.empty(Grid.covering(
            (b[:, 0].min(), b[:, 1].min(), b[:, 2].max(), b[:, 3].max()), cfg.res_deg))

    results = []
    jobs = [(s, cfg, available) for s in sites]
    with ProcessPoolExecutor(args.workers) as ex:
        for res in tqdm(ex.map(_worker, jobs), total=len(jobs), desc="Tracing beams", unit="radar"):
            composites[res.site.region].add(res, len(results), cfg.block_threshold)
            results.append(res)
    for comp in composites.values():
        comp.finish()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    for region, comp in composites.items():
        write_geotiff(args.out_dir / f"{region.lower().replace(' ', '_')}_beam_height_agl.tif", comp)
    write_summary(args.out_dir / "radar_summary.csv", results)
    for ext in (["kmz", "kml"] if args.format == "both" else [args.format]):
        out = args.out_dir / f"nexrad_beam_height_agl.{ext}"
        write_google_earth(out, results, composites, cfg, args.units, not args.no_individual,
                           args.image_base_url if ext == "kml" else None)
        print(f"Wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
    print(f"GeoTIFFs and radar_summary.csv in {args.out_dir}")


if __name__ == "__main__":
    main()
