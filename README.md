# NEXRAD beam height AGL

Computes the lowest-tilt beam height above ground level (AGL) for every U.S.
NEXRAD (WSR-88D) radar against the Copernicus GLO-90 (90 m) digital surface
model. It also flags areas where terrain blocks the beam. Results are written
as a Google Earth KMZ, GeoTIFFs and a per-radar CSV.

## Quick start

```bash
conda env create -f environment.yml
conda activate nexrad_elevations

python nexrad_beam_height.py                     # all 157 operational U.S. radars
python nexrad_beam_height.py --sites KABX,KFSX   # just a few
```

The first full run downloads about 1,370 GLO-90 tiles (~6 GB) from the public
AWS bucket `s3://copernicus-dem-90m`. They are cached in `dem_cache/`, so later
runs are fast.

Open `output/nexrad_beam_height_agl.kmz` in **Google Earth Pro**.

### Google Earth Web

Earth Web can't read images packed in a KMZ or stored on your computer, so use
the KML published on this repo's GitHub Pages site. Its overlay images are
hosted there:

1. Download [nexrad_beam_height_agl.kml](https://ethaniebbert.github.io/nexrad_elevations/nexrad_beam_height_agl.kml).
2. In Earth Web: *Projects → New project → Import KML file from computer*.

To rebuild that KML after changing settings, point its image links at your
host and publish `output/nexrad_beam_height_agl.kml` with `output/files/` to
that URL (here, the `gh-pages` branch):

```bash
python nexrad_beam_height.py --format kml \
    --image-base-url https://ethaniebbert.github.io/nexrad_elevations/
```

## What's in the KMZ

| Layer | Default | Description |
|---|---|---|
| Radar sites | on | Placemarks. Click one for antenna elevation, % of coverage blocked, and % of area with the beam below 1k/3k/6k/10k ft AGL |
| Composite: lowest beam-center height AGL | on | Lowest unblocked beam among all radars at each point |
| Composite: lowest beam-bottom height AGL | off | Same radar, bottom of the half-power beam |
| Terrain-blocked | on | In range of at least one radar, but every radar is ≥ 50% blocked |
| Individual radars | off | Each radar's own beam-height map (blocked areas in gray) and range ring |

Use the transparency slider in Google Earth's Places panel to see terrain
underneath.

## Other outputs (`output/`)

- `<region>_beam_height_agl.tif`: GeoTIFFs (EPSG:4326, 0.01°) for CONUS,
  Alaska, Hawaii, Puerto Rico and Guam. Band 1 is the lowest beam-center height
  AGL (m), band 2 is the beam-bottom height AGL (m), and band 3 is the index of
  the radar providing it.
- `radar_summary.csv`: per-radar statistics. Its `index` column matches band 3.

## Method

1. **Sites**: NCEI's [`nexrad-stations.txt`](https://www.ncei.noaa.gov/access/homr/file/nexrad-stations.txt)
   (snapshot in `data/`). Its `ELEV` column is the *antenna* elevation in ft MSL
   and matches the NWS WSR-88D list, so no tower height is added. KCRI and KOUN
   (non-operational test and research radars in Norman, OK) are skipped unless
   named with `--sites`.
2. **Beam path**: effective-earth-radius model with k = 4/3
   (Doviak & Zrnić 1993, eq. 2.28). The polar grid has 0.5° azimuth and 250 m
   gates out to 230 km slant range.
3. **Blockage**: terrain is sampled at every gate. Partial beam blockage of a
   circular 0.95° half-power beam follows Bech et al. (2003). Cumulative
   blockage is the running maximum along each radial, and ≥ 50% counts as
   blocked. Terrain within 1.5 km is ignored: that is inside the antenna's near
   field, and the DSM there often contains the radome and tower themselves.
4. **Gridding**: each 0.01° output pixel gets the beam height at its ground
   range, minus the DEM at that pixel. Blockage comes from the matching polar
   gate. The composite keeps the lowest unblocked beam.

### Caveats

- GLO-90 is a **surface** model: it includes tree canopy and buildings. "AGL"
  here means above that surface, which also makes blockage slightly
  conservative.
- Every radar uses a 0.5° lowest tilt by default. Some sites run lower tilts.
  Set those with `--elev-override ICAO=DEG`, e.g. `--elev-override KMSX=0.2`.
- The calculation assumes standard refraction. Super- or sub-refraction can
  change beam heights at long range by hundreds of meters. Try `--k-factor`.
- NEXRAD reflectivity reaches 460 km. Use `--max-range-km 460` for the full
  extent; this needs more tiles and memory.

Run `python nexrad_beam_height.py --help` for every option, including units
(`--units m`), grid spacing, beamwidth and blockage threshold.

## References

- Doviak, R. J., and D. S. Zrnić, 1993: *Doppler Radar and Weather Observations*, 2nd ed. Academic Press.
- Bech, J., et al., 2003: The sensitivity of single polarization weather radar beam blockage correction to variability in the vertical refractivity gradient. *J. Atmos. Oceanic Technol.*, **20**, 845–855.

## Data attribution

Copernicus DEM GLO-90 © DLR e.V. 2010-2014 and © Airbus Defence and Space GmbH
2014-2018 provided under COPERNICUS by the European Union and ESA; all rights
reserved. Radar site metadata: NOAA NCEI.
