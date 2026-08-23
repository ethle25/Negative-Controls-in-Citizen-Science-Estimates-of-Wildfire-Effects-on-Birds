# data/

**Empty by design.** No data is committed to this repository — the upstream
sources total roughly 94 GB and every intermediate table rebuilds from the
scripts. See README section 2 for provenance and access.

## What to put here

Unpack the three downloads so that the layout looks like this, with `eBird/`
alongside them:

```
<EBIRD_PROJ_ROOT>/
  eBird/                       the code in this repository
  ebd_US-CA_.../               eBird Basic Dataset + Sampling Event Data
  mtbs_perim/                  MTBS fire perimeters (shapefile)
  mtbs/                        MTBS yearly severity rasters, 1984-2025
  MOD13A2_061/                 MODIS granules, tiles h08v04 h08v05 h09v04
```

Then point the pipeline at it:

```bash
export EBIRD_PROJ_ROOT=/path/to/that/folder
```

## The two eBird files

Download **both** the Basic Dataset (EBD) and the Sampling Event Data. They
arrive as separate files from the same request.

A zero can only be written where an observer certified that they reported
everything they saw, and that certification lives in the Sampling Event file.
Absences cannot be reconstructed from the sightings alone. With only the EBD the
design does not exist.

## Sizes, so you can plan

| | Approximate |
|---|---|
| eBird EBD, California, 2015–2026 | 33 GB unzipped |
| eBird Sampling Event Data | 1.4 GB |
| MTBS rasters and perimeters | ~5 GB |
| MODIS granules, 1,821 files | ~55 GB |
| Intermediate tables the build writes | 200–600 MB each, ~40 GB total |
