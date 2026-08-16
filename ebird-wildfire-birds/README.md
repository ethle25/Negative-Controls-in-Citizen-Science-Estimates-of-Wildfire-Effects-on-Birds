# Wildfire severity and bird occurrence in California

Code for a study joining **eBird** checklists, **MTBS** burn-severity data and
**MODIS** vegetation indices across California, 2015–2026, to ask whether
satellite-derived wildfire severity affects where birds are observed — and
whether citizen-science data can support that kind of causal claim at all.

Nineteen species across four habitat guilds, including **three negative-control
species** (Mallard, American Coot, Black Phoebe) that wildfire cannot plausibly
affect. Those controls are the point of the design: a chaparral fire cannot
reach a duck, so any effect the method reports for them measures the method
rather than the fire.

## This repository contains code only

No data is included. The upstream sources total ~94 GB and the intermediate
tables are 200–600 MB each, so everything data-shaped is `.gitignore`d. What is
here is every script needed to rebuild the analysis from the raw downloads, plus
the four registry CSVs that configure it.

| source | what to download | used for |
|---|---|---|
| **eBird** | Basic Dataset (EBD) **and** Sampling Event Data for `US-CA`, from [ebird.org/data/download](https://ebird.org/data/download) | observations, effort, and the complete-checklist flag that makes zero-filling possible |
| **MTBS** | burn severity mosaics + fire perimeters, [mtbs.gov](https://www.mtbs.gov/) | fire extent, severity class, ignition dates |
| **MODIS** | MOD13A2 v061 granules (h08v04, h08v05, h09v04), NASA LP DAAC | 1 km 16-day NDVI/EVI |

The Sampling Event Data is not optional. A zero can only be written where an
observer certified they reported everything they saw, and that certification
lives in the checklist file — absences cannot be reconstructed from the
sightings alone.

## Layout

```
eBird/            all analysis code, plus the four registries
  figures/        output directory for make_figures.R
docs/             data dictionary and pre-registration
```

`eBird/` is flat because the scripts resolve each other and their outputs
relative to one directory. **The folder name is load-bearing** — scripts compute
their working directory as `$EBIRD_PROJ_ROOT/eBird`.

## Running it

Requires **R ≥ 4.5** and GDAL with the **HDF4** driver (without HDF4 the MODIS
granules will not open).

```r
install.packages(c("data.table", "terra", "sf", "lubridate", "auk",
                   "lightgbm", "xgboost", "gbm",
                   "ggplot2", "patchwork", "ragg"))
```

`torch` is needed only for `deep_models.R`. Point the pipeline at wherever the
data lives:

```bash
export EBIRD_PROJ_ROOT=/path/to/project    # scripts read $EBIRD_PROJ_ROOT/eBird
```

Build, in order:

```bash
Rscript zerofill_from_registry.R    # raw EBD -> zero-filled per species (hours)
Rscript join_species_all_years.R    # ~9 min   -> checklist-grain table
Rscript regrid_cells.R              # ~6 min   -> 5 km x 16-day cell frame
```

Then any analysis. Every script takes its species and burn measure from the
registries, so nothing below has a bird hardcoded in it:

```bash
Rscript dose_burn_effect.R          # dose-response: effect scaled by burned fraction
Rscript event_study.R               # year-by-year effect + pre-trend check
Rscript severity_effect.R           # severity as a dose, separate from extent
Rscript effort_placebo.R            # did fire change how people bird? (it did)
Rscript wild_cluster_boot.R         # the significance verdict
Rscript guild_contrast.R            # between-guild contrasts + placebos
Rscript sign_test.R                 # predicted vs realised direction, all species
Rscript verify_published.R          # 41-check regression test

./run_species.sh                    # whole pipeline, every active species
```

Switch species or treatment with an environment variable:

```bash
EBIRD_PREFIX=olsfly DOSE_MEASURE=two_dose_sev Rscript dose_burn_effect.R
```

Robustness diagnostics write to their own files rather than overwriting a
published result. `guild_contrast.R` takes a list of species to exclude:

```bash
# refit the guild contrasts without the three rarest snag species
GC_DROP=bkbwoo,moublu,lazbun GC_TAG=snagcommon Rscript guild_contrast.R
```

Figures:

```bash
Rscript make_figures.R              # writes figures/fig[0-5]_*.png and .pdf
```

`make_figures.R` parses every number out of the result files rather than
restating it, and each parser is followed by a check that re-derives a value the
report prints for itself — so a refit moves the figures, and a changed report
layout stops the script instead of silently drawing the wrong thing. Figure 0 is
the exception and says so on its face: it is a schematic of the eBird join with
invented example checklists, drawing no estimated quantity.

## Configuration

Four CSVs drive everything; adding a species or a treatment is a row, not a code
change.

| file | role |
|---|---|
| `species_registry.csv` | the 19 species: guild, habitat, movement, and the direction each was **predicted to move before it was run** |
| `treatment_registry.csv` | 25 burn-impact measures across extent, severity, timing, reburn, spatial and fire attributes |
| `feature_registry.csv` | the 53 cell-frame columns by block and role. **Row order is load-bearing** — it is the design-matrix column order, and the tree models sample features off a seeded RNG |
| `species_catalog.csv` | all 773 California eBird species with complete-checklist counts, used to validate names |

## Reproducibility

`verify_published.R` pins 41 published numbers (51 with `VERIFY_FULL=1`) and
should be run **before and after** touching analysis code. A refactor that moves
a result does not crash — it silently returns a different answer.

Three things that repeatedly caused silent, wrong results, documented here so
nobody rediscovers them:

- **Never join across species on `cell_id`.** It is assigned in row-appearance
  order and is not stable between separately built species. Join on `(gx, gy)`.
- **Never name a variable after a `data.table` column.**
- **Cluster by fire, not by cell.** 3,125 burned cells come from 480 fires, and
  the 2021 Dixie Fire alone accounts for 231 of them. Cells inside one perimeter
  share one ignition and one post-fire access regime — one treatment draw, not
  231. Cell-clustered intervals are too narrow.

## What the analysis found

The headline is methodological. The extent-based occurrence design **fails its
own falsification test**: all three water-bird controls register significant
fire effects, larger than the focal species. The mechanism is measurable —
after a cell burns, birders visit **42% fewer distinct locations** within it
(p = 0.0001). Comparing a cell with itself removes anything permanent about that
cell, but it cannot remove a change in *which spots inside it* people visit.
Migration timing fails the same way: sedentary birds show larger apparent shifts
than actual migrants.

Burn history is also a strong predictor for the wrong reason — it adds +0.046
AUC on held-out years, and +0.002 once the model is told the cell's coordinates.
It was largely a location detector, because fires burn chaparral and the focal
species lives in chaparral.

A between-guild severity contrast survived longer than anything else, on the
logic that a confound shared across species cancels in a difference. It is now
qualified: the `control − snag` placebo is **not** silent, and restricting the
snag guild to its three commonest species drops `snag − shrub` from 2.19
(p = 0.006) to 0.58 (p = 0.067). The contrasts decompose exactly, and every one
involving snag reduces to `snag − control`, which does not reach significance
under the restricted guild. Guild membership was pre-registered; which contrasts
to test was chosen afterwards.

The strongest surviving pre-registered result is the **sign test**: 12 of 16
species moved in the direction committed to `species_registry.csv` before the
analysis ran (snag 6/6, forest 4/4, shrub 2/6).

## Citation

eBird Basic Dataset and Sampling Event Data, Cornell Lab of Ornithology —
cite the release you download. MTBS is a joint program of USGS EROS and the
USDA Forest Service. MODIS MOD13A2 v061 is distributed by NASA LP DAAC.
