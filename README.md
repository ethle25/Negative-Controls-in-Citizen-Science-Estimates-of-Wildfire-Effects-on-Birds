# Negative controls in citizen-science estimates of wildfire effects on birds

Replication package for a study joining **eBird** checklists, **MTBS** burn
severity and **MODIS** vegetation indices across California, 2015–2026.

The study asks whether citizen-science bird records can support a causal claim
about wildfire. It answers no, and shows why — which is why the negative
controls, not the fire estimates, are the centre of the design.

---

## 1. Overview

Nineteen bird species were analysed across four habitat groups. Three of them —
**Mallard, American Coot and Black Phoebe** — are *negative controls*: birds
that live on open water, which a chaparral fire cannot reach. If the method
reports a fire effect on them, the method is measuring something other than
fire.

It does. All three controls register larger effects than the focal species. The
reason is measurable: after a square burns, birders visit **42% fewer distinct
locations** inside it. Comparing a square with itself removes anything permanent
about that square, but it cannot remove a change in *which spots inside it*
people visit.

**Running everything from raw downloads takes roughly 12–15 hours**, most of it
in the zero-filling step. Section 5 gives the order and the timings.

---

## 2. Data availability and provenance

**This repository contains code only.** No data is included: the upstream
sources total roughly 94 GB and the intermediate tables are 200–600 MB each, so
everything data-shaped is `.gitignore`d. Every script needed to rebuild the
analysis from the raw downloads is here, along with the four registry CSVs that
configure it.

All three sources are public and free.

| Source | What to download | Version used | Used for |
|---|---|---|---|
| **eBird** | Basic Dataset (EBD) **and** Sampling Event Data, region `US-CA` — [ebird.org/data/download](https://ebird.org/data/download) | `EBD_relJun-2026`, accessed 2026-07-19 and 2026-07-24 | Observations, effort, and the complete-checklist flag |
| **MTBS** | Burn severity mosaics and fire perimeters — [mtbs.gov](https://www.mtbs.gov/) | 1,956 California perimeters, 1984–2025; accessed 2026-07-20 | Fire extent, severity class, ignition dates |
| **MODIS** | MOD13A2 v061 granules, tiles `h08v04`, `h08v05`, `h09v04` — NASA LP DAAC | 1,821 granules, 2000–2026 | 1 km, 16-day NDVI and EVI |

eBird access requires a free account and agreement to the eBird Data Access
Terms of Use. The download is requested by region and date range and arrives by
email, usually within a day.

> **The Sampling Event Data is not optional.** A zero can only be written where
> an observer certified that they reported everything they saw, and that
> certification lives in the checklist file. Absences cannot be reconstructed
> from the sightings alone — without this file the design does not exist.

---

## 3. Computational requirements

Developed and last run on **Windows 11** with **R 4.5.2**.

**System.** 24 GB RAM is the practical floor. The eBird filtering step streams a
33 GB text file, and the MODIS rasters have exhausted 24 GB when their
coordinates are materialised.

**GDAL with the HDF4 driver.** Without HDF4 the MODIS granules will not open.
Check with `terra::gdal(drivers = TRUE)` and look for `HDF4`.

```r
install.packages(c(
  "data.table", "terra", "sf", "lubridate", "auk",   # data and geospatial
  "lightgbm", "xgboost", "gbm",                      # tree models
  "ggplot2", "patchwork", "ragg"                     # figures
))
```

`torch` is required only for `deep_models.R` and can be skipped otherwise.

---

## 4. Repository structure

```
main.R                 one entry point for the whole replication
DEPENDENCIES.md        package versions the analysis was last run against
CITATION.cff           how to cite this package
LICENSE                MIT for the code; the datasets carry their own terms

eBird/                 all analysis code and the four registries
  figures/             output directory for make_figures.R
data/                  empty — where the three downloads go (see data/README.md)
results/               empty — sessionInfo.txt and anything else generated
docs/
  DATA_DICTIONARY.md       all 78 columns and the 11 join rules
  LIST_LENGTH_PRESPEC.md   pre-registered pass/fail rule for the list-length repair
```

`eBird/` is deliberately flat, because the scripts resolve each other and their
outputs relative to a single directory. **The folder name is load-bearing** —
scripts compute their working directory as `$EBIRD_PROJ_ROOT/eBird`.

---

## 5. Instructions to replicators

### The short version

Download the three datasets (section 2), then:

```bash
export EBIRD_PROJ_ROOT=/path/to/project
Rscript main.R
```

`main.R` checks your packages and GDAL first, then runs the build, the analyses,
the figures and the regression checks in order, printing a timing for each. It
takes 12–15 hours. Individual stages:

```bash
Rscript main.R build       # data build only
Rscript main.R analysis    # analyses, assumes the build exists
Rscript main.R figures     # figures only
Rscript main.R verify      # regression checks only
```

Open the repository folder directly in your editor. **Do not call `setwd()`** —
`main.R` handles the one directory change the scripts need.

### The long version

**Step 1 — point the pipeline at your data.**

```bash
export EBIRD_PROJ_ROOT=/path/to/project    # scripts read $EBIRD_PROJ_ROOT/eBird
```

**Step 2 — build the tables, in this order.** Each depends on the one above it.

| # | Command | Time | Produces |
|---|---|---|---|
| 1 | `Rscript zerofill_from_registry.R` | hours | One zero-filled table per species |
| 2 | `Rscript join_species_all_years.R` | ~9 min | One row per checklist, joined to fire and vegetation |
| 3 | `Rscript regrid_cells.R` | ~6 min | The 5 km × 16-day cell frame the models read |

**Step 3 — run any analysis.** These are independent of each other and can be
run in any order, but **not at the same time** — two concurrent jobs will
exhaust memory.

| Command | Time | Answers |
|---|---|---|
| `Rscript dose_burn_effect.R` | ~6 min | Does the effect scale with how much of the square burned? |
| `Rscript event_study.R` | ~6 min | What happened each year after the fire, and were the pre-fire years flat? |
| `Rscript severity_effect.R` | ~6 min | Does severity matter on top of extent? |
| `Rscript effort_placebo.R` | ~1 min | Did the fire change how people birded? |
| `Rscript wild_cluster_boot.R` | ~35 min | The significance verdict, clustered by fire |
| `Rscript guild_contrast.R` | ~35 min | Do bird groups differ from each other? |
| `Rscript sign_test.R` | ~1 min | Did species move the way they were predicted to? |
| `./run_species.sh` | ~12 h | The whole pipeline, every species |

**Step 4 — check nothing moved.**

```bash
Rscript verify_published.R          # 41 checks; VERIFY_FULL=1 gives 51
```

Expect **40 of 41**. Check 6c is a known false positive: it looks for scripts
that call `cfg_*()` without sourcing the registry, and does not recognise
`sys.source()`.

**Step 5 — draw the figures.**

```bash
Rscript make_figures.R              # writes figures/fig[0-6]_*.png and .pdf
```

### Changing what is analysed

Species and treatment come from the registries, so nothing has a bird hardcoded
in it:

```bash
EBIRD_PREFIX=olsfly DOSE_MEASURE=two_dose_sev Rscript dose_burn_effect.R
```

Diagnostics write to their own files rather than overwriting a published result.
`guild_contrast.R` takes a list of species to exclude; **unset, it reproduces
the published run exactly**:

```bash
# refit the group comparisons without the three rarest snag species
GC_DROP=bkbwoo,moublu,lazbun GC_TAG=snagcommon Rscript guild_contrast.R
```

---

## 6. Description of programs

40 R scripts and 4 shell runners. The ones that matter:

**Building the data**

| Script | What it does |
|---|---|
| `parse_ebd_california.R` | Filters the raw 33 GB eBird file down to California |
| `zerofill_from_registry.R` | Turns sightings into presence and absence for every species in the registry |
| `join_species_all_years.R` | Joins checklists to fire perimeters and MODIS vegetation |
| `regrid_cells.R` | Aggregates checklists into 5 km × 16-day squares |

**Estimating effects**

| Script | What it does |
|---|---|
| `dose_burn_effect.R` | Treats fire as a dose — the share of the square that burned |
| `event_study.R` | Separates the effect by year, and checks that the pre-fire years are flat |
| `severity_effect.R` | Severity, holding extent fixed |
| `guild_contrast.R` | Compares whole groups of birds against each other |
| `sign_test.R` | Counts how many species moved the predicted way |

**Checking the design**

| Script | What it does |
|---|---|
| `effort_placebo.R` | Tests seven measures of birding effort for a fire response |
| `wild_cluster_boot.R` | The significance verdict, resampling whole fires |
| `fire_cluster_se.R` | Standard errors clustered by fire rather than by square |
| `cs_estimator.R` | Callaway–Sant'Anna estimator, as a check on the two-way fixed effects design |
| `verify_published.R` | Regression test pinning 41 published numbers |

**Shared code.** `infer.R` holds the inference core that every published number
flows through. `panel_utils.R` and `project_config.R` hold the shared helpers
and the registry loader.

---

## 7. Figures and the code that makes them

All figures come from `make_figures.R`, which **parses every number out of the
result files rather than restating it**. Each parser that reads a formatted text
report is followed by a check that re-derives a value the report prints for
itself, so a changed report layout stops the script instead of silently drawing
the wrong thing.

| Figure | Shows | Built from |
|---|---|---|
| 0 | How the two eBird files combine, and why zero-filling needs both | *Schematic — the only figure with invented data, and it says so on its face* |
| 1 | The Wrentit against the three water birds, year by year | `event_study_summary.txt` |
| 2 | The negative controls all register a fire effect | `point_controls_results.csv` |
| 3 | Fire changes how people bird | `effort_placebo_summary.txt` |
| 4 | The sign test and the group comparisons | `sign_test_results.csv`, `guild_contrast_results.csv` |
| 5 | Burn history predicts well for the wrong reason | `burn_effect_results_wrentit_covariate.csv` |
| 6 | Fire as a dose rather than a switch | `burn_dose_summary_wrentit.txt` |

---

## 8. Configuration

Four CSVs drive everything. Adding a species or a treatment is a row, not a code
change.

| File | Role |
|---|---|
| `species_registry.csv` | The 19 species: group, habitat, movement, and **the direction each was predicted to move, recorded before it was run** |
| `treatment_registry.csv` | 25 burn-impact measures across extent, severity, timing, reburn, and spatial attributes |
| `feature_registry.csv` | The 53 cell-frame columns by block and role. **Row order is load-bearing** — it is the design-matrix column order, and the tree models sample features off a seeded generator |
| `species_catalog.csv` | All 773 California eBird species with checklist counts, used to validate names |

`species_registry.csv` is the pre-registration. The predicted direction for each
species was committed before any model ran, which is what makes the sign test a
prediction rather than a pattern noticed afterwards.

---

## 9. Three mistakes that produced silent, wrong answers

Documented so nobody rediscovers them.

**Never join across species on `cell_id`.** It is assigned in row-appearance
order and is not stable between separately built species. Join on `(gx, gy)`.
This bug appeared three times.

**Never name a variable after a `data.table` column.** Four occurrences in one
day.

**Cluster by fire, not by square.** 3,125 burned squares come from 480 fires,
and the 2021 Dixie Fire alone accounts for 231 of them. Squares inside one
perimeter share one ignition, one weather system and one post-fire access
regime — that is one draw, not 231. Intervals clustered by square are too
narrow.

---

## 10. What the analysis found

**The design fails its own falsification test.** All three water-bird controls
register larger fire effects than the focal species. Migration timing fails the
same way: sedentary birds show larger apparent shifts than actual migrants.

**The mechanism is measurable.** After a square burns, birders visit 42% fewer
distinct locations inside it (p = 0.0001) and file 7% fewer checklists
(p = 0.013).

**Burn history predicts well, for the wrong reason.** It adds +0.046 AUC on
held-out years, and +0.002 once the model is told the square's coordinates. It
was largely a location detector, because fires burn chaparral and the focal
species lives in chaparral.

**Group comparisons survived longest, and are now qualified.** The
`control − snag` comparison is not silent, at −2.19, and every comparison
involving the snag group reduces to the same quantity — snag against the water
birds. Group membership was pre-registered; which comparisons to run was chosen
afterwards.

**What survives is the sign test.** 12 of 16 species moved in the direction
committed to `species_registry.csv` before the analysis ran — snag 6 of 6,
forest 4 of 4, shrub 2 of 6.

---

## 11. References

eBird Basic Dataset and Sampling Event Data. Cornell Lab of Ornithology, Ithaca,
New York. Cite the release you download.

Monitoring Trends in Burn Severity (MTBS), a joint program of the USGS EROS
Center and the USDA Forest Service.

MODIS/Terra Vegetation Indices 16-Day L3 Global 1 km SIN Grid, MOD13A2 version
061, distributed by the NASA LP DAAC.

---

## 12. Acknowledgements

Structure follows the replication-package conventions used by the
[Review of Financial Studies R template](https://github.com/review-of-financial-studies/template-r).
