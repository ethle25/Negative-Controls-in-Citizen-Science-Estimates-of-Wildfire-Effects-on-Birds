# Dependencies

Versions recorded on the machine the analysis was last run on.

**R 4.5.2 (2025-10-31 ucrt)**, platform `x86_64-w64-mingw32`, Windows 11.

## Packages

| Package | Version | Needed for |
|---|---|---|
| `data.table` | 1.17.8 | Every script — all tables are data.tables |
| `terra` | 1.9.34 | MTBS rasters and MODIS granules |
| `sf` | 1.1.1 | Fire perimeters, spatial joins |
| `lubridate` | 1.9.4 | Dates and the 16-day composite calendar |
| `auk` | 0.9.1 | eBird filtering and zero-filling |
| `ggplot2` | 4.0.1 | Figures |
| `patchwork` | 1.3.2 | Multi-panel figure layout |
| `ragg` | 1.5.0 | PNG device — the default device drops the typography in the captions |
| `lightgbm` | **not currently installed** | `model_burn_effect.R`, prediction |
| `xgboost` | **not currently installed** | Prediction comparison |
| `gbm` | **not currently installed** | Prediction comparison |
| `torch` | **not currently installed** | `deep_models.R` only |

> **The four model packages are not installed on this machine.** The prediction
> and deep-learning scripts therefore cannot run here as it stands, and their
> versions are not recorded. The published prediction numbers were produced when
> those packages were present. **Reinstall them and record the versions before
> publishing this repository** — a replicator who follows the README will hit
> this immediately.

```r
install.packages(c("lightgbm", "xgboost", "gbm"))
# torch additionally downloads libtorch on first use:
install.packages("torch"); torch::install_torch()
```

## System requirements

**Memory.** 24 GB is the practical floor. The eBird filtering step streams a
33 GB text file, and materialising the MODIS raster coordinates has exhausted
24 GB.

**GDAL with the HDF4 driver.** Without it the MODIS granules will not open:

```r
any(grepl("HDF4", terra::gdal(drivers = TRUE)$name))   # must be TRUE
```

On Windows, the GDAL that ships with `terra` normally includes HDF4. On Linux,
`libgdal-dev` often does not — check before starting a multi-hour build.

## Recording your own environment

`main.R` writes `results/sessionInfo.txt` at the end of every run. Include that
file when reporting a result that does not reproduce.

## A note on pinning

This repository does not use `renv`. The analysis was developed against the
versions above over roughly two months, and adding a lockfile after the fact
would pin whatever happens to be installed now rather than what produced the
published numbers — which, given the four missing packages, would be misleading.

The honest safeguard here is `verify_published.R`, which pins 41 published
numbers directly. If a package version change moves a result, that script fails,
which is the thing a lockfile is meant to prevent. Run it before and after any
change.
