# Data Dictionary — eBird × MTBS × MODIS (Wrentit)

Generated 2026-07-26. Covers every produced table, every column, and — in
detail — **the join criteria that connect the three sources**.

Species: **Wrentit** (`wrenti`, *Chamaea fasciata*). Region: California.
Study window: **2015-01-01 → 2026-06-30**. Fire history reaches back to **1984**.

---

## 0. Table index

| file | grain | rows | cols | size | built by |
|---|---|---|---|---|---|
| `wrentit_2015_2026_joined.csv` | one checklist × Wrentit | 4,445,934 | 25 | 611 MB | `join_species_all_years.R` |
| `wrentit_2015_2026_joined_clean.csv` | same, QA-masked | 4,012,521 | 25 | 559 MB | same |
| `wrentit_cells_5km_16d.csv` | one 5 km cell × 16-day step | 622,219 | 53 | 201 MB | `regrid_cells.R` |
| `burn_effect_results.csv` | one model fit | 20 | 8 | 1.7 KB | `model_burn_effect.R` |
| `wrentit_2015_2026_joined_summary.txt` | report | — | — | 16 KB | `join_species_all_years.R` |
| `wrentit_cells_5km_16d_summary.txt` | report | — | — | 5.7 KB | `regrid_cells.R` |
| `burn_effect_summary.txt` | report | — | — | — | `model_burn_effect.R` |
| `burn_did_summary_<prefix>.txt` | report | — | — | — | `did_burn_effect.R` |

**Lineage:** raw EBD → `zerofill_focal_species*.csv` → **checklist table** →
**cell × time frame** → models. The cell frame is derived *from* the checklist
table, not rebuilt from source.

---

## 1. JOIN CRITERIA (read this before using either table)

Six joins produce the checklist table. Each is stated as its exact rule.

### J1 — eBird observations → eBird sampling events (zero-fill)

**Performed upstream** by `auk_zerofill()` in `zerofill_focal_species.R` and
`filter_zerofill_2015_2023.R`.

- **Key:** `SAMPLING EVENT IDENTIFIER`
- **Criterion:** every **complete** checklist (`ALL SPECIES REPORTED = 1`) in
  California emits one row per focal species. If the species was not reported,
  the row is written with `presence = 0`.
- **Why it matters:** a 0 is a *verified absence* — an observer reported
  everything they detected and this species was not among it. It is **not** a
  missing value and must never be dropped or imputed.
- Filter applied at this stage: `auk_state("US-CA")`, `auk_complete()`,
  species ∈ {Wrentit, Brown Creeper, Olive-sided Flycatcher}.

### J2 — zero-fill → `protocol_type` (the S-id / G-id subtlety)

- **Left:** zero-fill table, key `checklist_id`
- **Right:** `ebd_focal_sampling_filtered*.txt`, keyed on **both**
  `SAMPLING EVENT IDENTIFIER` (S-ids) **and** `GROUP IDENTIFIER` (G-ids)
- **Criterion:** the lookup is the *union* of (S-id → protocol) and
  (G-id → protocol), de-duplicated by id.
- **Why both:** `auk_zerofill()` **re-keys shared checklists to the group id**.
  A lookup on S-ids alone silently loses every checklist birded by more than one
  person — a large, non-random slice. `checklist_id` values beginning `G` are
  group checklists; those beginning `S` are individual.
- **Protocol classification is EXACT string match**, not `grepl`:
  `"Traveling"` → Traveling, `"Stationary"` → Stationary, `"Historical"` →
  Historical, `"Incidental"` → Incidental, **everything else → Other**.
  This deliberately routes the banding-station variants
  (`"Traveling (2 band, 25m)"`, `"Stationary (3 band, 30m+100m)"`) and
  `"eBird Pelagic Protocol"`, `"Area"`, `"Banding"`,
  `"Nocturnal Flight Call Count"` into **Other**, which is then dropped. A
  substring test would have folded the banding variants into the kept classes.

### J3 — the two EBD downloads → one table

- 2015-01-01…2023-12-31 and 2024-01-01…2026-06-30 are concatenated.
- **De-duplicated on `checklist_id`**, keeping the first occurrence.
  The windows do not overlap, but **4 checklist_ids genuinely appear in both**
  downloads and are collapsed. The result is asserted unique.

### J4 — checklist point → MTBS fire perimeter *(the critical spatial-temporal join)*

- **Left:** checklist coordinates, WGS84 → **EPSG:5070** (NAD83 / Conus Albers,
  metres, equal-area)
- **Right:** `mtbs_perim/mtbs_perims_DD.shp`, filtered
  `WHERE Event_ID LIKE 'CA%'` → **1,956 CA perimeters, 1984-01-26 … 2025-09-03**,
  passed through `st_make_valid()`, transformed to EPSG:5070
- **Spatial predicate:** `st_within` — the point must be **strictly inside** the
  polygon. Not nearest-neighbour, not buffered.
- **Temporal predicate:** `Ig_Date < observation_date`, **STRICT**. A fire that
  ignited on or after the checklist date can never attach.
- **Tie-break:** among all perimeters satisfying both, take the one with the
  **maximum `Ig_Date`** — the most recent prior fire.
- **Result:** `fire_event_id` = that fire's `Event_ID`, or the literal string
  `"NONE"` if no perimeter qualifies.
- **Consequences to internalise:**
  - `days_since_fire ≥ 1` always; 0 and negative are impossible by construction.
  - A point inside a perimeter whose fire came *later* is a **control**
    (`"NONE"`), not a match. 29,169 such rows exist — genuine pre-fire
    observations at locations that burn afterwards.
  - **Only the single most recent fire is retained.** Reburn history is *not*
    represented at checklist grain. (The cell frame fixes this — see J8.)
  - Performance: candidate perimeters depend only on location, so `st_within`
    runs once per **unique coordinate** (498,566) and the date-ordered pick is
    then resolved per row.

### J5 — matched fire → burn severity class

- **Right:** `mtbs/mtbs_CA_<fire_year>/mtbs_CA_<fire_year>.tif` — the statewide
  30 m thematic mosaic **for the matched fire's own year**, already EPSG:5070.
- **Criterion:** point sample (`terra::extract`) at the checklist coordinate,
  from that specific year's raster. Sampled once per unique
  (location, fire_year) pair.
- **NA handling:** inside a perimeter but on an unmapped pixel → **0**.
  This is why `burn_severity_class == 0` does **not** mean "no fire" —
  8,034 matched rows carry class 0. **Always test `fire_event_id != "NONE"`.**
- **No `dnbr` column exists.** Verified from the raster itself: `INT1U`, carries
  a colour table, holds only values 1–6. It is the thematic *dnbr6* product, not
  continuous dNBR (which would be a scaled float ≈ −500…1300). Continuous dNBR
  requires per-fire `*_dnbr.tif` bundles that are not in this download.

### J6 — checklist → MODIS vegetation

- **Right:** `MOD13A2_061/*.hdf` — MOD13A2 v061, 1 km, 16-day, tiles
  h08v04 / h08v05 / h09v04, 2000-02-18 … 2026-06-26 (1,821 granules, no gaps).
- **Temporal criterion:** each checklist is assigned the **16-day composite
  window containing its date**:
  `comp_doy = 16 × ((day_of_year − 1) %/% 16) + 1`, capped at 353.
  Windows are DOY 1, 17, 33, …, 353 — 23 per year.
- **Spatial criterion:** checklist coordinates projected WGS84 → MODIS
  sinusoidal, then the **containing 1 km pixel** is sampled. Tiles are read in
  turn and coalesced (a point falls in exactly one tile; ocean fill and
  out-of-tile both return NA, so coalescing equals mosaicking).
- **Scaling:** raw integer × `1e-4`. GDAL reports MOD13A2's scale as `10000`
  because MODIS documents it as a **divisor**; `scoff()` is therefore stripped
  and scaling applied manually. Trusting terra's own scaling would inflate NDVI
  by 1e8.
- **Fill:** values < −2000 → NA (fill is −3000); `pixel_reliability` outside
  0–3 → NA; composite DOY outside 1–366 → NA.

### J7 — NDVI baseline (year-over-year reference)

- **Criterion:** the **same 16-day composite window, one calendar year earlier**,
  at the **same pixel**. `ndvi_delta = ndvi − ndvi_baseline`.
- This is a same-season comparison. The superseded 2024 table used a single
  2023 DOY-353 winter composite as everyone's baseline, which confounded season
  with fire effect.
- 2014 granules exist, so every 2015 row has a valid baseline window.

### J8 — checklist → 5 km cell × 16-day step *(cell frame only)*

- **Spatial:** coordinates → EPSG:5070, then
  `gx = floor(x / 5000)`, `gy = floor(y / 5000)`. Equal-area, so a cell is a
  true 5 km × 5 km square anywhere in the state.
- **Temporal:** the same `comp_doy` binning as J6, so the vegetation join is
  **exact rather than interpolated**.
- **Emission rule:** only cell-times **that received birding effort** are
  written. 622,219 of a possible 3,451,890 (13,026 cells × 265 steps) = **18.0%**.
  Unvisited cell-times carry no label and are absent by design.
- **Source table:** the **flagged** checklist table, not the clean one — the
  clean table drops rows on MODIS pixel quality, which is a vegetation-QA
  decision and must not delete birding effort from a cell's tally.

### J9 — cell → burn history *(cell frame only; supersedes J4/J5 at cell grain)*

- **Right:** all 42 `mtbs_CA_<yr>.tif` rasters, 1984–2025.
- **Criterion:** each raster is **streamed in row blocks**; every burned
  (non-NA, > 0) 30 m pixel is mapped to a cell by
  `gx = floor(x / 5000)`, `gy = floor(y / 5000)` — exact integer arithmetic,
  since raster and grid share EPSG:5070. Pixels are tallied by
  (cell, fire_year, severity class) and converted to area
  (`npx × 900 m² / 25,000,000 m²`).
- **Why not zonal statistics:** the 42 rasters have **42 different origins**
  (each is cropped to that year's fires), so no common lattice exists to
  aggregate onto; and materialising coordinates for a 789 M-cell raster exhausts
  memory.
- **Time gate:** a fire-year contributes to a step only when the **last ignition
  date among that cell's fires that year** is strictly earlier than the step
  start. Conservative — it can delay crediting an early fire in a cell that
  burned twice in one calendar year — but it can never credit a future fire.
  This preserves the J4 invariant at cell grain.
- **Reburns are retained.** Unlike J4, every fire year a cell has ever seen
  contributes. 99,559 units sit in cells with ≥ 2 fire years; the maximum is 11.
- **Discretisation:** 5000 / 30 = 166.67, so a cell holds 166 *or* 167 pixel
  centres per side depending on how the grid falls against each raster's origin.
  Counting each at a full 900 m² can overshoot nominal cell area by at most
  (167/166.67)² = **1.004**. `burn_frac_max` is clamped to 1; `burn_frac_cum` is
  left uncapped because it is cumulative across years anyway.

### J10 — cell → MODIS *(cell frame only)*

- **Criterion:** the 1 km pixel-centre → 5 km cell mapping is computed **once**;
  then each composite is read and averaged over the pixels whose centres fall in
  the cell (`mean` of valid NDVI; NA pixels excluded).
- A 5 km cell contains ~25 MODIS pixels, which is why cell-level NDVI
  missingness (0.4%) is far below point-level (5.5%).
- **Lags** are integer offsets on a sequential composite index
  `ci = (year − 2000) × 23 + ((comp_doy − 1) / 16) + 1`, so `l1` = the previous
  composite, `l3` = three back (~64 days ≈ "2 months up to this time"), and
  `ndvi_prevyr` = `ci − 23`.

### J11 — cell → neighbouring cells *(cell frame only)*

- **Criterion:** the 8 cells at `(gx ± 1, gy ± 1)` excluding self;
  `burn_frac_nbr` / `burn_decay_nbr` are the **mean over those neighbours at the
  same time step**.
- ⚠ **Known limitation:** only neighbours **that were themselves visited in that
  same 16-day window** contribute. `n_nbr` records how many did (0–8). The
  neighbourhood term is therefore an effort-biased sample of the surroundings,
  not a true spatial average. Treat `n_nbr` as a reliability weight.
- ⚠⚠ **`burn_frac_nbr = 0` is ambiguous.** 61,085 units (9.8%) have
  `n_nbr = 0` — no neighbour was visited that step — and those receive 0, which
  is indistinguishable from "neighbours were observed and none had burned".
  Only 32,398 units (5.2%) have the full 8. **Filter on `n_nbr > 0`, or
  interact the neighbour terms with `n_nbr`, before reading anything into
  them.** This is the weakest feature block in the frame.

---

## 2. QC filters (applied once, in `join_species_all_years.R`)

Rows must satisfy **all** of:

| filter | rule | dropped |
|---|---|---|
| protocol | `protocol_type ∈ {Traveling, Stationary}` | 60,654 |
| duration | `duration_minutes ≤ 300` or NA | 80,134 |
| distance | `effort_distance_km ≤ 5` or NA | 360,952 |
| observers | `number_observers ≤ 10` or NA | 80,045 |
| coordinates | non-missing | 0 |

Total dropped 521,755 of 4,967,689 (10.5%); criteria overlap, so the total is
their union. NA effort values are **kept** — NA means "not recorded", and for
Stationary checklists a missing distance is structurally correct.

**Deliberate spec deviation:** the written spec says keep `Historical`. The
later review instruction to drop it is followed instead — Historical rows are
retroactively entered records with 86.9% missing `effort_distance_km` and 60.9%
missing duration, and effort is a required covariate here.

### Flagged vs clean

`*_clean.csv` additionally requires: `ndvi` non-NA (this is also the **land
mask** — MODIS returns fill over ocean), `pixel_reliability ∉ {2,3}`, and
`burn_severity_class ≠ 6`. 433,413 rows removed.

---

## 3. `wrentit_2015_2026_joined.csv` — checklist grain (25 columns)

### Identity

| # | column | type | notes |
|---|---|---|---|
| 1 | `checklist_id` | string | **Primary key**, unique across all rows. `S…` individual, `G…` group (see J2). Never a feature. |
| 2 | `species_code` | string | Constant `wrenti`. Present so the schema survives changing `SPECIES`. Zero information — drop. |

### Location & time

| # | column | type | range | notes |
|---|---|---|---|---|
| 3 | `latitude` | float | 30.775 – 42.006 | WGS84 |
| 4 | `longitude` | float | −126.029 – −114.133 | WGS84. The western extreme is a pelagic checklist. |
| 5 | `observation_date` | date | 2015-01-01 – 2026-06-30 | |
| 6 | `year` | int | 2015–2026 | derived; 2026 is a half-year |
| 7 | `day_of_year` | int | 1–366 | derived; drives seasonality/detectability |

### Effort covariates

| # | column | type | range | notes |
|---|---|---|---|---|
| 8 | `protocol_type` | categorical | Traveling (3,090,519), Stationary (1,355,415) | only these two survive QC |
| 9 | `duration_minutes` | float | 0–300 | mean 57.5; 23 NA |
| 10 | `effort_distance_km` | float | 0–5 | mean 1.48; **30.48% NA = Stationary (structural)** |
| 11 | `number_observers` | int | 1–10 | mean 1.40 |

### Response

| # | column | type | notes |
|---|---|---|---|
| 12 | `presence` | binary | **The target.** 10.32% positive (10.63% clean). 0 = verified absence (see J1). |
| 13 | `observation_count` | int 0–82 | Individuals reported; eBird "X" → NA (5,073). **⚠ PERFECT LEAKAGE against `presence`** — exactly 0 on all 3,586,156 absences. Valid only as an abundance response, never as a predictor. |

### Fire (MTBS)

| # | column | type | notes |
|---|---|---|---|
| 14 | `fire_event_id` | string | MTBS `Event_ID` or `"NONE"`. 1,194 distinct fires. **The correct matched/control test.** High-cardinality identifier — not a feature. |
| 15 | `ignition_date` | date | Matched fire's `Ig_Date`. NA on controls (95.36%). Always < `observation_date`. |
| 16 | `fire_year` | int 1984–2025 | year of `ignition_date`; NA on controls |
| 17 | `burn_severity_class` | int 0–6 | 0 = no prior fire **or** unmapped pixel inside a perimeter (8,034 rows); 1 unburned-low, 2 low, 3 moderate, 4 high, **5 increased greenness**, 6 mask. **Not ordinal** — 5 is regrowth, outside the severity ladder. Treat as categorical. |
| 18 | `distance_to_perimeter_km` | float, **signed** | min −18.80, median 7.98, max 398.73. **Negative = INSIDE a perimeter**; magnitude is distance to the nearest perimeter *boundary* either way, so −5 and +5 are both "5 km from an edge". Boundaries simplified to 30 m ⇒ error ≤ 0.03 km. 12 rows are exactly 0 (on a boundary). |
| 19 | `days_since_fire` | int | 1 – 15,416 (42.2 yr), median 4,019. **NA = never burned, not unknown.** Never ≤ 0. |

### Vegetation (MODIS MOD13A2)

| # | column | type | notes |
|---|---|---|---|
| 20 | `ndvi` | float | −0.2 – 0.995, median 0.410. Raw × 1e-4. |
| 21 | `evi` | float | −0.153 – 0.871, mean 0.230. Less saturation-prone than NDVI in dense canopy. |
| 22 | `modis_composite_date` | date | The day the composite pixel actually represents; drifts from `observation_date` by up to ~16 days. Use to audit that drift. |
| 23 | `pixel_reliability` | int 0–3 | 0 good, 1 marginal, 2 snow/ice, 3 cloud. Clean keeps {0,1,NA}. |
| 24 | `ndvi_baseline` | float | −0.2 – 1.0, mean 0.420. Same window, previous year (J7). |
| 25 | `ndvi_delta` | float | `ndvi − ndvi_baseline`, median +0.004. 7.26% NA (union of inputs). |

**Missingness (flagged):** `ignition_date`/`fire_year`/`days_since_fire` 95.36%
(controls) · `effort_distance_km` 30.48% (Stationary) · MODIS block 5.52%
(ocean fill) · `ndvi_baseline` 5.50% · `ndvi_delta` 7.26% ·
`observation_count` 0.11% ("X") · everything else 0.00%.

---

## 4. `wrentit_cells_5km_16d.csv` — cell × time grain (53 columns)

### Keys (1–9)

| # | column | type | notes |
|---|---|---|---|
| 1 | `cell_id` | int | Surrogate id for a (gx, gy) pair. **Use as the CV grouping key.** Not a feature. |
| 2–3 | `gx`, `gy` | int | `floor(easting/5000)`, `floor(northing/5000)` in EPSG:5070. Cell SW corner = (gx×5000, gy×5000). |
| 4 | `year` | int | 2015–2026 |
| 5 | `comp_doy` | int | 1, 17, …, 353 — MODIS composite window start |
| 6 | `ci` | int | Sequential composite index (J10). Lags are offsets on this. |
| 7 | `step_start` | date | `Jan 1 of year + (comp_doy − 1)` — first day of the step |
| 8–9 | `date_first`, `date_last` | date | Actual first/last checklist date within the unit |

**Grain is (`cell_id`, `year`, `comp_doy`)** — asserted unique. `ci` is a
one-to-one recoding of (`year`, `comp_doy`).

### Target & effort (10–19)

| # | column | type | notes |
|---|---|---|---|
| 10 | `any_detection` | binary | 1 if Wrentit on **any** checklist in the unit. 22.1% positive. **⚠ Rises mechanically with effort** — 0.096 at 1 checklist → 0.489 at 20+. |
| 11 | `n_detections` | int | Checklists in the unit reporting Wrentit. **Leaks the target.** |
| 12 | `n_checklists` | int | Checklists in the unit. Feature *and* the effort denominator. |
| 13 | `det_rate` | float | `n_detections / n_checklists`, mean 0.107. **Effort-normalised target.** Also leaks if used as a feature. |
| 14 | `n_locations` | int | Distinct lat/lon within the unit |
| 15 | `effort_hours` | float | `Σ duration_minutes / 60` |
| 16 | `dur_median` | float | Median checklist duration |
| 17 | `dist_median` | float | Median distance; **16.51% NA = every checklist in the unit was Stationary** |
| 18 | `obs_median` | float | Median party size |
| 19 | `frac_traveling` | float 0–1 | Share of checklists on the Traveling protocol |

> **Target choice.** `any_detection` is the binary task as specified, but the
> table above shows it is not effort-neutral. `det_rate` already divides effort
> out. Use `det_rate` (or a fixed-effort subsample) when the question is
> ecological; use `any_detection` with effort covariates when the question is
> operational forecasting.

### Burn history (20–42) — all strictly time-gated by J9

| # | column | type | notes |
|---|---|---|---|
| 20 | `ever_burned` | binary | 1 if any prior fire in-cell. 33.7% of units. |
| 21 | `burn_frac_cum` | float ≥ 0 | **Cumulative** burned fraction, summed over all prior fire years. **> 1 means reburn, not an error** — max 4.47; 27,729 units exceed 1. Median 0.287 where > 0. |
| 22 | `burn_frac_max` | float 0–1 | Largest **single fire-year** burned fraction. Bounded (clamped, see J9). |
| 23 | `burn_sev_wt` | float | `Σ frac × sev` over prior fires, where `sev` counts only classes 1–4 (5 = greenness and 6 = mask contribute 0 to severity but **do** contribute to `burn_frac_*`). |
| 24 | `burn_decay` | float | `Σ frac × sev × exp(−yrs_since / 10)` — recency-weighted severity, τ = 10 yr. |
| 25 | `yrs_since_burn` | float | Years since the most recent prior in-cell fire. **66.30% NA = never burned.** Median 8.84, max 42.07. |
| 26 | `max_sev_prior` | int 0–4 | Highest severity class ever recorded in-cell. Strongest single burn feature in permutation importance. |
| 27 | `n_fire_epochs` | int 0–11 | Distinct prior fire **years** — the reburn counter. |
| 28–33 | `burn_frac_y0_1` … `y20p` | float | Burned fraction by epoch: 0–1, 1–3, 3–5, 5–10, 10–20, 20+ years before the step. Cumulative within epoch. |
| 34–39 | `burn_sev_y0_1` … `y20p` | float | Same epochs, severity-weighted. |
| 40 | `burn_frac_nbr` | float | Mean `burn_frac_cum` of the 8 neighbours **present at this step** (J11). |
| 41 | `burn_decay_nbr` | float | Same for `burn_decay`. |
| 42 | `n_nbr` | int 0–8 | How many neighbours contributed. **Reliability weight for cols 40–41.** |

All burn columns are **0** (not NA) where no prior fire applies; only
`yrs_since_burn` is NA, because 0 would be a meaningful value there.

### Vegetation sequence (43–53)

| # | column | type | notes |
|---|---|---|---|
| 43–46 | `ndvi_l0` … `ndvi_l3` | float | Cell-mean NDVI for this step and the 3 preceding composites (~64 days). Medians 0.433, 0.432, 0.429, 0.425. **0.42% NA.** |
| 47–50 | `evi_l0` … `evi_l3` | float | Same lags, EVI |
| 51 | `ndvi_prevyr` | float | Same composite window one year earlier (`ci − 23`) |
| 52 | `ndvi_delta_yr` | float | `ndvi_l0 − ndvi_prevyr` — year-over-year change |
| 53 | `ndvi_trend` | float | `ndvi_l0 − ndvi_l3` — change across the 2-month run-up |

---

## 5. `burn_effect_results.csv` — model results (8 columns)

| column | notes |
|---|---|
| `split` | `temporal` (train ≤ 2023, test 2024–2026) or `spatial` (cells 75/25, no cell in both) |
| `algo` | `lightgbm`, `xgboost`, `gbm`, `glm`, or `lightgbm+geo` (the geography probe, which adds `gx`/`gy`) |
| `variant` | `full` (43 features) or `noburn` (21 — burn block removed) |
| `n_feat` | Feature count |
| `auc` | Test ROC-AUC |
| `logloss` | Test log loss |
| `brier` | Test Brier score |
| `secs` | Fit seconds (NA for the geo probe) |

**Δ-AUC = `full` − `noburn` within a (split, algo)** is the headline: what the
burn block buys once effort, season and vegetation are present.

### Reading the headline results

| | Δ-AUC (temporal) |
|---|---|
| LightGBM | +0.0460 |
| XGBoost | +0.0434 |
| gbm | +0.0340 |
| logistic | +0.0338 |
| **LightGBM + geography** | **+0.0021** |

Burn history is worth ~+0.04 AUC **until the model is told where the cell is**,
at which point it adds ~0. The burn block was largely proxying for location,
i.e. for chaparral.

`gbm` trains on a 150,000-row subsample (single-threaded, far slower) — it is a
reference implementation, not a tuned competitor.

---

## 6. Cross-cutting cautions

1. **Never split at random.** Cells and locations repeat heavily — 50% of
   checklist rows sit on 0.48% of locations. Use `cell_id` grouping or a
   temporal split.
2. **Exclude from any feature matrix:** `checklist_id`, `species_code`,
   `observation_count`, `fire_event_id` (checklist table); `cell_id`, `ci`,
   `step_start`, `date_first`, `date_last`, `n_detections`, `det_rate`
   (cell table, when `any_detection` is the target).
3. **`burn_severity_class == 0` ≠ unburned.** Test `fire_event_id != "NONE"`.
4. **Effective sample size is cells/locations, not rows** — 36,519 distinct
   locations carry the 206,379 fire-matched checklist rows; 1,009 cells carry
   the within-cell fire estimate. Uncertainty should be computed at that level.
5. **Treatment is diluted at 5 km.** Median burned fraction among burned cells
   is 0.29 — a typical "treated" cell is ~70% unburned. Effects are attenuated
   roughly 3× relative to point grain.
6. **MTBS under-maps 2023–2025** (5 CA fires mapped for 2025 vs a 2015–2021 mean
   of 53/yr). Recent fires are missing, so some "controls" burned unrecorded.
7. **Detection ≠ occupancy.** Every target here is *observation*, jointly
   determined by presence, availability and detectability.

---

## 7. Reproduction

```bash
Rscript join_species_all_years.R    # ~9 min  -> checklist table
Rscript regrid_cells.R      # ~6 min  -> cell x time frame
Rscript model_burn_effect.R         # ~9 min  -> predictive results
Rscript did_burn_effect.R           # ~2 min  -> within-cell estimates
```

Paths resolve from `EBIRD_PROJ_ROOT` (default
`/Users/ethanle/Downloads/eBird_project`). Checkpoints live in
`.join_cache_all_years/` and `.regrid_cache/`; delete to force a full recompute.
`EBIRD_SAMPLE_N` and `REGRID_SAMPLE_YEARS` run fast subset smoke tests.

Requires R ≥ 4.5 with `terra`, `sf`, `data.table`, `lubridate` (and
`lightgbm`, `xgboost`, `gbm` for the models). **GDAL must have the HDF4 driver**
or the MODIS granules cannot be opened.
