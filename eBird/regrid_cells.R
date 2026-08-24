# =============================================================================
# Wrentit cell x time-step prediction frame
#
# Converts the checklist-grain joined table into the N x T prediction setting:
# California on an equal-area GRID, one row per (grid cell, 16-day time step)
# THAT ACTUALLY RECEIVED BIRDING EFFORT. Unvisited cell-times are not emitted --
# with no observer there is no observation, so they carry no label.
#
#   grid  : 5 km equal-area cells, EPSG:5070 (NAD83 / Conus Albers)
#   time  : 16-day steps aligned to the MODIS composite calendar
#           (DOY 1, 17, 33, ... 353; 23 per year), so the vegetation join is
#           exact rather than interpolated
#   span  : 2015-01-01 .. 2026-06-30  (~598k cell-time units, 12.5k cells)
#
# TARGET
#   any_detection    1 if Wrentit was reported on ANY checklist in the unit
#   n_detections / n_checklists  -> detection rate, if a rate target is wanted
#   CAUTION: any_detection rises mechanically with effort (more checklists ->
#   more chances). Effort aggregates are emitted alongside so this is modelled,
#   not ignored; det_rate and the fixed-effort flag are the alternatives.
#
# FEATURE FAMILIES
#   effort_*   n_checklists, observer-hours, distance, party size, protocol mix
#   ndvi_l0..l3 / evi_l0  vegetation SEQUENCE: the cell-mean composite for this
#              step and the 3 preceding ones (~64 days = the "2 months up to
#              this time"), plus the same window one year earlier
#   burn_*     cell burn HISTORY from the 30 m severity rasters: what fraction
#              of the cell burned, at what severity, how long ago -- bucketed
#              into epochs, plus a distance/time-decayed weight and the same
#              quantities for the 8 neighbouring cells
#
# WHY CELL-GRAIN BURN FEATURES ARE BETTER THAN THE CHECKLIST ONES
#   The checklist table asks "is this POINT inside a perimeter", which (a) sees
#   only the single most recent containing fire, so REBURNS are invisible, and
#   (b) reads 0 for a point 200 m outside a huge high-severity fire. Aggregating
#   the 30 m rasters over a cell instead gives burned FRACTION per severity
#   class per fire-year -- every fire the cell has ever seen, weighted by how
#   much of the cell each one took. The neighbour terms then supply the
#   "weighted distance to burn severity" the prediction setting asks for.
#
# TIME-ORDER INVARIANT (carried over from the checklist join)
#   A cell-time step may only see fires that finished burning BEFORE the step
#   starts. Enforced with the per-cell, per-year LAST ignition date: a fire-year
#   contributes to a step only once max(Ig_Date) for that cell-year is strictly
#   earlier than the step start. Conservative -- it can delay crediting an early
#   fire in a cell that burned twice in one calendar year, which is rare -- but
#   it can never credit a fire that has not happened yet.
#
# OUTPUT
#   wrentit_cells_5km_16d.csv         the prediction frame
#   wrentit_cells_5km_16d_summary.txt row counts, coverage, missingness
#
# Checkpoints in .regrid_cache/ so an interrupted run resumes.
# Set REGRID_SAMPLE_YEARS to run on a subset of years for a smoke test.
# =============================================================================
suppressMessages({
  library(terra); library(sf); library(data.table); library(lubridate)
})
sf_use_s2(FALSE)
terraOptions(progress = 0)
setDTthreads(0)

t_start <- Sys.time()
log <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

# --- config -------------------------------------------------------------------
GRID_M    <- 5000L        # cell edge, metres (equal-area)
N_LAG     <- 3L           # extra composites back (l0..l3 = 4 steps ~ 64 days)
COMP_PER_YR <- 23L        # MODIS composites per year
DECAY_YRS <- 10           # tau for the exponential burn-recency weight
EPOCHS    <- list(y0_1 = c(0, 1), y1_3 = c(1, 3), y3_5 = c(3, 5),
                  y5_10 = c(5, 10), y10_20 = c(10, 20), y20p = c(20, 999))

proj_root <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
ebird_dir <- file.path(proj_root, "eBird")
mtbs_dir  <- file.path(proj_root, "mtbs")
modis_dir <- file.path(proj_root, "MOD13A2_061")
perim_shp <- file.path(ebird_dir, "mtbs_perim", "mtbs_perims_DD.shp")
# No species is hardcoded: identity comes from species_registry.csv via
# cfg_current_species(), which honours EBIRD_PREFIX when set and otherwise takes
# the first ACTIVE registry row. The old default was the literal "wrentit".
source(file.path(ebird_dir, "project_config.R"))
PREFIX    <- Sys.getenv("EBIRD_PREFIX", cfg_current_species()$prefix[1])
in_csv    <- file.path(ebird_dir, sprintf("%s_2015_2026_joined.csv", PREFIX))

out_csv <- file.path(ebird_dir, sprintf("%s_cells_%dkm_16d.csv", PREFIX, GRID_M / 1000))
out_sum <- file.path(ebird_dir, sprintf("%s_cells_%dkm_16d_summary.txt", PREFIX, GRID_M / 1000))
cache   <- file.path(ebird_dir, ".regrid_cache")
dir.create(cache, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(in_csv), file.exists(perim_shp),
          dir.exists(mtbs_dir), dir.exists(modis_dir))
log("grid ", GRID_M / 1000, " km | 16-day steps | lags l0..l", N_LAG)

# =============================================================================
# 1. cell x time frame from the checklist table: target + effort
#    The FLAGGED table is used, not the clean one: the clean table drops rows
#    for MODIS pixel quality, which is a vegetation-QA decision and must not
#    silently delete birding effort from a cell's tally.
# =============================================================================
# The three focal species share the SAME checklists (the zero-fill emits one row
# per checklist per species) and QC filters on checklist-level fields only, so
# the surviving cell set, effort and geometry are IDENTICAL across species --
# only the detection columns differ. veg.rds and burn.rds therefore stay shared
# (they are the expensive stages); only the frame is per-species.
frame_rds <- file.path(cache, sprintf("frame_%s.rds", PREFIX))
if (file.exists(frame_rds)) {
  log("Loading cached cell-time frame"); fr <- readRDS(frame_rds)
  cells <- fr$cells; ct <- fr$ct
} else {
  log("Reading checklist table")
  d <- fread(in_csv, showProgress = FALSE,
             select = c("checklist_id","latitude","longitude","observation_date",
                        "year","day_of_year","protocol_type","duration_minutes",
                        "effort_distance_km","number_observers","presence"))
  yrs_keep <- Sys.getenv("REGRID_SAMPLE_YEARS", "")
  if (nzchar(yrs_keep)) {
    keep <- as.integer(strsplit(yrs_keep, ",")[[1]])
    d <- d[year %in% keep]; log("SMOKE TEST years: ", yrs_keep, " -> ", fmt(nrow(d)), " rows")
  }

  # project unique coordinates once, then assign the equal-area cell
  log("Assigning ", GRID_M / 1000, " km cells")
  uc <- unique(d[, .(longitude, latitude)])
  xy <- st_coordinates(st_transform(
          st_as_sf(uc, coords = c("longitude","latitude"), crs = 4326), 5070))
  uc[, `:=`(gx = as.integer(floor(xy[,1] / GRID_M)),
            gy = as.integer(floor(xy[,2] / GRID_M)))]
  d[uc, on = .(longitude, latitude), `:=`(gx = i.gx, gy = i.gy)]

  # 16-day step aligned to the MODIS composite calendar
  d[, comp_doy := 16L * ((day_of_year - 1L) %/% 16L) + 1L]
  d[comp_doy > 353L, comp_doy := 353L]

  log("Aggregating effort and target per cell-time")
  ct <- d[, .(
    n_checklists   = .N,
    n_detections   = sum(presence),
    any_detection  = as.integer(any(presence == 1L)),
    det_rate       = round(mean(presence), 6),
    n_locations    = uniqueN(paste(latitude, longitude)),
    # effort: total observer-hours is the natural cell-level effort unit
    effort_hours   = round(sum(duration_minutes, na.rm = TRUE) / 60, 4),
    dur_median     = as.numeric(median(duration_minutes, na.rm = TRUE)),
    dist_median    = as.numeric(median(effort_distance_km, na.rm = TRUE)),
    obs_median     = as.numeric(median(number_observers, na.rm = TRUE)),
    frac_traveling = round(mean(protocol_type == "Traveling"), 4),
    date_first     = min(observation_date), date_last = max(observation_date)
  ), by = .(gx, gy, year, comp_doy)]
  ct[is.nan(dist_median), dist_median := NA_real_]

  cells <- unique(ct[, .(gx, gy)])
  # cell_id is assigned by ROW-APPEARANCE order, so it depends on the order of
  # the source joined table. Species built from the same combined zero-fill share
  # a numbering; one built from a per-species zero-fill does NOT (calthr matches
  # wrentit on 0 of 13,026 cells). Never join on cell_id across species -- use
  # (gx, gy), as attach_fire() does.
  #
  # CELL_ID_SORTED=1 numbers by sorted (gx, gy) instead, making it identical for
  # every species. NOT the default, and the reason matters: control cells receive
  # a RANDOM pseudo event year drawn in cell order, so renumbering reshuffles that
  # draw and shifts every published estimate -- without making any of them more
  # correct. Flip it only as part of a deliberate re-baseline, then re-pin every
  # expected value in verify_published.R.
  if (Sys.getenv("CELL_ID_SORTED", "0") != "0") {
    setorder(cells, gx, gy)
    log("  CELL_ID_SORTED=1: numbering cells by sorted (gx,gy) -- ",
        "published numbers WILL move; re-pin verify_published.R")
  }
  cells[, cell_id := .I]
  ct[cells, on = .(gx, gy), cell_id := i.cell_id]
  setorder(ct, cell_id, year, comp_doy)
  saveRDS(list(cells = cells, ct = ct), frame_rds)
  rm(d, uc); invisible(gc())
}
n_cells <- nrow(cells); n_ct <- nrow(ct)
# order-sensitive fingerprint of the cell_id -> (gx,gy) mapping
CELLFP <- c(GRID_M, n_cells,
            sum(as.numeric(cells$gx) * cells$cell_id),
            sum(as.numeric(cells$gy) * cells$cell_id))

# Short tag identifying THIS numbering. The caches are keyed by cell_id, so a
# species with a different numbering must not overwrite another's cache -- doing
# so silently re-keys it and every cache-sourced burn measure then mis-joins for
# every other species, with no error. (Happened 2026-07-31: the calthr regrid
# re-keyed burn.rds while wrentit/brncre/olsfly still expected their own.)
# One cache file per numbering; species sharing a numbering share the cache, which
# is what makes the wrentit/brncre/olsfly builds cheap.
CELLTAG <- substr(sprintf("%08x", abs(sum(CELLFP) %% .Machine$integer.max)), 1, 8)
log("cell-time units: ", fmt(n_ct), " on ", fmt(n_cells), " cells")

# sequential composite index, so a lag is just an integer offset
ct[, ci := (year - 2000L) * COMP_PER_YR + ((comp_doy - 1L) %/% 16L) + 1L]
ci_need <- sort(unique(c(outer(ct$ci, 0:N_LAG, "-"), ct$ci - COMP_PER_YR)))
ci_need <- ci_need[ci_need > 0]
ci2yd <- function(ci) list(year = 2000L + (ci - 1L) %/% COMP_PER_YR,
                           doy  = 16L * ((ci - 1L) %% COMP_PER_YR) + 1L)

# =============================================================================
# 2. cell-level vegetation: aggregate 1 km MOD13A2 to the grid, per composite
#    The 1km-pixel -> grid-cell mapping is fixed, so it is built ONCE and then
#    every composite is a read + a group-mean over that mapping.
# =============================================================================
veg_rds <- file.path(cache, sprintf("veg_%s.rds", CELLTAG))
if (!file.exists(veg_rds) && file.exists(file.path(cache, "veg.rds")) &&
    identical(readRDS(file.path(cache, "veg.rds"))$key, c(CELLFP, length(ci_need))))
  file.rename(file.path(cache, "veg.rds"), veg_rds)   # migrate the legacy unnamed cache
if (file.exists(veg_rds) && identical(readRDS(veg_rds)$key, c(CELLFP, length(ci_need)))) {
  log("Loading cached cell vegetation"); veg <- readRDS(veg_rds)$veg
} else {
  MG <- "MODIS_Grid_16DAY_1km_VI"
  sds <- function(h, l) sprintf('HDF4_EOS:EOS_GRID:"%s":%s:"%s"', h, MG, l)
  L_NDVI <- "1 km 16 days NDVI"; L_EVI <- "1 km 16 days EVI"

  hdfs <- list.files(modis_dir, pattern = "^MOD13A2\\.A\\d{7}\\.h\\d{2}v\\d{2}\\..*\\.hdf$",
                     full.names = TRUE)
  bn <- basename(hdfs)
  meta <- data.table(path = hdfs, ayear = as.integer(substr(bn, 10, 13)),
                     adoy = as.integer(substr(bn, 14, 16)), tile = substr(bn, 18, 23))
  setkey(meta, ayear, adoy, tile)
  tiles <- sort(unique(meta$tile))

  log("Building the 1 km pixel -> ", GRID_M / 1000, " km cell map")
  pxmap <- rbindlist(lapply(seq_along(tiles), function(i) {
    r  <- rast(sds(meta[tile == tiles[i]]$path[1], L_NDVI))
    xy <- xyFromCell(r, seq_len(ncell(r)))
    p  <- project(vect(xy, type = "points", crs = crs(r)), "EPSG:5070")
    g  <- crds(p)
    data.table(tile_i = i, cell = seq_len(ncell(r)),
               gx = as.integer(floor(g[,1] / GRID_M)),
               gy = as.integer(floor(g[,2] / GRID_M)))
  }))
  pxmap <- pxmap[cells, on = .(gx, gy), nomatch = NULL]   # occupied cells only
  log("  1 km pixels landing in occupied cells: ", fmt(nrow(pxmap)))

  read_layer <- function(f, l) { r <- rast(sds(f, l)); scoff(r) <- NULL
                                 as.vector(values(r)) }
  veg <- vector("list", length(ci_need))
  for (j in seq_along(ci_need)) {
    yd <- ci2yd(ci_need[j]); yy <- yd$year; dd <- yd$doy
    parts <- list()
    for (i in seq_along(tiles)) {
      f <- meta[.(yy, dd, tiles[i]), path, nomatch = NULL]
      if (!length(f)) next
      pm <- pxmap[tile_i == i]
      if (!nrow(pm)) next
      nd <- read_layer(f[1], L_NDVI)[pm$cell]; nd[nd < -2000] <- NA
      ev <- read_layer(f[1], L_EVI )[pm$cell]; ev[ev < -2000] <- NA
      parts[[i]] <- data.table(cell_id = pm$cell_id, nd = nd * 1e-4, ev = ev * 1e-4)
    }
    if (!length(parts)) next
    v <- rbindlist(parts)[!is.na(nd), .(ndvi = mean(nd), evi = mean(ev, na.rm = TRUE),
                                        npx = .N), by = cell_id]
    v[, ci := ci_need[j]]
    veg[[j]] <- v
    if (j %% 25 == 0 || j == length(ci_need))
      log("  composite ", j, "/", length(ci_need), " (", yy, " DOY ", dd, ")")
  }
  veg <- rbindlist(veg)
  setkey(veg, cell_id, ci)
  saveRDS(list(key = c(CELLFP, length(ci_need)), veg = veg), veg_rds)
}
log("cell-composite vegetation rows: ", fmt(nrow(veg)))

# =============================================================================
# 3. cell burn history from the 30 m severity rasters
#    The rasters are already EPSG:5070, so a burned pixel's grid cell is exact
#    integer division -- no resampling, no polygon zonal stats.
#    Per-cell-year LAST ignition date comes from the perimeters and enforces the
#    time-order invariant in step 4.
# =============================================================================
burn_rds <- file.path(cache, sprintf("burn_%s.rds", CELLTAG))
if (!file.exists(burn_rds) && file.exists(file.path(cache, "burn.rds")) &&
    identical(readRDS(file.path(cache, "burn.rds"))$grid, CELLFP))
  file.rename(file.path(cache, "burn.rds"), burn_rds)  # migrate the legacy unnamed cache
if (file.exists(burn_rds) && identical(readRDS(burn_rds)$grid, CELLFP)) {
  log("Loading cached cell burn history"); bh <- readRDS(burn_rds)
  cellburn <- bh$cellburn; cellig <- bh$cellig
} else {
  PX_AREA <- 900                       # 30 m x 30 m
  CELL_AREA <- as.numeric(GRID_M)^2
  yrs <- sort(as.integer(sub(".*mtbs_CA_(\\d{4})\\.tif$", "\\1",
        list.files(mtbs_dir, pattern = "mtbs_CA_\\d{4}\\.tif$", recursive = TRUE))))

  # Each year's mosaic is cropped to that year's fires, so the 42 rasters have
  # 42 different origins -- there is no common lattice to aggregate onto. Stream
  # each raster in row blocks instead and map burned pixels to grid cells by
  # arithmetic (both are EPSG:5070, so this is exact). Materialising the whole
  # raster's coordinates would need far more memory than the machine has:
  # the 2020 mosaic alone is 789M cells.
  tally_burn <- function(tif) {
    r  <- rast(tif)
    nc <- ncol(r); nr <- nrow(r)
    x0 <- xmin(r); y1 <- ymax(r); rx <- xres(r); ry <- yres(r)
    rows_per <- max(1L, as.integer(1e7 %/% nc))
    acc <- list()
    readStart(r); on.exit(readStop(r), add = TRUE)
    for (s in seq(1L, nr, by = rows_per)) {
      nrows <- min(rows_per, nr - s + 1L)
      v <- readValues(r, s, nrows, 1, nc)
      idx <- which(!is.na(v) & v > 0)
      if (!length(idx)) { rm(v, idx); next }
      rr <- s + (idx - 1L) %/% nc                 # global row (1-based)
      cc <- ((idx - 1L) %% nc) + 1L               # global col
      acc[[length(acc) + 1L]] <- data.table(
        gx  = as.integer(floor((x0 + (cc - 0.5) * rx) / GRID_M)),
        gy  = as.integer(floor((y1 - (rr - 0.5) * ry) / GRID_M)),
        sev = as.integer(v[idx]))[, .N, by = .(gx, gy, sev)]
      rm(v, idx, rr, cc)
    }
    if (!length(acc)) return(NULL)
    rbindlist(acc)[, .(npx = sum(N)), by = .(gx, gy, sev)]
  }

  log("Tallying burned pixels per cell for ", length(yrs), " fire years")
  cb <- vector("list", length(yrs))
  for (k in seq_along(yrs)) {
    yr  <- yrs[k]
    tif <- file.path(mtbs_dir, sprintf("mtbs_CA_%d", yr), sprintf("mtbs_CA_%d.tif", yr))
    a <- tally_burn(tif)
    if (is.null(a)) { log("  ", yr, ": no burned pixels"); next }
    a <- a[cells, on = .(gx, gy), nomatch = NULL]
    if (nrow(a)) {
      a[, `:=`(fire_year = yr, frac = (npx * PX_AREA) / CELL_AREA)]
      cb[[k]] <- a[, .(cell_id, fire_year, sev, frac)]
    }
    log("  ", yr, ": ", fmt(if (is.null(cb[[k]])) 0L else nrow(cb[[k]])),
        " (cell,severity) tallies")
    invisible(gc())
  }
  cellburn <- rbindlist(cb); rm(cb); invisible(gc())

  # per-cell-year last ignition date, from the perimeters
  log("Per-cell-year last ignition date")
  perim <- st_read(perim_shp, quiet = TRUE, query =
    "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
  perim <- st_make_valid(perim); perim$Ig_Date <- as.Date(perim$Ig_Date)
  p5070 <- st_transform(perim, 5070)
  # rasterise each perimeter onto the grid by intersecting with cell polygons
  cg <- copy(cells)
  geom <- st_sfc(lapply(seq_len(nrow(cg)), function(i) {
    x0 <- cg$gx[i] * GRID_M; y0 <- cg$gy[i] * GRID_M
    st_polygon(list(cbind(c(x0, x0 + GRID_M, x0 + GRID_M, x0, x0),
                          c(y0, y0, y0 + GRID_M, y0 + GRID_M, y0))))
  }), crs = 5070)
  cellpoly <- st_sf(cell_id = cg$cell_id, geometry = geom)
  hit <- st_intersects(cellpoly, p5070)
  cellig <- rbindlist(lapply(which(lengths(hit) > 0), function(i) {
    p <- hit[[i]]
    data.table(cell_id = cellpoly$cell_id[i], fire_year = year(p5070$Ig_Date[p]),
               ig = p5070$Ig_Date[p])
  }))[, .(last_ig = max(ig)), by = .(cell_id, fire_year)]
  saveRDS(list(grid = CELLFP, cellburn = cellburn, cellig = cellig), burn_rds)
}
log("cell-burn tallies: ", fmt(nrow(cellburn)), " | cell-years with ignition: ",
    fmt(nrow(cellig)))

# =============================================================================
# 4. roll the burn history up to each cell-time, respecting time order
# =============================================================================
log("Building time-relative burn features")
cb <- merge(cellburn, cellig, by = c("cell_id", "fire_year"), all.x = TRUE)
cb <- cb[!is.na(last_ig)]
# severity 1-4 carry the burn ladder; 5 is regrowth and 6 is a mask
cb[, sev_w := fifelse(sev >= 1L & sev <= 4L, as.numeric(sev), 0)]

ct[, step_start := as.Date(sprintf("%d-01-01", year)) + (comp_doy - 1L)]
setkey(cb, cell_id)

steps <- unique(ct[, .(cell_id, year, comp_doy, step_start)])
setkey(steps, cell_id)
j <- cb[steps, on = "cell_id", allow.cartesian = TRUE][last_ig < step_start]
j[, yrs_since := as.numeric(step_start - last_ig) / 365.25]

agg <- j[, .(
  burn_frac_cum   = sum(frac),          # cumulative; > 1 is a REBURN, not a bug
  burn_sev_wt       = sum(frac * sev_w),
  burn_decay        = sum(frac * sev_w * exp(-yrs_since / DECAY_YRS)),
  yrs_since_burn    = min(yrs_since),
  max_sev_prior     = max(sev_w),
  n_fire_epochs     = uniqueN(fire_year)
), by = .(cell_id, year, comp_doy)]

# largest single-fire-year burned fraction -- bounded to [0,1], unlike the
# cumulative sum, so it answers "how much of this cell did its biggest fire take"
# NOTE 5000/30 = 166.67, so a cell holds 166 or 167 pixel-centres per side
# depending on how the grid falls against each raster's own origin. Counting
# every centre at a full 900 m2 can therefore overshoot the nominal cell area by
# at most (167/166.67)^2 = 1.004. That 0.4% is pixel discretisation, not double
# counting; it is clamped to 1 here so the column is a true fraction.
fy <- j[, .(fy_frac = sum(frac)), by = .(cell_id, year, comp_doy, fire_year)]
agg <- merge(agg, fy[, .(burn_frac_max = max(fy_frac)),
                     by = .(cell_id, year, comp_doy)],
             by = c("cell_id","year","comp_doy"), all.x = TRUE)
agg[burn_frac_max > 1, burn_frac_max := 1]
rm(fy); invisible(gc())

for (nm in names(EPOCHS)) {
  lo <- EPOCHS[[nm]][1]; hi <- EPOCHS[[nm]][2]
  e <- j[yrs_since >= lo & yrs_since < hi,
         .(v = sum(frac), s = sum(frac * sev_w)), by = .(cell_id, year, comp_doy)]
  setnames(e, c("v","s"), c(paste0("burn_frac_", nm), paste0("burn_sev_", nm)))
  agg <- merge(agg, e, by = c("cell_id","year","comp_doy"), all.x = TRUE)
}
rm(j); invisible(gc())

ct <- merge(ct, agg, by = c("cell_id","year","comp_doy"), all.x = TRUE)
burn_cols <- grep("^burn_", names(ct), value = TRUE)
for (cc in burn_cols) set(ct, which(is.na(ct[[cc]])), cc, 0)
ct[is.na(max_sev_prior), max_sev_prior := 0]
ct[is.na(n_fire_epochs), n_fire_epochs := 0L]
ct[, ever_burned := as.integer(burn_frac_cum > 0)]

# --- neighbourhood burn: the "weighted distance to burn severity" term --------
log("Neighbourhood burn terms (8 adjacent cells)")
nb_src <- unique(ct[, .(cell_id, year, comp_doy, burn_frac_cum, burn_decay)])
nb_src <- merge(nb_src, cells, by = "cell_id")
off <- CJ(dx = -1:1, dy = -1:1)[!(dx == 0 & dy == 0)]
nb <- rbindlist(lapply(seq_len(nrow(off)), function(i) {
  z <- copy(nb_src); z[, `:=`(gx = gx + off$dx[i], gy = gy + off$dy[i])]; z
}))
nb <- merge(nb[, .(gx, gy, year, comp_doy, burn_frac_cum, burn_decay)],
            cells, by = c("gx","gy"))
nb <- nb[, .(burn_frac_nbr = mean(burn_frac_cum), burn_decay_nbr = mean(burn_decay),
             n_nbr = .N), by = .(cell_id, year, comp_doy)]
ct <- merge(ct, nb, by = c("cell_id","year","comp_doy"), all.x = TRUE)
for (cc in c("burn_frac_nbr","burn_decay_nbr")) set(ct, which(is.na(ct[[cc]])), cc, 0)
ct[is.na(n_nbr), n_nbr := 0L]
rm(nb, nb_src); invisible(gc())

# =============================================================================
# 5. attach the vegetation sequence
# =============================================================================
log("Attaching vegetation sequence l0..l", N_LAG, " + previous year")
for (L in 0:N_LAG) {
  v <- veg[, .(cell_id, ci = ci + L, ndvi, evi)]
  setnames(v, c("ndvi","evi"), c(paste0("ndvi_l", L), paste0("evi_l", L)))
  ct <- merge(ct, v, by = c("cell_id","ci"), all.x = TRUE)
}
vp <- veg[, .(cell_id, ci = ci + COMP_PER_YR, ndvi_prevyr = ndvi)]
ct <- merge(ct, vp, by = c("cell_id","ci"), all.x = TRUE)
ct[, `:=`(ndvi_delta_yr = ndvi_l0 - ndvi_prevyr,
          ndvi_trend    = ndvi_l0 - ndvi_l3)]     # change across the 2-month run
rm(veg); invisible(gc())

# =============================================================================
# 6. verification
# =============================================================================
log("Verification")
stopifnot(
  nrow(ct) == uniqueN(ct[, .(cell_id, year, comp_doy)]),   # grain
  all(ct$n_detections <= ct$n_checklists),
  all(ct$any_detection %in% c(0L, 1L)),
  all(ct$any_detection == as.integer(ct$n_detections > 0)),
  all(ct$burn_frac_cum >= 0),                              # >1 == reburn
  all(ct$burn_frac_max >= 0 & ct$burn_frac_max <= 1),       # single year is bounded
  all(is.na(ct$ndvi_l0) | (ct$ndvi_l0 >= -1 & ct$ndvi_l0 <= 1)),
  all(ct$yrs_since_burn[ct$ever_burned == 1] > 0)          # time order held
)
log("  all checks passed")

# =============================================================================
# 7. write
# =============================================================================
col_order <- c(
  "cell_id","gx","gy","year","comp_doy","ci","step_start","date_first","date_last",
  "any_detection","n_detections","n_checklists","det_rate",
  "n_locations","effort_hours","dur_median","dist_median","obs_median","frac_traveling",
  "ever_burned","burn_frac_cum","burn_frac_max","burn_sev_wt","burn_decay","yrs_since_burn",
  "max_sev_prior","n_fire_epochs",
  paste0("burn_frac_", names(EPOCHS)), paste0("burn_sev_", names(EPOCHS)),
  "burn_frac_nbr","burn_decay_nbr","n_nbr",
  paste0("ndvi_l", 0:N_LAG), paste0("evi_l", 0:N_LAG),
  "ndvi_prevyr","ndvi_delta_yr","ndvi_trend")
# ct already carries gx/gy from the original grouping -- merging `cells` back in
# would collide and suffix them
out <- ct[, ..col_order]
setorder(out, cell_id, year, comp_doy)
for (cc in c("burn_frac_cum","burn_frac_max","burn_sev_wt","burn_decay","yrs_since_burn",
             "burn_frac_nbr","burn_decay_nbr",
             paste0("burn_frac_", names(EPOCHS)), paste0("burn_sev_", names(EPOCHS))))
  set(out, NULL, cc, round(out[[cc]], 6))
log("Writing ", basename(out_csv))
fwrite(out, out_csv, na = "", dateTimeAs = "ISO")

# =============================================================================
# 8. summary
# =============================================================================
sink(out_sum)
cat("Wrentit cell x time-step prediction frame\n")
cat("=========================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Script    : regrid_cells.R\n")
cat("Source    : ", basename(in_csv), "\n", sep = "")
cat("Grid      : ", GRID_M / 1000, " km equal-area cells (EPSG:5070)\n", sep = "")
cat("Time step : 16 days, aligned to the MODIS composite calendar\n")
cat("Span      : ", format(min(out$step_start)), " .. ", format(max(out$date_last)), "\n", sep = "")
cat("Runtime   : ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
    " min\n", sep = "")

cat("\n-- FRAME SIZE --------------------------------------------------------\n")
cat("  cell-time units (rows)   : ", fmt(nrow(out)), "\n", sep = "")
cat("  distinct cells           : ", fmt(uniqueN(out$cell_id)), "\n", sep = "")
cat("  distinct time steps      : ", fmt(uniqueN(out$ci)), "\n", sep = "")
cat("  checklists represented   : ", fmt(sum(out$n_checklists)), "\n", sep = "")
cat("  columns                  : ", ncol(out), "\n", sep = "")
cat("  NOTE only cell-times WITH birding effort are emitted. The full N x T\n")
cat("  lattice would be ", fmt(uniqueN(out$cell_id) * uniqueN(out$ci)),
    " cells x steps; ", sprintf("%.1f%%", 100 * nrow(out) /
    (uniqueN(out$cell_id) * uniqueN(out$ci))), " of it was ever visited.\n", sep = "")

cat("\n-- TARGET ------------------------------------------------------------\n")
cat("  any_detection = 1 : ", fmt(sum(out$any_detection)),
    sprintf("  (%.1f%%)\n", 100 * mean(out$any_detection)), sep = "")
cat("  mean det_rate     : ", round(mean(out$det_rate), 4), "\n", sep = "")
cat("  EFFORT WARNING -- any_detection by checklists-per-unit:\n")
print(as.data.frame(out[, .(units = .N, any_det = round(mean(any_detection), 3)),
        by = .(checklists = cut(n_checklists, c(0,1,2,5,10,20,Inf),
               labels = c("1","2","3-5","6-10","11-20","20+")))][order(checklists)]),
      row.names = FALSE)
cat("  This is why det_rate or a fixed-effort subsample is the safer target:\n")
cat("  'any detection' climbs with effort regardless of the birds.\n")

cat("\n-- BURN COVERAGE -----------------------------------------------------\n")
cat("  units with a prior burn in-cell   : ", fmt(sum(out$ever_burned)),
    sprintf("  (%.1f%%)\n", 100 * mean(out$ever_burned)), sep = "")
cat("  units with burn in a NEIGHBOUR    : ", fmt(sum(out$burn_frac_nbr > 0)),
    sprintf("  (%.1f%%)\n", 100 * mean(out$burn_frac_nbr > 0)), sep = "")
cat("  cells that ever burned            : ",
    fmt(uniqueN(out[ever_burned == 1]$cell_id)), " of ", fmt(uniqueN(out$cell_id)), "\n", sep = "")
cat("  REBURNS -- units whose cell saw >1 fire year:\n")
print(as.data.frame(out[, .N, by = n_fire_epochs][order(n_fire_epochs)]), row.names = FALSE)
cat("  (the checklist-grain table could only ever see ONE fire per row)\n")
cat("\n  CUMULATIVE burned fraction, where > 0 (sums across fire years, so a\n")
cat("  value above 1 means the cell reburned):\n")
print(summary(out[burn_frac_cum > 0]$burn_frac_cum))
cat("  units with burn_frac_cum > 1 (reburn): ", fmt(sum(out$burn_frac_cum > 1)), "\n", sep="")
cat("\n  LARGEST SINGLE fire-year burned fraction, where > 0 (clamped to 0-1;\n")
cat("  see the note in the script on 0.4% pixel-discretisation overshoot):\n")
print(summary(out[burn_frac_max > 0]$burn_frac_max))
cat("\n  years since most recent burn:\n"); print(summary(out[ever_burned == 1]$yrs_since_burn))

cat("\n-- DETECTION vs BURN (raw, CONFOUNDED -- see note) --------------------\n")
print(as.data.frame(out[, .(units = .N, any_det = round(mean(any_detection), 3),
                            det_rate = round(mean(det_rate), 4),
                            chk_per_unit = round(mean(n_checklists), 1)),
        by = .(burned = ifelse(ever_burned == 1, "cell burned", "never burned"))]),
      row.names = FALSE)
cat("  NOT a fire effect: fires concentrate in chaparral, which is Wrentit\n")
cat("  habitat. Use within-cell before/after or matched cells to remove it.\n")

cat("\n-- VEGETATION SEQUENCE -----------------------------------------------\n")
for (L in 0:N_LAG) {
  v <- out[[paste0("ndvi_l", L)]]
  cat(sprintf("  ndvi_l%d : median %6.3f | NA %s (%.1f%%)\n", L, median(v, na.rm = TRUE),
              fmt(sum(is.na(v))), 100 * mean(is.na(v))))
}
cat("  ndvi_trend (l0 - l3, 2-month change):\n"); print(summary(out$ndvi_trend))

cat("\n-- MISSING VALUES ----------------------------------------------------\n")
mv <- sapply(out, function(x) sum(is.na(x)))
mv <- mv[mv > 0]
if (length(mv)) print(data.frame(column = names(mv), na = as.integer(mv),
      pct = round(100 * mv / nrow(out), 2), row.names = NULL)) else cat("  none\n")
cat("\n  dist_median NA  : every checklist in the unit was Stationary\n")
cat("  ndvi_l* NA      : MODIS fill (ocean/edge) or the lag predates the archive\n")

cat("\n-- SUGGESTED MODELLING SETUP -----------------------------------------\n")
cat("  target    : any_detection (binary) or det_rate (rate)\n")
cat("  EXCLUDE   : cell_id, gx, gy, ci, step_start, date_first, date_last,\n")
cat("              n_detections, det_rate  <- leak the target if used as features\n")
cat("  effort    : n_checklists, effort_hours, dur_median, dist_median,\n")
cat("              obs_median, frac_traveling, n_locations\n")
cat("  splitting : group by cell_id, or train year <= 2023 / test year >= 2024.\n")
cat("              NEVER split at random -- cells repeat across time steps.\n")
cat("\nCOLUMNS (", ncol(out), ")\n", sep = ""); print(names(out))
sink()
cat(readLines(out_sum), sep = "\n")
log(sprintf("DONE in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
