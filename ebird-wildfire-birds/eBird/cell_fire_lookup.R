# =============================================================================
# Which FIRE burned each cell?  (the id the cell frame throws away)
#
# WHY THIS EXISTS
#   Every CI in this project bootstraps or clusters by CELL. But 2,910 treated
#   cells come from only ~340 distinct fires with eBird coverage -- ~8.6 cells
#   per fire. Cells inside one perimeter share a SINGLE treatment draw: one
#   ignition, one weather system, one fuel type, one set of road closures. They
#   are not independent observations, so cell-level resampling counts the same
#   event ~9 times over and reports intervals that are too narrow.
#
#   regrid_cells.R already reads Event_ID out of the perimeter
#   shapefile (line ~302) and then DISCARDS it, keeping only cell_id +
#   fire_year + last ignition date. This script recovers it.
#
# WHAT IT PRODUCES
#   cell_fire_lookup.csv : cell_id, fire_year, fire_event_id, overlap_frac,
#                          n_fires_yr
#   One row per (cell, fire year). `fire_event_id` is the DOMINANT fire -- the
#   one covering the most of that cell that year -- because a cell clipped by
#   two fires needs a single cluster label and the bigger overlap is the better
#   claim on it. `n_fires_yr` records how often that choice was even contested.
#
# DESIGN NOTES
#   * The grid is NOT rebuilt from the checklists. cell_id is assigned inside
#     regrid_cells.R by row order of unique (gx,gy), which cannot be
#     reproduced from outside. Instead gx/gy/cell_id are read straight out of
#     the written cell frame, so the mapping is exact by construction.
#   * Cells are 5 km squares on the EPSG:5070 lattice: x0 = gx * 5000. Same
#     arithmetic as the regrid, so the polygons are identical.
#   * Only cells that are ever_burned somewhere in the frame are considered.
#   * Only perimeters igniting 2014+ are read: a DiD-treated cell must flip from
#     unburned to burned INSIDE the 2015-2026 window, so its causing fire lit in
#     that window. Older fires cannot define treatment here. Widen MIN_IG_YEAR
#     if you ever need the full 1984+ history.
# =============================================================================
suppressMessages({ library(data.table); library(sf) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))   # cfg_default_prefix()
PREFIX <- cfg_default_prefix()
GRID_M <- 5000L                       # must match regrid_cells.R
MIN_IG_YEAR <- as.integer(Sys.getenv("MIN_IG_YEAR", "2014"))
CELL_AREA <- as.numeric(GRID_M)^2

perim_shp <- file.path(dir_, "mtbs_perim", "mtbs_perims_DD.shp")
in_csv    <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX))
out_csv   <- file.path(dir_, "cell_fire_lookup.csv")
out_sum   <- file.path(dir_, "cell_fire_lookup_summary.txt")
stopifnot(file.exists(perim_shp), file.exists(in_csv))

# --- cells, straight from the written frame ----------------------------------
say("Reading cell grid from ", basename(in_csv))
cf <- fread(in_csv, select = c("cell_id","gx","gy","ever_burned"), showProgress = FALSE)
burned <- cf[, .(any_burn = max(ever_burned)), by = .(cell_id, gx, gy)][any_burn == 1]
say("cells in frame: ", fmt(uniqueN(cf$cell_id)), " | ever-burned: ", fmt(nrow(burned)))

geom <- st_sfc(lapply(seq_len(nrow(burned)), function(i) {
  x0 <- burned$gx[i] * GRID_M; y0 <- burned$gy[i] * GRID_M
  st_polygon(list(cbind(c(x0, x0 + GRID_M, x0 + GRID_M, x0, x0),
                        c(y0, y0, y0 + GRID_M, y0 + GRID_M, y0))))
}), crs = 5070)
cellpoly <- st_sf(cell_id = burned$cell_id, geometry = geom)

# --- perimeters --------------------------------------------------------------
say("Reading CA perimeters, ignition ", MIN_IG_YEAR, "+")
perim <- st_read(perim_shp, quiet = TRUE, query =
  "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
perim$Ig_Date <- as.Date(perim$Ig_Date)
perim <- perim[!is.na(perim$Ig_Date) &
               as.integer(format(perim$Ig_Date, "%Y")) >= MIN_IG_YEAR, ]
perim <- st_make_valid(perim)
p5070 <- st_transform(perim, 5070)
p5070$fire_year <- as.integer(format(p5070$Ig_Date, "%Y"))
say("perimeters: ", fmt(nrow(p5070)), " across ", uniqueN(p5070$fire_year), " fire years")

# --- overlap area, year by year ----------------------------------------------
# sf-on-sf st_intersection returns ONE ROW PER INTERSECTING PAIR carrying the
# attributes of both inputs, so cell_id and Event_ID come back already aligned.
# Do NOT hand-pair via st_intersects + geometry subscripts: st_intersection has
# no `by_feature` argument, so passing one is silently swallowed by `...` and
# the returned areas do not line up with the pair list (this bug produced a
# minimum n_fires_yr of 8 and cells "touched by" 216 fires in one year).
# Doing it per year keeps peak memory small; sf indexes internally.
say("Intersecting")
acc <- list()
for (yr in sort(unique(p5070$fire_year))) {
  py <- p5070[p5070$fire_year == yr, "Event_ID"]
  inter <- suppressWarnings(st_intersection(cellpoly, py))
  if (!nrow(inter)) next
  acc[[length(acc) + 1L]] <- data.table(
    cell_id       = inter$cell_id,
    fire_year     = yr,
    fire_event_id = inter$Event_ID,
    overlap_m2    = as.numeric(st_area(inter)))
  say("  ", yr, ": ", fmt(nrow(inter)), " cell-fire pairs")
}
ov <- rbindlist(acc)
ov <- ov[overlap_m2 > 0]
# a fire can appear as several polygon features; sum them before ranking
ov <- ov[, .(overlap_m2 = sum(overlap_m2)), by = .(cell_id, fire_year, fire_event_id)]
say("cell-fire overlaps: ", fmt(nrow(ov)))

# --- dominant fire per cell-year --------------------------------------------
ov[, n_fires_yr := .N, by = .(cell_id, fire_year)]
setorder(ov, cell_id, fire_year, -overlap_m2)
dom <- ov[, .SD[1], by = .(cell_id, fire_year)]
dom[, overlap_frac := overlap_m2 / CELL_AREA]
dom <- dom[, .(cell_id, fire_year, fire_event_id,
               overlap_frac = round(overlap_frac, 6), n_fires_yr)]

# gx/gy ARE REQUIRED OUTPUT, not decoration.
# attach_fire() and .perim_attrs() join this lookup on (gx, gy) because cell_id
# is not stable across separately built species. This script wrote only cell_id;
# the gx/gy columns present in the file on disk had been added out-of-band on
# 2026-07-31 when attach_fire() was fixed, and the script was never updated to
# match. So regenerating the lookup produced a file attach_fire() rejects
# outright -- "cell_fire_lookup.csv has no gx/gy columns" -- and every
# fire-clustered analysis stopped. Found 2026-07-31 by re-running it.
dom <- merge(dom, burned[, .(cell_id, gx, gy)], by = "cell_id")
setcolorder(dom, c("cell_id", "gx", "gy"))
setorder(dom, cell_id, fire_year)
stopifnot(all(c("gx","gy") %in% names(dom)))
fwrite(dom, out_csv)
say("wrote ", basename(out_csv), ": ", fmt(nrow(dom)), " rows")

# =============================================================================
# report
# =============================================================================
cells_per_fire <- dom[, .(cells = uniqueN(cell_id)), by = fire_event_id]
by_year <- dom[, .(cells = uniqueN(cell_id), fires = uniqueN(fire_event_id)),
               by = fire_year][order(fire_year)]
by_year[, cells_per_fire := round(cells / fires, 1)]

sink(out_sum)
cat("Cell -> dominant fire lookup\n")
cat("============================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Grid      : ", GRID_M/1000, " km cells, EPSG:5070, read from ",
    basename(in_csv), "\n", sep = "")
cat("Perimeters: ", basename(perim_shp), ", CA only, ignition ", MIN_IG_YEAR,
    "+\n", sep = "")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  All bootstraps in this project resample CELLS. Cells inside one fire\n")
cat("  perimeter share one ignition, one weather system, one fuel type -- they\n")
cat("  are not independent draws. Resampling them separately overstates how\n")
cat("  much information the data holds, so every CI is too narrow. This lookup\n")
cat("  supplies the cluster label needed to fix that.\n")

cat("\n-- 1. COVERAGE -------------------------------------------------------\n")
cat("  ever-burned cells in frame : ", fmt(nrow(burned)), "\n", sep = "")
cat("  cells matched to a fire    : ", fmt(uniqueN(dom$cell_id)), "\n", sep = "")
cat("  distinct fires             : ", fmt(uniqueN(dom$fire_event_id)), "\n", sep = "")
cat("  cell-year rows             : ", fmt(nrow(dom)), "\n", sep = "")

cat("\n-- 2. CELLS PER FIRE (the whole point) -------------------------------\n")
cat("  mean   : ", sprintf("%.1f", mean(cells_per_fire$cells)), "\n", sep = "")
cat("  median : ", sprintf("%.0f", median(cells_per_fire$cells)), "\n", sep = "")
cat("  max    : ", fmt(max(cells_per_fire$cells)), "\n", sep = "")
cat("\n  If this is much above 1, cell-level resampling is double-counting.\n")
cat("\n  largest fires by cells covered:\n")
print(as.data.frame(head(cells_per_fire[order(-cells)], 10)), row.names = FALSE)

cat("\n-- 3. BY FIRE YEAR ---------------------------------------------------\n\n")
print(as.data.frame(by_year), row.names = FALSE)

cat("\n-- 4. CONTESTED CELLS ------------------------------------------------\n")
cat("  cell-years touched by >1 fire: ", fmt(nrow(dom[n_fires_yr > 1])),
    " of ", fmt(nrow(dom)), " (",
    sprintf("%.1f", 100*nrow(dom[n_fires_yr > 1])/nrow(dom)), "%)\n", sep = "")
cat("  These take the fire with the largest overlap. Where the split is near\n")
cat("  50/50 the cluster label is arbitrary, but it only has to be CONSISTENT\n")
cat("  to fix the double-counting.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * overlap_frac is perimeter overlap, NOT burned fraction. MTBS\n")
cat("    perimeters include unburned islands; burn_frac_max in the cell frame\n")
cat("    comes from the severity raster and is the right dose. Do not swap them.\n")
cat("  * Only ignitions ", MIN_IG_YEAR, "+ are read, which covers any cell that\n", sep="")
cat("    flips from unburned to burned inside the 2015-2026 study window.\n")
cat("  * MTBS under-maps 2023-2025, so recent fires are missing entirely here\n")
cat("    just as they are everywhere else in the project.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
