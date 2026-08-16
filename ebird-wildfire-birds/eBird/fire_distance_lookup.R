# =============================================================================
# Cell -> ALL nearby fires, with distance. The variable the project does not have.
#
# WHY THIS EXISTS
#   cell_fire_lookup.R records only fires that OVERLAP a cell. Verified at
#   checklist grain: 100% of the 206,379 fire-matched rows sit INSIDE a
#   perimeter, so a point 500 m from a 100,000-acre scar is encoded identically
#   to one 200 km away. Everything outward-facing is missing.
#
#   Step 0 measured what that costs: at 25 km the median birding location has
#   NINE prior fires, and 92% of locations that currently carry no fire
#   information at all have one nearby. All 1,956 CA fires lie within 10 km of
#   some birded cell, against the 225 fires behind the current DiD.
#
# WHAT IT WRITES
#   cell_fire_distance.csv -- one row per (cell, fire) pair within MAX_KM:
#     cell_id, gx, gy, fire_event_id, fire_year, ig_date, dist_km
#
#   dist_km is 0 for a fire whose perimeter contains the cell centre, and > 0
#   for one the cell is outside of. THAT SECOND CASE IS THE NEW INFORMATION.
#
#   gx/gy are written deliberately: cell_id is not stable across separately
#   built species, and cell_fire_lookup.R shipped for weeks without the gx/gy
#   its own consumers required. Join this on (gx, gy).
#
# TIME-INVARIANCE
#   Distance never changes, so it is computed ONCE for the whole grid and reused.
#   Recency is applied downstream, where the observation date is known.
#
# ENV
#   FD_MAX_KM   widest radius to record (default 50)
# =============================================================================
suppressMessages({ library(data.table); library(sf) })
setDTthreads(0)
t0 <- Sys.time()
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))
PREFIX <- cfg_default_prefix()
MAX_KM <- as.numeric(Sys.getenv("FD_MAX_KM", "50"))
GRID_M <- 5000L                      # must match regrid_cells.R

in_csv  <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX))
out_csv <- file.path(dir_, "cell_fire_distance.csv")
out_sum <- file.path(dir_, "cell_fire_distance_summary.txt")
stopifnot(file.exists(in_csv))

# --- cells: centre of each 5 km square ---------------------------------------
say("Reading cell grid from ", basename(in_csv))
cf <- unique(fread(in_csv, select = c("cell_id","gx","gy"), showProgress = FALSE))
say("  cells: ", fmt(nrow(cf)))
ctr <- st_sfc(lapply(seq_len(nrow(cf)), function(i)
  st_point(c(cf$gx[i] * GRID_M + GRID_M/2, cf$gy[i] * GRID_M + GRID_M/2))), crs = 5070)

# --- perimeters ---------------------------------------------------------------
say("Reading CA perimeters")
p <- st_read(file.path(dir_, "mtbs_perim", "mtbs_perims_DD.shp"), quiet = TRUE,
             query = "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
p$Ig_Date <- as.Date(p$Ig_Date)
p <- p[!is.na(p$Ig_Date), ]
p <- st_make_valid(st_transform(p, 5070))
p$fire_year <- as.integer(format(p$Ig_Date, "%Y"))
say("  perimeters: ", fmt(nrow(p)), " (", min(p$fire_year), "-", max(p$fire_year), ")")

# --- candidate pairs ----------------------------------------------------------
say("Spatial query at ", MAX_KM, " km")
t1 <- Sys.time()
hit <- st_is_within_distance(ctr, p, dist = MAX_KM * 1000)
pair <- data.table(ci = rep(seq_along(hit), lengths(hit)), pi = unlist(hit))
say("  ", fmt(nrow(pair)), " candidate pairs in ",
    round(as.numeric(difftime(Sys.time(), t1, units = "secs")), 1), " s")
if (!nrow(pair)) stop("no cell-fire pairs found -- check the CRS of the perimeters")

# --- exact distances, chunked so peak memory stays flat -----------------------
say("Exact distances (chunked)")
t2 <- Sys.time()
CH <- 200000L
pair[, dist_km := NA_real_]
for (a in seq(1L, nrow(pair), by = CH)) {
  b <- min(a + CH - 1L, nrow(pair))
  pair[a:b, dist_km := as.numeric(
    st_distance(ctr[ci], p[pi, ], by_element = TRUE)) / 1000]
  say("    ", fmt(b), " / ", fmt(nrow(pair)))
}
say("  done in ", round(as.numeric(difftime(Sys.time(), t2, units = "mins")), 1), " min")

out <- data.table(cell_id = cf$cell_id[pair$ci], gx = cf$gx[pair$ci], gy = cf$gy[pair$ci],
                  fire_event_id = p$Event_ID[pair$pi], fire_year = p$fire_year[pair$pi],
                  ig_date = as.character(p$Ig_Date[pair$pi]),
                  dist_km = round(pair$dist_km, 4))
setorder(out, cell_id, dist_km)
stopifnot(all(c("gx","gy") %in% names(out)))     # the cell_fire_lookup.R lesson
fwrite(out, out_csv)
say("wrote ", basename(out_csv), ": ", fmt(nrow(out)), " rows")

# --- report -------------------------------------------------------------------
n_by <- out[, .N, by = cell_id]
sink(out_sum)
cat("Cell -> nearby fires, with distance\n===================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Grid from : ", basename(in_csv), "\n", sep = "")
cat("Radius    : ", MAX_KM, " km\n\n", sep = "")
cat(sprintf("cells            : %s\n", fmt(nrow(cf))))
cat(sprintf("CA perimeters    : %s\n", fmt(nrow(p))))
cat(sprintf("cell-fire pairs  : %s\n", fmt(nrow(out))))
cat(sprintf("distinct fires   : %s of %s CA fires\n",
            fmt(uniqueN(out$fire_event_id)), fmt(nrow(p))))
cat("\n-- fires per cell, by radius ----------------------------------------\n\n")
for (R in c(5, 10, 25, 50)) {
  if (R > MAX_KM) next
  s <- out[dist_km <= R]; n <- s[, .N, by = cell_id]
  cat(sprintf("  %2d km : %5.1f%% of cells | median %2.0f | mean %5.2f | max %3d | %s distinct fires\n",
      R, 100*uniqueN(s$cell_id)/nrow(cf), median(n$N), mean(n$N), max(n$N),
      fmt(uniqueN(s$fire_event_id))))
}
cat("\n-- how much is OUTWARD-facing (the new information) ------------------\n")
cat(sprintf("  pairs with dist_km == 0 (cell centre inside perimeter) : %s\n",
            fmt(sum(out$dist_km == 0))))
cat(sprintf("  pairs with dist_km >  0 (currently unrecorded)         : %s\n",
            fmt(sum(out$dist_km > 0))))
cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat("  * Distance is from the CELL CENTRE, so it is accurate to about half a\n")
cat("    cell (2.5 km). Point grain would be sharper.\n")
cat("  * A fire being near does not mean it affected the cell. This file is the\n")
cat("    raw material for a treatment definition, not a treatment.\n")
cat("  * MTBS under-maps 2023-2025, so recent fires are missing here too.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
