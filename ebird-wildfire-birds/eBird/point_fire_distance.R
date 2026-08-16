# =============================================================================
# BIRDING SPOT -> nearby fires, with distance. Step 2's prerequisite.
#
# WHY SNAPPING
#   There are 498,566 distinct (lat, lon) pairs in the checklist table -- NOT the
#   36,519 CLAUDE.md claims, which is distinct named localities. Running all of
#   them against 1,956 perimeters was still in the spatial query after 11 minutes
#   and the exact-distance pass would have been far worse.
#
#   Most of those coordinates are near-duplicate personal locations a few metres
#   apart. Snapping to a SNAP_M grid collapses them. At 500 m the loss is
#   negligible against a signal measured in kilometres -- and it is 10x finer
#   than the 5 km cells everything else uses.
#
# WHAT IT WRITES
#   point_fire_distance.csv -- one row per (snapped spot, fire) within MAX_KM:
#     px, py (snapped grid indices, EPSG:5070), fire_event_id, fire_year,
#     ig_date, dist_km
#   point_grid_map.csv -- latitude, longitude -> px, py, so checklists can be
#     joined back on. Written separately so the big table stays small.
#
# Distance is time-invariant: computed once, recency applied downstream.
#
# ENV
#   PFD_SNAP_M   snap grid in metres (default 500)
#   PFD_MAX_KM   widest radius      (default 50)
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
SNAP_M <- as.numeric(Sys.getenv("PFD_SNAP_M", "500"))
MAX_KM <- as.numeric(Sys.getenv("PFD_MAX_KM", "50"))

in_csv   <- file.path(dir_, sprintf("%s_2015_2026_joined.csv", PREFIX))
out_csv  <- file.path(dir_, "point_fire_distance.csv")
out_map  <- file.path(dir_, "point_grid_map.csv")
out_sum  <- file.path(dir_, "point_fire_distance_summary.txt")
stopifnot(file.exists(in_csv))

say("Reading checklist coordinates")
d <- unique(fread(in_csv, select = c("latitude","longitude"), showProgress = FALSE))
n_raw <- nrow(d)
say("  distinct coordinates: ", fmt(n_raw))

say("Projecting and snapping to ", SNAP_M, " m")
xy <- st_coordinates(st_transform(
        st_as_sf(d, coords = c("longitude","latitude"), crs = 4326), 5070))
d[, `:=`(px = as.integer(floor(xy[,1] / SNAP_M)),
         py = as.integer(floor(xy[,2] / SNAP_M)))]
fwrite(d, out_map)
g <- unique(d[, .(px, py)])[, pid := .I]
say("  snapped spots: ", fmt(nrow(g)), sprintf("  (%.1fx fewer)", n_raw / nrow(g)))

ctr <- st_sfc(lapply(seq_len(nrow(g)), function(i)
  st_point(c(g$px[i] * SNAP_M + SNAP_M/2, g$py[i] * SNAP_M + SNAP_M/2))), crs = 5070)

say("Reading CA perimeters")
p <- st_read(file.path(dir_, "mtbs_perim", "mtbs_perims_DD.shp"), quiet = TRUE,
             query = "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
p$Ig_Date <- as.Date(p$Ig_Date); p <- p[!is.na(p$Ig_Date), ]
p <- st_make_valid(st_transform(p, 5070))
p$fire_year <- as.integer(format(p$Ig_Date, "%Y"))
say("  perimeters: ", fmt(nrow(p)))

say("Spatial query at ", MAX_KM, " km")
t1 <- Sys.time()
hit <- st_is_within_distance(ctr, p, dist = MAX_KM * 1000)
pair <- data.table(ci = rep(seq_along(hit), lengths(hit)), pi = unlist(hit))
say("  ", fmt(nrow(pair)), " candidate pairs in ",
    round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), " min")

say("Exact distances (chunked)")
t2 <- Sys.time()
CH <- 200000L
pair[, dist_km := NA_real_]
for (a in seq(1L, nrow(pair), by = CH)) {
  b <- min(a + CH - 1L, nrow(pair))
  pair[a:b, dist_km := as.numeric(st_distance(ctr[ci], p[pi, ], by_element = TRUE)) / 1000]
  if (((a - 1L) %/% CH) %% 5 == 0) say("    ", fmt(b), " / ", fmt(nrow(pair)))
}
say("  done in ", round(as.numeric(difftime(Sys.time(), t2, units = "mins")), 1), " min")

out <- data.table(px = g$px[pair$ci], py = g$py[pair$ci],
                  fire_event_id = p$Event_ID[pair$pi], fire_year = p$fire_year[pair$pi],
                  ig_date = as.character(p$Ig_Date[pair$pi]),
                  dist_km = round(pair$dist_km, 4))
setorder(out, px, py, dist_km)
fwrite(out, out_csv)
say("wrote ", basename(out_csv), ": ", fmt(nrow(out)), " rows")

sink(out_sum)
cat("Birding spot -> nearby fires, with distance\n")
cat("===========================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Snap grid : ", SNAP_M, " m | radius ", MAX_KM, " km\n\n", sep = "")
cat(sprintf("distinct coordinates in the data : %s\n", fmt(n_raw)))
cat(sprintf("after snapping to %g m           : %s (%.1fx fewer)\n", SNAP_M, fmt(nrow(g)), n_raw/nrow(g)))
cat(sprintf("CA perimeters                    : %s\n", fmt(nrow(p))))
cat(sprintf("spot-fire pairs                  : %s\n", fmt(nrow(out))))
cat(sprintf("distinct fires reached           : %s\n", fmt(uniqueN(out$fire_event_id))))
cat("\n-- fires per spot, by radius ----------------------------------------\n\n")
for (R in c(5, 10, 25, 50)) {
  if (R > MAX_KM) next
  s <- out[dist_km <= R]; n <- s[, .N, by = .(px, py)]
  cat(sprintf("  %2d km : %5.1f%% of spots | median %2.0f | mean %5.2f | max %3d\n",
      R, 100*nrow(n)/nrow(g), median(n$N), mean(n$N), max(n$N)))
}
cat("\n-- inside vs outside ------------------------------------------------\n")
cat(sprintf("  dist_km == 0 (spot inside a perimeter)  : %s pairs\n", fmt(sum(out$dist_km == 0))))
cat(sprintf("  dist_km >  0 (currently unrecorded)     : %s pairs\n", fmt(sum(out$dist_km > 0))))
cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat(sprintf("  * Snapped to %g m, so distances carry that much slop. Immaterial\n", SNAP_M))
cat("    against a signal measured in kilometres, and 10x finer than the 5 km grid.\n")
cat("  * Join checklists via point_grid_map.csv (latitude, longitude -> px, py).\n")
cat("  * Nearness is not exposure. This is raw material for a treatment\n")
cat("    definition, not a treatment.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
