# =============================================================================
# STEP 0 -- is there any distance information to work with?
#
# The join currently marks a checklist as fire-affected only if the point sits
# INSIDE a burn perimeter (verified: 100% of 206,379 matched rows have
# distance_to_perimeter_km <= 0). So a checklist 500 m from a 100,000-acre scar
# is encoded identically to one 200 km away.
#
# Before building anything, answer one question with a number: how many past
# fires actually sit near a birding location, at several radii? If even at 25 km
# the answer is "about one", a tau = 10 nearest-fires feature has nothing to fill
# it with and the idea dies here for the cost of this script.
#
# READ-ONLY. Writes one summary file; touches no analysis input or output.
# =============================================================================
suppressMessages({ library(data.table); library(sf) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))
PREFIX  <- cfg_default_prefix()
RADII_KM <- as.numeric(strsplit(Sys.getenv("R_KM", "5,10,25,50"), ",")[[1]])
out_sum  <- file.path(dir_, "fire_distance_feasibility.txt")

# --- birding locations --------------------------------------------------------
# Unique (lat, lon) from the CLEAN table -- the rows analysis actually uses.
# Locations repeat heavily (50% of checklist rows sit on 0.48% of locations), so
# the distance work is per LOCATION, not per row: distance is time-invariant.
say("Reading locations")
d <- fread(file.path(dir_, sprintf("%s_2015_2026_joined_clean.csv", PREFIX)),
           select = c("latitude","longitude"), showProgress = FALSE)
locs <- unique(d[!is.na(latitude) & !is.na(longitude)])
say("  checklist rows ", fmt(nrow(d)), " -> unique locations ", fmt(nrow(locs)))
rm(d); invisible(gc())

pts <- st_transform(st_as_sf(locs, coords = c("longitude","latitude"), crs = 4326), 5070)

# --- cell centroids (what a cell-grain version would use) ---------------------
say("Reading cell centroids")
cf <- unique(fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX)),
                   select = c("cell_id","gx","gy"), showProgress = FALSE))
cpt <- st_as_sf(data.table(x = cf$gx * 5000 + 2500, y = cf$gy * 5000 + 2500),
                coords = c("x","y"), crs = 5070)
say("  cells ", fmt(nrow(cf)))

# --- MTBS California perimeters ----------------------------------------------
say("Reading MTBS perimeters")
p <- st_read(file.path(dir_, "mtbs_perim", "mtbs_perims_DD.shp"), quiet = TRUE,
             query = "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
p$ig_year <- as.integer(format(as.Date(p$Ig_Date), "%Y"))
p <- p[!is.na(p$ig_year), ]
p <- st_transform(st_make_valid(p), 5070)
say("  CA perimeters ", fmt(nrow(p)), " (", min(p$ig_year), "-", max(p$ig_year), ")")

# In-window fires are the ones that can act as TREATMENT: a cell must be
# observed both before and after its fire, so only 2015+ ignitions create the
# before/after contrast the causal design needs.
p_win <- p[p$ig_year >= 2015, ]
say("  of which ignited 2015+ (usable as treatment): ", fmt(nrow(p_win)))

# --- counts within each radius ------------------------------------------------
count_within <- function(x, poly, r_m) {
  lengths(st_is_within_distance(x, poly, dist = r_m))
}
tab <- function(x, poly, label) {
  rbindlist(lapply(RADII_KM, function(rk) {
    n <- count_within(x, poly, rk * 1000)
    data.table(unit = label, radius_km = rk,
               mean_fires = mean(n), median_fires = median(n),
               pct_with_0 = 100 * mean(n == 0), pct_with_ge2 = 100 * mean(n >= 2),
               pct_with_ge5 = 100 * mean(n >= 5), pct_with_ge10 = 100 * mean(n >= 10),
               max_fires = max(n))
  }))
}

say("Counting fires within radius -- locations x all fires")
res <- tab(pts, p, "location / all fires 1984+")
say("Counting -- locations x in-window fires (2015+)")
res <- rbind(res, tab(pts, p_win, "location / fires 2015+"))
say("Counting -- cell centroids x in-window fires (2015+)")
res <- rbind(res, tab(cpt, p_win, "cell centroid / fires 2015+"))

# --- how much is NEW information? --------------------------------------------
# Today a location only registers a fire if it is INSIDE the perimeter. Anything
# at positive distance is invisible. So: what share of locations currently have
# nothing, but would gain a signal at each radius?
say("Measuring how much of this is information the join does not already have")
inside <- lengths(st_intersects(pts, p_win)) > 0
gain <- rbindlist(lapply(RADII_KM, function(rk) {
  n <- count_within(pts, p_win, rk * 1000)
  data.table(radius_km = rk,
             pct_inside_today   = 100 * mean(inside),
             pct_near_not_inside = 100 * mean(n > 0 & !inside),
             pct_still_nothing  = 100 * mean(n == 0))
}))

# =============================================================================
sink(out_sum)
cat("STEP 0 -- is there distance information to work with?\n")
cat("=====================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, " (locations are species-independent)\n", sep = "")
cat("Locations : ", fmt(nrow(locs)), " unique birding sites\n", sep = "")
cat("Cells     : ", fmt(nrow(cf)), " 5 km cells\n", sep = "")
cat("Fires     : ", fmt(nrow(p)), " CA perimeters, ", fmt(nrow(p_win)),
    " ignited 2015+\n\n", sep = "")

cat("-- THE QUESTION -----------------------------------------------------\n")
cat("  Francisco's proposal takes the tau = 10 most recent fires near a box.\n")
cat("  If a typical site has ~1 fire nearby even at a wide radius, tau = 10 has\n")
cat("  nothing to fill it with. These counts decide that.\n\n")

cat("-- FIRES WITHIN RADIUS ----------------------------------------------\n\n")
print(as.data.frame(res[, .(unit, radius_km,
                            mean = round(mean_fires, 2), median = median_fires,
                            `%none` = round(pct_with_0, 1),
                            `%>=2` = round(pct_with_ge2, 1),
                            `%>=5` = round(pct_with_ge5, 1),
                            `%>=10` = round(pct_with_ge10, 1),
                            max = max_fires)]), row.names = FALSE)

cat("\n-- HOW MUCH IS GENUINELY NEW ----------------------------------------\n")
cat("  'inside today'  = the join already sees this fire (point in perimeter)\n")
cat("  'near, not in'  = INVISIBLE today; this is what the proposal would add\n\n")
print(as.data.frame(gain[, .(radius_km,
                             `%inside today` = round(pct_inside_today, 2),
                             `%near not inside` = round(pct_near_not_inside, 2),
                             `%still nothing` = round(pct_still_nothing, 2))]),
      row.names = FALSE)

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  * '%>=10' is the share of sites that could actually support tau = 10.\n")
cat("    If it is near zero at every radius, cap tau at what the data holds.\n")
cat("  * '%near not inside' is the payoff: sites that currently register NO\n")
cat("    fire but would gain one. If small, distance adds little.\n")
cat("  * Wider radius is not free -- a fire 50 km away is unlikely to change\n")
cat("    whether a bird is seen, so a large count at 50 km is not automatically\n")
cat("    useful signal. Judge 5-10 km hardest.\n")
cat("  * NONE OF THIS ADDS FIRES. 225 independent fires remain the effective\n")
cat("    sample size for anything fire-related; distance can sharpen a measure,\n")
cat("    not manufacture power.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE -> ", basename(out_sum))
