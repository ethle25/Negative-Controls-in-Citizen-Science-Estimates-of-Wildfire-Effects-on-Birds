# =============================================================================
# STEP 1 GATE -- does distance-based treatment lift the usable FIRE count?
#
# THE QUESTION THIS ANSWERS
#   The project's binding constraint is not the burn measure, it is 225: the
#   number of independent fires behind the 1,009 cells in the difference-in-
#   differences. Cells inside one perimeter share one ignition, one weather
#   system, one access regime -- one treatment draw, not eight. Reaching the
#   conventional t = 1.96 for Wrentit needs roughly 285 fires.
#
#   Today a cell is treated only if a fire BURNED it. If "treated" instead means
#   "a fire happened within R km", many more fires enter. Step 0 showed all 1,956
#   CA fires lie within 10 km of some birded cell. The question is how many
#   survive the DiD's actual requirements: the fire must ignite inside the study
#   window, and the cell must be observed both before and after it.
#
#   That is computable WITHOUT changing any shared code, which is why it runs
#   first: it decides whether Step 1 is worth building at all.
#
# READ-ONLY. Writes one report and nothing else.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
PREFIX  <- cfg_default_prefix()
out_sum <- file.path(dir_, "step1_gate_summary.txt")
RADII   <- c(0, 5, 10, 25, 50)

fd <- file.path(dir_, "cell_fire_distance.csv")
if (!file.exists(fd)) stop("cell_fire_distance.csv not found -- run fire_distance_lookup.R first")

say("Reading distance lookup")
D <- fread(fd, showProgress = FALSE)

say("Reading cell frame (", PREFIX, ")")
f <- fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX)),
           select = c("cell_id","gx","gy","year","ever_burned"), showProgress = FALSE)

# --- the CURRENT baseline, straight from the shared panel code ----------------
say("Current treatment (cells that actually burned)")
p_now  <- attach_fire(build_panel(PREFIX, dose = "frac_max"), PREFIX)
cur_cells <- uniqueN(p_now[treat == 1]$cell_id)
cur_fires <- uniqueN(p_now[treat == 1]$fire)
say("  ", fmt(cur_cells), " cells / ", fmt(cur_fires), " fires")

# --- observation window per cell ---------------------------------------------
# A cell can only contribute to a before/after comparison if it was birded both
# before and after the fire. This is the filter that stops the count being a
# meaningless "every fire is near something".
obs <- f[, .(y_min = min(year), y_max = max(year)), by = .(cell_id, gx, gy)]
YR_LO <- min(f$year); YR_HI <- max(f$year)
say("study window: ", YR_LO, "-", YR_HI)

D <- merge(D, obs[, .(gx, gy, y_min, y_max)], by = c("gx","gy"))   # join on geography

res <- rbindlist(lapply(RADII, function(R) {
  s <- D[dist_km <= R & fire_year > y_min & fire_year <= y_max & fire_year >= YR_LO]
  data.table(radius_km = R,
             cells = uniqueN(s$cell_id),
             fires = uniqueN(s$fire_event_id),
             pairs = nrow(s),
             median_fires_per_cell = if (nrow(s)) as.numeric(median(s[, .N, by = cell_id]$N)) else 0)
}))
res[, vs_current := sprintf("%+.0f%%", 100 * (fires - cur_fires) / cur_fires)]
res[, clears_285 := fifelse(fires >= 285, "YES", "no")]

sink(out_sum)
cat("STEP 1 GATE -- does distance-based treatment add independent FIRES?\n")
cat("==================================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, " (geometry is species-independent; any built frame gives the same answer)\n", sep = "")
cat("Window    : ", YR_LO, "-", YR_HI, "\n\n", sep = "")

cat("-- WHY THIS NUMBER IS THE ONE THAT MATTERS --------------------------\n")
cat("  225 fires is the effective sample size for everything fire-related in\n")
cat("  this project -- not 4.4M checklists, not 622k cell-times, not 1,009\n")
cat("  cells. Wrentit needs roughly 285 to reach t = 1.96. If a distance-based\n")
cat("  treatment cannot clear that, the project stays power-limited no matter\n")
cat("  how the burn measure is engineered.\n\n")

cat(sprintf("  CURRENT (cell actually burned) : %s cells / %s FIRES\n\n",
            fmt(cur_cells), fmt(cur_fires)))

cat("-- IF 'TREATED' MEANT 'A FIRE WITHIN R km' --------------------------\n")
cat("   Counting only fires that ignite inside the window AND after the cell\n")
cat("   was first birded -- i.e. fires that could support a before/after test.\n\n")
print(as.data.frame(res), row.names = FALSE)

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  radius 0 km  = fires whose perimeter contains the cell centre. This is\n")
cat("                 the closest thing to the current definition and is the\n")
cat("                 sanity check: it should land near ", cur_fires, ".\n", sep = "")
cat("  Larger radii trade RELEVANCE for POWER. A fire 25 km away plausibly does\n")
cat("  very little, so those are weakly-treated units: more fires, weaker signal\n")
cat("  per fire. This table does not say which radius is right -- it says how\n")
cat("  much power is on the table at each.\n")

cat("\n-- VERDICT ----------------------------------------------------------\n")
r10 <- res[radius_km == 10]; r25 <- res[radius_km == 25]
cat(sprintf("  at 10 km : %s fires (%s vs current)  -- clears 285: %s\n",
            fmt(r10$fires), r10$vs_current, r10$clears_285))
cat(sprintf("  at 25 km : %s fires (%s vs current)  -- clears 285: %s\n",
            fmt(r25$fires), r25$vs_current, r25$clears_285))
if (max(res$fires) >= 285) {
  cat("\n  => THE POWER CEILING IS NOT FIXED AT 225. A distance-based treatment\n")
  cat("     reaches the range where the Wrentit result could resolve. Step 1 is\n")
  cat("     worth building.\n")
} else {
  cat("\n  => EVEN WITH DISTANCE THE FIRE COUNT STAYS BELOW 285. Distance may\n")
  cat("     still sharpen the measure, but it will not fix the power problem.\n")
}

cat("\n-- WHAT THIS DOES *NOT* SHOW ----------------------------------------\n")
cat("  * That distant fires actually affect birds. More fires only helps if the\n")
cat("    added units carry real signal; otherwise it is noise with a bigger n.\n")
cat("  * Anything about effect size or direction. That needs the fit, not this.\n")
cat("  * Independence. Two fires in the same year and county are not two clean\n")
cat("    draws, so the effective count is somewhere below these numbers.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
