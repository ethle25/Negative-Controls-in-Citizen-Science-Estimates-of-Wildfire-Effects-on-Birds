# =============================================================================
# Dose-response CIs clustered by FIRE instead of by CELL
#
# WHY THIS EXISTS
#   dose_burn_effect.R bootstraps cells. cell_fire_lookup.R shows that is wrong:
#   3,125 burned cells come from 480 fires -- mean 8.0 cells per fire, max 231
#   (the 2021 Dixie fire). Cells inside one perimeter share ONE ignition, one
#   weather system, one fuel type, one set of road closures and one post-fire
#   access regime. They are one treatment draw, not eight. Resampling them
#   independently pretends the data holds ~8x more information about fire than
#   it does.
#
# READ THIS BEFORE QUOTING ANYTHING FROM HERE
#   These are PERCENTILE intervals and they are NOT the project's verdict.
#   wild_cluster_boot.R supersedes them. Percentile intervals at ~225 clusters
#   got the WIDTH right and the POSITION wrong: Wrentit's implied SE was 0.0174
#   against an analytic 0.0171, but the interval centred on -0.0367 rather than
#   the actual -0.0298 because the resampled distribution is skewed, dragging
#   the bound clear of zero. OSFL's interval even flipped verdict on rep count
#   (B=200 fails, B=2000 clears). This script is kept because the width ratio it
#   reports is genuinely informative -- it measures how much false precision
#   cell-level resampling was buying -- not because its verdicts are right.
#
# WHAT CHANGES vs dose_burn_effect.R
#   ONLY the resampling unit. Same panel, same estimator, same point estimates.
#   Stratified cluster bootstrap: treated resamples FIRES (every cell in a drawn
#   fire travels with it, preserving within-fire correlation); controls resample
#   CELLS (never-burned cells share no fire -- see CAVEATS for what that misses).
#   Point estimates are unchanged by construction; the report prints both so it
#   can be checked.
#
# REFACTORED 2026-07-30 onto panel_utils.R + infer.R.
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   FC_BOOT_B      replications                              (default 200)
#   FC_DOSE        burn measure from treatment_registry.csv  (default frac_max)
#   FC_SUBSETS=0   headline dose only, skip the threshold table
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)
set.seed(1)

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
source(file.path(dir_, "infer.R"))

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

PREFIX  <- cfg_default_prefix()
SPECIES <- cfg_default_species()
DOSE    <- Sys.getenv("FC_DOSE", "frac_max")
B_REPS  <- as.integer(Sys.getenv("FC_BOOT_B", "200"))
DO_SUBS <- Sys.getenv("FC_SUBSETS", "1") != "0"
out_sum <- file.path(dir_, sprintf("fire_cluster_summary_%s.txt", PREFIX))

# --- panel: two lines, was ~45 ----------------------------------------------
say("Building panel (", PREFIX, ", dose = ", DOSE, ")")
cy <- attach_fire(build_panel(PREFIX, dose = DOSE), PREFIX)
terms <- attr(cy, "dose_terms")
stopifnot(length(terms) == 1L)          # the threshold table needs a single dose
DCOL <- terms[1]
cy[, treatpost := as.numeric(treat) * as.numeric(post)]
n_unmatched <- uniqueN(cy[treat == 1 & grepl("^CTRL_", fire)]$cell_id)
say("panel ", fmt(nrow(cy)), " | treated cells ", fmt(uniqueN(cy[treat==1]$cell_id)),
    " | fire clusters ", fmt(uniqueN(cy[treat == 1]$fire)))

# boot_cell(), boot_fire() and boot_ci() now live in infer.R -- dose_burn_effect.R
# needs the same cell bootstrap, so keeping a private copy here would recreate
# exactly the duplication this refactor is removing.
ci <- boot_ci

# --- 1. headline: continuous dose -------------------------------------------
XD <- c("post_d", paste0("post_", DCOL))
say("Continuous dose: point estimate")
pt_dose <- twfe(cy, XD)[[2]]
say("  ", sprintf("%+.4f", pt_dose))
say("Continuous dose: cell bootstrap (", B_REPS, ")")
cd <- ci(boot_cell(cy, XD, B_REPS), 2)
say("Continuous dose: FIRE bootstrap (", B_REPS, ")")
fd <- ci(boot_fire(cy, XD, B_REPS), 2)
n_fire_all <- uniqueN(cy[treat == 1]$fire)

# --- 2. heavy-burn subsets ---------------------------------------------------
THR <- c(0, 0.25, 0.50, 0.70, 0.90)
if (DO_SUBS) say("Subsets") else say("Subsets skipped (FC_SUBSETS=0)")
subs <- if (!DO_SUBS) NULL else rbindlist(lapply(THR, function(thr) {
  keep <- cy[treat == 0 | get(DCOL) >= thr]
  XS <- c("post_d", "treatpost")
  est <- twfe(keep, XS)[["treatpost"]]
  say("  >= ", round(100*thr), "%: cell boot");  cc <- ci(boot_cell(keep, XS, B_REPS), "treatpost")
  say("  >= ", round(100*thr), "%: fire boot");  cf <- ci(boot_fire(keep, XS, B_REPS), "treatpost")
  data.table(threshold = sprintf(">= %d%%", round(100*thr)),
             cells = uniqueN(keep[treat == 1]$cell_id),
             fires = uniqueN(keep[treat == 1]$fire), est = est,
             cell_lo = cc[1], cell_hi = cc[2], fire_lo = cf[1], fire_hi = cf[2],
             cell_w = cc[2]-cc[1], fire_w = cf[2]-cf[1])
}))

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Fire-clustered CIs vs cell-clustered  [", SPECIES, "]\n", sep = "")
cat("=====================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Burn measure: ", DOSE, " (", attr(cy, "dose_label"), ")\n", sep = "")
cat("Panel     : shared build_panel(); only the BOOTSTRAP UNIT changes here\n")
cat("Reps      : ", B_REPS, " per interval\n", sep = "")

cat("\n-- THESE ARE NOT THE PROJECT'S VERDICT -------------------------------\n")
cat("  Percentile intervals at ~225 clusters get the WIDTH right and the\n")
cat("  POSITION wrong. wild_cluster_boot.R supersedes every significance call\n")
cat("  below. What IS informative here is the width RATIO -- it measures how\n")
cat("  much false precision cell-level resampling was buying.\n")

cat("\n-- WHY CLUSTER BY FIRE -----------------------------------------------\n")
cat("  3,125 burned cells come from 480 fires. Mean 8.0 cells per fire, max 231\n")
cat("  (2021 Dixie). Cells in one perimeter share one ignition, one weather\n")
cat("  system, one fuel type, one access regime -- one treatment draw, not eight.\n")
cat("\n  Treated cells here: ", fmt(uniqueN(cy[treat==1]$cell_id)), " across ",
    fmt(n_fire_all), " clusters\n", sep = "")
if (n_unmatched) cat("  (unmatched treated cells kept as singletons: ",
                     fmt(n_unmatched), ")\n", sep = "")

cat("\n-- 1. CONTINUOUS DOSE (effect at 100% burned) ------------------------\n\n")
cat("  point estimate            : ", sprintf("%+.4f", pt_dose), "\n", sep = "")
cat("  cell-clustered 95% CI     : [", sprintf("%+.4f", cd[1]), ", ",
    sprintf("%+.4f", cd[2]), "]   width ", sprintf("%.4f", cd[2]-cd[1]), "\n", sep = "")
cat("  FIRE-clustered 95% CI     : [", sprintf("%+.4f", fd[1]), ", ",
    sprintf("%+.4f", fd[2]), "]   width ", sprintf("%.4f", fd[2]-fd[1]), "\n", sep = "")
cat("  width ratio (fire / cell) : ", sprintf("%.2fx", (fd[2]-fd[1])/(cd[2]-cd[1])),
    "   <- the informative number\n", sep = "")

cat("\n-- 2. HEAVY-BURN SUBSETS ---------------------------------------------\n")
if (is.null(subs)) cat("  skipped (FC_SUBSETS=0)\n") else {
cat("  'fires' is the effective sample size for anything fire-related.\n\n")
print(as.data.frame(subs[, .(threshold, cells, fires, est = round(est, 4),
                             cell_ci = sprintf("[%+.4f,%+.4f]", cell_lo, cell_hi),
                             fire_ci = sprintf("[%+.4f,%+.4f]", fire_lo, fire_hi),
                             ratio = sprintf("%.2fx", fire_w/cell_w))]),
      row.names = FALSE) }

cat("\n-- HOW TO READ -------------------------------------------------------\n")
cat("  * Point estimates are IDENTICAL to dose_burn_effect.R by construction.\n")
cat("    If one differs, the panel was built wrong.\n")
cat("  * A wider fire CI is not a worse result -- it is the honest one.\n")
cat("  * For SIGNIFICANCE, read wild_cluster_boot.R, not this file.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * Controls are still resampled as CELLS. Never-burned cells share no\n")
cat("    fire but are spatially autocorrelated, so even these are not fully\n")
cat("    conservative. A spatial block bootstrap on gx/gy would go further.\n")
cat("  * Cluster label is the DOMINANT fire (largest perimeter overlap) at the\n")
cat("    cell's flip year. 2.0% of cell-years are touched by more than one fire.\n")
cat("  * Width ratios below ~1.0 in the subset table are Monte-Carlo noise at\n")
cat("    B=200, not evidence that fire clustering tightens anything.\n")
cat("  * MTBS under-maps 2023-2025, so recent fires are missing here too.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
