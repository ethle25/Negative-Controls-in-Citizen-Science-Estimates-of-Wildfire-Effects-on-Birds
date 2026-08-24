# =============================================================================
# Does the fire effect scale with HOW MUCH of the cell burned?
#
# did_burn_effect.R treated fire as yes/no: a cell that burned 5% counted the
# same as one that burned 95%. Since the median treated cell is only ~29%
# burned, that binary contrast is diluted -- the leading explanation for why it
# found nothing. This script uses burn EXTENT as the treatment instead, which
# does three useful things:
#
#   1. DOSE-RESPONSE. "Does detection fall further where more of the cell
#      burned?" A dose-response is much stronger evidence than a yes/no
#      difference, because confounding would have to scale with burn extent to
#      fake it.
#   2. AN UNDILUTED ESTIMATE, without leaving the 5 km grid. 290 cells burned
#      >= 70% and 189 burned >= 90%; in those, "burned cell" really does mean
#      "burned land".
#   3. A SEVERITY CHECK among heavily burned cells.
#
# Design: every comparison is WITHIN a cell, before vs after its own fire,
# against never-burned controls carrying pseudo event years. Outcome is
# det_rate (effort-normalised).
#
# READ THIS BEFORE QUOTING THE INTERVALS
#   The CIs here are CELL-level percentile bootstraps. Both the unit and the
#   instrument are now known to be wrong: cells inside one fire perimeter are
#   not independent (480 fires behind 3,125 burned cells, max 231 in one fire),
#   and percentile intervals mis-POSITION at ~225 clusters even when their width
#   is right. wild_cluster_boot.R supersedes every significance call below --
#   under it NO species reaches significance (Wrentit p = 0.110, BRCR 0.251,
#   OSFL 0.080). The point estimates here remain the project's headline numbers;
#   the intervals do not.
#
# REFACTORED 2026-07-30 onto panel_utils.R + infer.R. The burn measure is no
# longer hardcoded -- set DOSE_MEASURE to any name in treatment_registry.csv.
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   DOSE_MEASURE    burn measure (default frac_max = largest single fire-year
#                   burned fraction, bounded 0-1; burn_frac_cum is NOT a
#                   fraction and exceeds 1 on reburns)
#   DOSE_SUBSETS=0  skip section 3 (5 x 200 bootstrapped fits, ~45 min)
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
DOSE    <- Sys.getenv("DOSE_MEASURE", "frac_max")
DO_SUBS <- Sys.getenv("DOSE_SUBSETS", "1") != "0"
out_sum <- file.path(dir_, sprintf("burn_dose_summary_%s.txt", PREFIX))

# --- panel: one line, was ~30 ------------------------------------------------
say("Building panel (", PREFIX, ", dose = ", DOSE, ")")
cy <- build_panel(PREFIX, dose = DOSE)
terms <- attr(cy, "dose_terms")
stopifnot(length(terms) == 1L)     # bands and subsets need a single dose
DCOL <- terms[1]; PCOL <- paste0("post_", DCOL)
cy[, treatpost := as.numeric(treat) * as.numeric(post)]
say("treated ", fmt(uniqueN(cy[treat==1]$cell_id)),
    " | control ", fmt(uniqueN(cy[treat==0]$cell_id)),
    " | cell-year panel: ", fmt(nrow(cy)))

# worst severity class ever recorded in the cell -- for section 4 only, never a
# regressor. NB it SATURATES (285 of 290 heavily burned cells return class 4),
# which is exactly why section 4 fails to find a gradient.
sev <- fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX)),
             select = c("cell_id","ever_burned","max_sev_prior"), showProgress = FALSE
       )[ever_burned == 1, .(sev = max(max_sev_prior)), by = cell_id]
cy <- merge(cy, sev, by = "cell_id", all.x = TRUE)
cy[is.na(sev), sev := 0]

# =============================================================================
# 1. before/after change by DOSE BAND
# =============================================================================
say("Dose bands")
cy[, band := cut(get(DCOL), c(-0.001, 0.10, 0.25, 0.50, 0.70, 0.90, 1.001),
     labels = c("1-10%","10-25%","25-50%","50-70%","70-90%","90-100%"))]
bands <- cy[treat == 1, .(cells = uniqueN(cell_id), rate = sum(det)/sum(w)),
            by = .(band, period = ifelse(post == 1, "after", "before"))]
bw <- dcast(bands, band ~ period, value.var = c("rate","cells"))
bw[, change := rate_after - rate_before]
ctrl_b <- cy[treat == 0, .(rate = sum(det)/sum(w)),
             by = .(period = ifelse(post == 1, "after", "before"))]
ctrl_change <- ctrl_b[period == "after"]$rate - ctrl_b[period == "before"]$rate
bw[, did := change - ctrl_change]

# =============================================================================
# 2. CONTINUOUS dose-response
#    rate ~ cell FE + year FE + post + post x dose
#    The dose is time-invariant so its main effect is absorbed by the cell FE;
#    the interaction is the dose-response slope.
# =============================================================================
say("Continuous dose-response")
XD <- c("post_d", PCOL)
b_dose <- twfe(cy, XD)
ci_dose <- apply(boot_cell(cy, XD, B = 200), 2, quantile, c(0.025, 0.975), na.rm = TRUE)

# =============================================================================
# 3. UNDILUTED estimates: heavily burned cells + all controls
#    Regress on treat x post, NOT post alone -- controls carry a pseudo post
#    too, so post alone estimates the common trend and shrinks every threshold
#    toward zero (this silently flattened all of them to ~-0.003 once).
# =============================================================================
say("Heavily-burned subsets")
sub_est <- function(thr) {
  keep <- cy[treat == 0 | get(DCOL) >= thr]
  XS  <- c("post_d","treatpost")
  est <- twfe(keep, XS)[["treatpost"]]
  ci  <- boot_ci(boot_cell(keep, XS, B = 200), "treatpost")
  base <- keep[treat == 1 & post == 0, sum(det)/sum(w)]
  data.table(threshold = sprintf(">= %d%%", round(100*thr)),
             cells = uniqueN(keep[treat == 1]$cell_id),
             base_rate = base, beta = est, lo = ci[1], hi = ci[2],
             pct = 100 * est / base)
}
subs <- if (DO_SUBS) rbindlist(lapply(c(0, 0.25, 0.50, 0.70, 0.90), sub_est)) else NULL

# =============================================================================
# 4. severity check within heavily burned cells
# =============================================================================
say("Severity within heavily burned cells")
hi <- cy[treat == 0 | get(DCOL) >= 0.70]
hi[, sev_band := fifelse(treat == 0, "control",
                 fifelse(sev >= 4, "high (4)",
                 fifelse(sev == 3, "moderate (3)", "low (1-2)")))]
sw <- dcast(hi[treat == 1, .(cells = uniqueN(cell_id), rate = sum(det)/sum(w)),
               by = .(sev_band, period = ifelse(post == 1, "after", "before"))],
            sev_band ~ period, value.var = c("rate","cells"))
sw[, change := rate_after - rate_before]

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Does the fire effect scale with HOW MUCH of the cell burned?  [", SPECIES, "]\n", sep="")
cat("===========================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Outcome   : det_rate = detections / checklists (effort-normalised)\n")
cat("Treatment : ", DOSE, " -- ", attr(cy, "dose_label"), "\n", sep = "")
cat("Design    : within-cell before/after vs never-burned controls, cell + year\n")
cat("            fixed effects, checklist-weighted, cell-level bootstrap\n")

cat("\n-- INTERVALS HERE ARE SUPERSEDED -------------------------------------\n")
cat("  These CIs bootstrap CELLS. Cells inside one fire are not independent\n")
cat("  (480 fires behind 3,125 burned cells, max 231 in one fire), and\n")
cat("  percentile intervals mis-position at ~225 clusters. Under the wild\n")
cat("  cluster bootstrap NO species reaches significance: Wrentit p = 0.110,\n")
cat("  Brown Creeper 0.251, OSFL 0.080. The POINT ESTIMATES below stand; the\n")
cat("  intervals do not. See wild_cluster_boot.R.\n")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  The yes/no analysis found -0.0024 [-0.0066, +0.0020] -- no detectable\n")
cat("  effect. But the median treated cell is only ~29% burned, so that test\n")
cat("  largely compared mostly-unburned land with unburned land. If fire\n")
cat("  matters, the effect should GROW with burn extent. That is testable here.\n")

cat("\n-- 1. BEFORE/AFTER BY BURN EXTENT ------------------------------------\n")
cat("  'did' subtracts the never-burned control trend over the same period.\n\n")
print(as.data.frame(bw[, .(burn_extent = band, cells = cells_after,
                           rate_before = round(rate_before, 4),
                           rate_after = round(rate_after, 4),
                           change = round(change, 4), did = round(did, 4))]),
      row.names = FALSE)
cat("\n  control cells, same period split: ", sprintf("%+.4f", ctrl_change), "\n", sep = "")
cat("  READ THIS COLUMN-WISE: if fire matters, `did` should get MORE NEGATIVE\n")
cat("  as burn extent rises. A flat column means no dose-response.\n")

cat("\n-- 2. CONTINUOUS DOSE-RESPONSE ---------------------------------------\n")
cat("  rate ~ cell FE + year FE + post + post x dose\n\n")
cat("  post (effect at ~0 dose)           : ", sprintf("%+.4f", b_dose[["post_d"]]),
    "   95% CI [", sprintf("%+.4f", ci_dose[1,"post_d"]), ", ",
    sprintf("%+.4f", ci_dose[2,"post_d"]), "]\n", sep = "")
cat("  post x dose (dose-response)        : ", sprintf("%+.4f", b_dose[[PCOL]]),
    "   95% CI [", sprintf("%+.4f", ci_dose[1,PCOL]), ", ",
    sprintf("%+.4f", ci_dose[2,PCOL]), "]\n", sep = "")
cat("  => treatment effect at full dose   : ", sprintf("%+.4f", b_dose[[PCOL]]),
    "   (the interaction ALONE; post_d is the common trend controls share)\n", sep = "")

cat("\n-- 3. UNDILUTED ESTIMATES (restrict to heavily burned cells) ---------\n")
if (is.null(subs)) cat("  skipped (DOSE_SUBSETS=0)\n") else
print(as.data.frame(subs[, .(treated_cells_burned = threshold, cells,
                             base_rate = round(base_rate, 4),
                             effect = round(beta, 4),
                             ci_low = round(lo, 4), ci_high = round(hi, 4),
                             pct_change = round(pct, 1))]), row.names = FALSE)
cat("\n  '>= 0%' is the whole treated group. Moving down trades sample size for a\n")
cat("  cleaner treatment: at >= 90%, 'this cell burned' really does mean 'this\n")
cat("  land burned', so the dilution objection no longer applies.\n")

cat("\n-- 4. SEVERITY WITHIN HEAVILY BURNED CELLS (>= 70% burned) -----------\n")
print(as.data.frame(sw[, .(severity = sev_band, cells = cells_after,
                           rate_before = round(rate_before, 4),
                           rate_after = round(rate_after, 4),
                           change = round(change, 4))]), row.names = FALSE)
cat("\n  This finds nothing because max_sev_prior SATURATES -- it records the\n")
cat("  worst patch, and 285 of 290 heavily burned cells contain SOME class-4.\n")
cat("  For a real severity dose use DOSE_MEASURE=two_dose_sev, which splits\n")
cat("  extent from fraction-at-high-severity.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * Sample size falls as the threshold rises: 1,009 treated cells at >= 0%\n")
cat("    down to 189 at >= 90%. Intervals widen accordingly.\n")
cat("  * Heavily burned cells came from BIG fires, which may sit in more remote\n")
cat("    terrain. The within-cell design controls for that, but the RESULT then\n")
cat("    applies to big-fire landscapes.\n")
cat("  * Burn extent is measured over the 5 km cell, not the exact spots birders\n")
cat("    walked. Point-grain analysis would be sharper -- and effort_placebo.R\n")
cat("    shows people visit 43% fewer locations per cell after fire, which the\n")
cat("    cell fixed effect cannot absorb.\n")
cat("  * Detection is not occupancy.\n")
cat("  * MTBS under-maps 2023-2025, so some controls may have burned unrecorded.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
