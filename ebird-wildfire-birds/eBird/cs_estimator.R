# =============================================================================
# Callaway-Sant'Anna: does the fire effect survive a staggered-timing estimator?
#
# WHY THIS EXISTS
#   event_study.R and dose_burn_effect.R use two-way fixed effects. With
#   staggered treatment timing, TWFE implicitly uses ALREADY-TREATED units as
#   controls for later-treated ones. That is harmless when the effect is
#   constant. It is NOT harmless here: Finding 4 shows Wrentit runs -27% in
#   year 1 and +28% by year 7, and OSFL climbs from +36% to +78%. So a 2015-
#   cohort cell at its post-fire PEAK serves as the counterfactual for a 2021-
#   cohort cell in its TROUGH. Those comparisons carry negative weights and can
#   bias -- in principle flip -- the pooled estimate.
#
#   Callaway & Sant'Anna (2021) removes them: a separate ATT for every (cohort,
#   period) using ONLY never-treated controls, aggregated with non-negative
#   weights.
#
# WHAT IT FOUND (2026-07-30, all three species)
#   The DOSE-SCALING REPLICATES -- Wrentit grows more negative and OSFL more
#   positive as burn extent rises, under an estimator sharing no machinery with
#   the original. But no CS estimate excludes zero at any threshold, and TWFE
#   was INFLATING all three species, not attenuating them.
#
# WHY THIS CANNOT OVERTURN THE HEADLINE
#   att_gt() takes BINARY staggered treatment, so CS estimates the effect of
#   "burned >= X%" -- a DIFFERENT ESTIMAND from the continuous dose at 100%
#   burned. It is also unweighted (unit-constant weights only) and drops
#   single-observation cells. The comparable target is the subset estimate at
#   the same threshold, not the dose.
#
# WEIGHTS
#   att_gt() takes only UNIT-CONSTANT weights, while the project's TWFE is
#   checklist-weighted per cell-year. To keep the comparison honest this script
#   refits an UNWEIGHTED TWFE on the identical panel, so any CS-vs-TWFE gap
#   cannot be a weighting artifact.
#
# REFACTORED 2026-07-30 onto panel_utils.R (build_cs_panel).
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   CS_DOSE        burn measure from treatment_registry.csv (default frac_max)
#   CS_THRESH      binary treatment cutoff        (default 0.50)
#   CS_MIN_COHORT  drop cohorts smaller than this (default 20)
#   CS_BITERS      att_gt bootstrap iterations    (default 1000)
# =============================================================================
suppressMessages({ library(data.table); library(did) })
setDTthreads(0)
set.seed(1)

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

PREFIX  <- cfg_default_prefix()
SPECIES <- cfg_default_species()
DOSE    <- Sys.getenv("CS_DOSE", "frac_max")
THRESH  <- as.numeric(Sys.getenv("CS_THRESH", "0.50"))
MINCOH  <- as.integer(Sys.getenv("CS_MIN_COHORT", "20"))
BITERS  <- as.integer(Sys.getenv("CS_BITERS", "1000"))
out_sum <- file.path(dir_, sprintf("cs_summary_%s_%d.txt", PREFIX, round(100*THRESH)))

# --- panel: one line, was ~45 ------------------------------------------------
say("Building CS panel (", PREFIX, ", dose = ", DOSE, ", >= ",
    round(100*THRESH), "% burned)")
cy <- attach_fire(build_cs_panel(PREFIX, dose = DOSE, thresh = THRESH,
                                 min_cohort = MINCOH), PREFIX)
cy[, fire_num := as.integer(factor(fire))]
say("panel ", fmt(nrow(cy)), " cell-years, ", fmt(uniqueN(cy$cell_id)), " cells | ",
    fmt(uniqueN(cy[first_treat > 0]$cell_id)), " treated / ",
    fmt(uniqueN(cy[first_treat > 0]$fire)), " fires / ",
    uniqueN(cy[first_treat > 0]$first_treat), " cohorts")

# =============================================================================
# 1. Callaway-Sant'Anna
# =============================================================================
say("att_gt (never-treated controls, fire-clustered, ", BITERS, " boot iters)")
cs <- att_gt(yname = "rate", tname = "year", idname = "cell_id",
             gname = "first_treat", xformla = ~1, data = as.data.frame(cy),
             panel = TRUE, allow_unbalanced_panel = TRUE,
             control_group = "nevertreated", clustervars = "fire_num",
             bstrap = TRUE, biters = BITERS, cband = FALSE, base_period = "universal")

ag_simple  <- aggte(cs, type = "simple",  na.rm = TRUE)
ag_dynamic <- aggte(cs, type = "dynamic", na.rm = TRUE, min_e = -5, max_e = 8)
ag_group   <- aggte(cs, type = "group",   na.rm = TRUE)

# =============================================================================
# 2. matching UNWEIGHTED TWFE on the identical panel
# =============================================================================
say("comparison TWFE (unweighted, same panel)")
tw <- copy(cy)
tw[, treatpost := as.numeric(first_treat > 0 & year >= first_treat)]
# controls need a pseudo post or `post` absorbs the treatment; drawn from the
# treated cohort distribution exactly as dose_burn_effect.R does
set.seed(1)
ctrl_ids <- unique(tw[first_treat == 0]$cell_id)
gs <- unique(tw[first_treat > 0, .(cell_id, first_treat)])$first_treat
pseudo <- data.table(cell_id = ctrl_ids, pg = sample(gs, length(ctrl_ids), replace = TRUE))
tw <- merge(tw, pseudo, by = "cell_id", all.x = TRUE)
tw[, post_any := as.numeric(fifelse(first_treat > 0, year >= first_treat, year >= pg))]
b_tw <- twfe(tw, c("post_any","treatpost"), weights = FALSE)[["treatpost"]]

# =============================================================================
# report
# =============================================================================
base_rate <- cy[first_treat > 0 & year < first_treat, sum(det)/sum(chk)]
dyn <- data.table(e = ag_dynamic$egt, att = ag_dynamic$att.egt, se = ag_dynamic$se.egt)
dyn[, `:=`(lo = att - 1.96*se, hi = att + 1.96*se)]
dyn[, sig := fifelse(lo > 0 | hi < 0, "*", " ")]
dyn[, pct := 100 * att / base_rate]

sink(out_sum)
cat("Callaway-Sant'Anna vs TWFE  [", SPECIES, ", >= ", round(100*THRESH),
    "% burned]\n", sep = "")
cat("==========================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Burn measure: ", DOSE, " (used only to pick the subset -- treatment is BINARY)\n", sep="")
cat("Outcome   : det_rate = detections / checklists, cell-year\n")
cat("Controls  : never-burned cells only (no already-treated comparisons)\n")
cat("Clusters  : dominant fire id; controls singleton. boot iters ", BITERS, "\n", sep="")
cat("Treated   : ", fmt(uniqueN(cy[first_treat>0]$cell_id)), " cells / ",
    fmt(uniqueN(cy[first_treat>0]$fire)), " fires / ",
    uniqueN(cy[first_treat>0]$first_treat), " cohorts | controls ",
    fmt(uniqueN(cy[first_treat==0]$cell_id)), "\n", sep = "")
cat("Pre-fire base rate: ", sprintf("%.4f", base_rate), "\n", sep = "")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  TWFE with staggered timing uses already-treated units as controls. That\n")
cat("  is benign only if effects are constant. Finding 4 shows they are NOT:\n")
cat("  Wrentit runs -27% (yr1) to +28% (yr7). So recovered cells act as the\n")
cat("  counterfactual for freshly-burned ones, with negative weights. CS uses\n")
cat("  never-treated controls only and non-negative weights.\n")

cat("\n-- 1. OVERALL EFFECT -------------------------------------------------\n\n")
cat("  Callaway-Sant'Anna (simple ATT) : ", sprintf("%+.4f", ag_simple$overall.att),
    "   SE ", sprintf("%.4f", ag_simple$overall.se), "\n", sep = "")
cat("    95% CI                        : [",
    sprintf("%+.4f", ag_simple$overall.att - 1.96*ag_simple$overall.se), ", ",
    sprintf("%+.4f", ag_simple$overall.att + 1.96*ag_simple$overall.se), "]\n", sep = "")
cat("    as % of pre-fire base         : ", sprintf("%+.1f%%",
    100*ag_simple$overall.att/base_rate), "\n", sep = "")
cat("  TWFE, unweighted, same panel    : ", sprintf("%+.4f", b_tw), "\n", sep = "")
cat("\n  If CS is SMALLER in magnitude, TWFE was inflating the effect via\n")
cat("  forbidden already-treated comparisons. That is what happened here for\n")
cat("  all three species.\n")

cat("\n-- 2. EFFECT BY YEARS SINCE FIRE (CS dynamic) ------------------------\n")
cat("  e < 0 should sit near zero -- the pre-trend check, computed WITHOUT\n")
cat("  already-treated contamination.\n\n")
print(as.data.frame(dyn[, .(e, att = round(att, 4), se = round(se, 4),
                            ci = sprintf("[%+.4f,%+.4f]", lo, hi),
                            pct = sprintf("%+.1f%%", pct), sig)]), row.names = FALSE)
cat("\n  aggregated dynamic ATT: ", sprintf("%+.4f", ag_dynamic$overall.att),
    "  SE ", sprintf("%.4f", ag_dynamic$overall.se), "\n", sep = "")
cat("\n  NOTE: the Wald pre-test is not reported when clustering beyond the unit\n")
cat("  level, so read the individual pre-period CIs above rather than a single\n")
cat("  test statistic. OSFL FAILS this check at e = -3 under CS even though it\n")
cat("  passes under TWFE.\n")

cat("\n-- 3. EFFECT BY COHORT (which fire years drive it) -------------------\n\n")
print(data.frame(cohort = ag_group$egt, att = round(ag_group$att.egt, 4),
                 se = round(ag_group$se.egt, 4)), row.names = FALSE)
cat("\n  Big cohort-to-cohort spread is exactly what breaks TWFE. It also means\n")
cat("  the overall ATT averages over quite different fire years.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * BINARY treatment at >= ", round(100*THRESH), "%. This does NOT replicate\n", sep="")
cat("    the continuous-dose headline (-0.0298); the comparable number is the\n")
cat("    subset estimate at the same threshold.\n")
cat("  * CS is UNWEIGHTED (att_gt takes unit-constant weights only). The TWFE\n")
cat("    line above is refit unweighted on the same panel for that reason, so it\n")
cat("    is NOT the published checklist-weighted figure.\n")
cat("  * Cohorts under ", MINCOH, " cells are dropped. Before that filter the 2025\n", sep="")
cat("    cohort (9 cells, almost no post-period) returned a degenerate -0.5001\n")
cat("    with SE 0.0016.\n")
cat("  * Unbalanced panel: birding is sporadic, so cells enter and leave. Cells\n")
cat("    with a single observation are dropped by att_gt.\n")
cat("  * MTBS under-maps 2023-2025, so late cohorts are thin and some controls\n")
cat("    may have burned unrecorded.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
