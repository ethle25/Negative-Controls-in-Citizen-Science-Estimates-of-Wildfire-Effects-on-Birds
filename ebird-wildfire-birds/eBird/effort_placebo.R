# =============================================================================
# Does fire change HOW PEOPLE BIRD?  (the confound det_rate does not remove)
#
# WHY THIS EXISTS
#   det_rate = detections / checklists divides out HOW MANY checklists a
#   cell-time received. It does NOT divide out how people birded: how far they
#   walked, how long they stayed, how many were in the party, how many distinct
#   spots inside the cell they visited.
#
#   Every one of those plausibly changes after a fire. Roads and trails close.
#   Burned ground is open, so sightlines lengthen. Burn scars draw curious
#   birders. And a checklist covering more ground detects more species for
#   reasons that have nothing to do with how many birds live there.
#
#   THIS IS A PLACEBO. Finding nothing is the good outcome.
#
# IT DID NOT FIND NOTHING (2026-07-30)
#   Distinct locations visited per cell falls 43% (p = 0.0001) and checklists
#   fall 7% (p = 0.015) at 100% burned; hours and distance trend down at
#   p ~ 0.07. Cell fixed effects compare a cell to itself but CANNOT absorb a
#   shift in which spots INSIDE the cell get sampled. If post-fire access
#   concentrates on roadsides and open ground, the habitat being sampled changes
#   and the species estimates are contaminated. Point-grain analysis -- same
#   location before and after -- is the principled fix and has not been run.
#
# DESIGN
#   Identical panel to dose_burn_effect.R, with the OUTCOME swapped from
#   detections to each effort measure in turn. det_rate is included as a
#   reference row. Inference is the wild cluster bootstrap by fire.
#
# WEIGHTS
#   UNWEIGHTED. The detection model weights by checklist count, but when the
#   outcome IS checklist count that is circular. Note the consequence: the
#   det_rate row here (-0.0401, p = 0.0003) is NOT the published weighted figure
#   (-0.0298, p = 0.110). That specification dependence is a real finding and
#   both numbers belong in any write-up.
#
# REFACTORED 2026-07-30 onto panel_utils.R + infer.R.
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   EP_DOSE   burn measure from treatment_registry.csv (default frac_max)
#   EP_B      wild cluster replications (default 9999)
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
DOSE    <- Sys.getenv("EP_DOSE", "frac_max")
B_REPS  <- as.integer(Sys.getenv("EP_B", "9999"))
out_sum <- file.path(dir_, "effort_placebo_summary.txt")

# --- panel: effort measures carried through, then aggregated -----------------
# effort_hours is a TOTAL so it is summed; the per-checklist medians and the
# location/protocol shares are averaged over the 16-day steps in the year.
say("Building panel (", PREFIX, ", dose = ", DOSE, ")")
cy <- attach_fire(
  build_panel(PREFIX, dose = DOSE,
              extra     = c("dur_median","dist_median","obs_median",
                            "frac_traveling","n_locations"),
              extra_sum = "effort_hours"), PREFIX)
terms <- attr(cy, "dose_terms")
stopifnot(length(terms) == 1L)
DCOL <- terms[1]

cy[, `:=`(det_rate  = rate,
          log_chk   = log1p(chk),
          log_hours = log1p(effort_hours),
          dur_med   = dur_median,
          dist_med  = dist_median,
          obs_med   = obs_median,
          frac_trav = frac_traveling,
          n_loc     = n_locations)]
say("panel ", fmt(nrow(cy)), " | treated fire clusters ",
    fmt(uniqueN(cy[treat == 1]$fire)))

XCOLS <- c("post_d", paste0("post_", DCOL)); KI <- 2L

# --- one outcome: fit unweighted, fire-clustered SE, wild cluster p ----------
run_outcome <- function(yname) {
  z <- cy[is.finite(get(yname))]
  b <- twfe(z, XCOLS, yname = yname, weights = FALSE)
  P <- .prep(copy(z)[, rate := get(yname)], XCOLS, weights = FALSE)
  bb <- .fit(P, P$y)
  se <- sqrt(crve(P, as.vector(P$y - P$X %*% bb))[KI, KI])
  data.table(outcome = yname, n = nrow(z),
             mean_pre = cy[treat == 1 & post == 0, mean(get(yname), na.rm = TRUE)],
             est = bb[KI], se = se, t = bb[KI]/se,
             p = wcb_p(P, KI, bb[KI], se, B = B_REPS))
}

OUTS <- c("det_rate","log_chk","log_hours","dur_med","dist_med","obs_med",
          "frac_trav","n_loc")
res <- rbindlist(lapply(OUTS, function(o) { say("  outcome: ", o); run_outcome(o) }))
res[, pct := 100 * est / mean_pre]

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Does fire change HOW PEOPLE BIRD?  (placebo test)  [", SPECIES, "]\n", sep="")
cat("=================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Burn measure: ", DOSE, " (", attr(cy, "dose_label"), ")\n", sep = "")
cat("Model     : outcome ~ cell FE + year FE + post + post x burn dose\n")
cat("Inference : wild cluster bootstrap by FIRE, B = ", B_REPS, ", null imposed\n", sep="")
cat("Weights   : UNWEIGHTED (weighting by checklists is circular when the\n")
cat("            outcome IS checklists)\n")
cat("Clusters  : ", fmt(uniqueN(cy[treat==1]$fire)), " fires + singleton controls\n", sep="")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  det_rate divides out HOW MANY checklists, not HOW people birded. Roads\n")
cat("  close after fires; burned ground has longer sightlines; burn scars draw\n")
cat("  visitors. A checklist covering more ground finds more species for reasons\n")
cat("  unrelated to bird abundance. If burn dose shifts any effort measure\n")
cat("  within-cell, every species result is contaminated.\n")
cat("\n  THIS IS A PLACEBO. Finding NOTHING is the good outcome.\n")

cat("\n-- RESULT: effect of going from 0% to 100% burned ---------------------\n")
cat("  'pct' is the effect as a share of the pre-fire treated mean.\n")
cat("  log_chk / log_hours are log1p, so est is roughly a proportional change.\n\n")
print(as.data.frame(res[, .(outcome, n, pre_mean = round(mean_pre, 4),
                            est = round(est, 4), se = round(se, 4),
                            t = round(t, 2), p = round(p, 4), pct = round(pct, 1),
                            flag = fifelse(p < 0.05, "*** SHIFTS ***", ""))]),
      row.names = FALSE)

cat("\n-- HOW TO READ -------------------------------------------------------\n")
cat("  * Any effort row flagged at p < 0.05 is a live confound. det_rate cannot\n")
cat("    absorb it. n_loc is the dangerous one: cell FE compares a cell to\n")
cat("    itself but cannot absorb a change in WHICH SPOTS INSIDE IT are visited.\n")
cat("  * dist_med matters too -- distance walked scales almost mechanically with\n")
cat("    species detected per checklist.\n")
cat("  * The det_rate row is UNWEIGHTED and so differs from the published\n")
cat("    weighted figure. That the verdict changes with weighting is itself a\n")
cat("    finding; report both, do not pick the favourable one.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * NULL IS NOT PROOF. With ", fmt(uniqueN(cy[treat==1]$fire)),
    " fire clusters this test has the same\n", sep = "")
cat("    limited power as the species analyses. Report point estimates too.\n")
cat("  * Effort measures are aggregated over the 16-day steps within each\n")
cat("    cell-year, so within-year timing shifts are invisible here.\n")
cat("  * frac_trav and n_loc are compositional: a change could reflect which\n")
cat("    locations inside the cell were visited, not how anyone birded.\n")
cat("  * Controls are singleton clusters, so spatial correlation among\n")
cat("    never-burned cells is unmodelled.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
