# =============================================================================
# Event study: how does the fire effect evolve, year by year, per species?
#
# Everything so far collapses "after the fire" into ONE number. That throws away
# the two things we most need:
#
#   VALIDITY. Are the pre-fire years flat? The Wrentit result currently depends on
#     the year-effect correction (its raw before/after is positive while the fixed
#     -effects estimate is negative), and "trust the correction" is a weak defence.
#     If treated and control cells were ALREADY diverging before anything burned,
#     the comparison was never fair. This is a check the result can FAIL.
#
#   DURATION. Does the Olive-sided Flycatcher's boost fade as burnt snags rot and
#     fall (typically 5-15 years)? Does Wrentit recover as chaparral regrows
#     (roughly a decade)? The 11.5-year span exists precisely to answer this, and
#     a single before/after average discards it.
#
# SPECIFICATION
#   rate_it = a_i + b_t + sum_k delta_k * D_it^k + e_it
#   D_it^k = 1 if cell i is treated and (t - its fire year) == k, for k in
#   -5..+8, with k = -1 omitted as the reference. NEVER-BURNED cells carry all
#   D = 0 and identify the year effects b_t. delta_k is therefore the difference
#   from the year immediately before the fire, net of the statewide trend.
#
#   Treated cells are restricted to those that burned >= DOSE_MIN of their area,
#   so "treated" means the land actually burned (the dilution fix from
#   dose_burn_effect.R carried into an event-time design).
#
# INFERENCE
#   Cluster-robust (by cell) analytic standard errors on the within-transformed
#   model -- cells, not cell-years, are the independent unit. This replaces the
#   bootstrap used elsewhere: with 13 event-time coefficients a bootstrap would
#   cost ~30 min per species for no extra rigour.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
out_sum <- file.path(dir_, "event_study_summary.txt")

source(file.path(dir_, "panel_utils.R"))
DOSE_MIN <- as.numeric(Sys.getenv("ES_DOSE_MIN", "0.50"))   # >= 50% of cell burned
DOSE     <- Sys.getenv("ES_DOSE", "frac_max")   # any name in treatment_registry.csv
K_LO <- -5L; K_HI <- 8L; K_REF <- -1L

# Species and their habitat labels come from species_registry.csv, not a literal
# list. The hardcoded version silently excluded any species added after it was
# written -- California Thrasher was built end-to-end on 2026-07-31 and never
# appeared here. Only species with a BUILT frame are kept, so a registry row
# without data is skipped rather than erroring.
sp_reg <- cfg_species(active_only = TRUE)
sp_reg <- sp_reg[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]
if (!nrow(sp_reg)) stop("no active species in species_registry.csv has a built frame")
SPP <- setNames(as.list(sprintf("%s (%s)", sp_reg$common_name, sp_reg$habitat)),
                sp_reg$prefix)

# --- within estimator with cluster-robust SEs ---------------------------------
twfe_cl <- function(dt, xcols, iters = 40) {
  z <- copy(dt)[, c("cell_id","year","rate","w", xcols), with = FALSE]
  cl <- z$cell_id                                  # keep clusters before demeaning
  for (cc in c("rate", xcols)) set(z, NULL, cc, as.numeric(z[[cc]]))
  for (i in seq_len(iters)) {                      # two-way within transform
    for (cc in c("rate", xcols)) {
      z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
      z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
    }
  }
  X <- as.matrix(z[, ..xcols]); y <- z$rate; w <- z$w
  keep <- apply(X, 2, function(v) sum(abs(v)) > 1e-12)   # drop empty event bins
  X <- X[, keep, drop = FALSE]
  XtWX <- t(X * w) %*% X
  if (abs(det(XtWX)) < 1e-24) return(NULL)
  Xi <- solve(XtWX)
  b  <- as.vector(Xi %*% (t(X * w) %*% y))
  e  <- as.vector(y - X %*% b)
  U  <- X * (w * e)                                # score contributions
  Ug <- rowsum(U, cl)                              # sum within cluster
  V  <- Xi %*% (t(Ug) %*% Ug) %*% Xi               # cluster-robust sandwich
  se <- sqrt(diag(V))
  data.table(term = colnames(X)[keep][seq_along(b)], beta = b, se = se)
}

results <- list(); meta <- list()
for (pre in names(SPP)) {
  say("=== ", SPP[[pre]], " ===")
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  if (!file.exists(f)) { say("  missing ", basename(f), " -- skipped"); next }
  d <- fread(f, showProgress = FALSE,
             select = c("cell_id","year","ever_burned","n_checklists","n_detections"))

  # Treated/control/dose comes from the shared definition so this agrees with
  # dose_burn_effect.R by construction, and so DOSE_MEASURE reaches this script.
  # Controls keep ev_year = NA here (not the pseudo year build_panel assigns) --
  # an event study needs them carrying all-zero event dummies, not a fake k.
  ct <- cell_treatment(pre, dose = DOSE)
  if (length(attr(ct, "dose_terms")) > 1L)
    stop("dose '", DOSE, "' declares ", length(attr(ct, "dose_terms")),
         " terms; this script thresholds on a SINGLE dose. Using only the first\n",
         "  would silently threshold on an arbitrary component (e.g. s1 of\n",
         "  class_doses). Use a single-term measure, or analyze()/sweep() which\n",
         "  handle multi-term doses properly.")
  DCOL <- attr(ct, "dose_terms")[1]
  treated <- ct[treat == 1 & get(DCOL) >= DOSE_MIN, .(cell_id, ev_year)]
  never   <- ct[treat == 0]$cell_id
  say(sprintf("  treated (>= %.0f%% burned): %d | never-burned controls: %d",
              100 * DOSE_MIN, nrow(treated), length(never)))

  # cell-year panel: treated (with event time) + never-burned controls
  cy <- rbind(
    merge(d[cell_id %in% treated$cell_id], treated, by = "cell_id"
          )[, .(det = sum(n_detections), chk = sum(n_checklists)),
            by = .(cell_id, year, ev_year)][, treat := 1L],
    d[cell_id %in% never, .(det = sum(n_detections), chk = sum(n_checklists)),
      by = .(cell_id, year)][, `:=`(ev_year = NA_integer_, treat = 0L)],
    use.names = TRUE)
  cy <- cy[chk > 0]
  cy[, `:=`(rate = det / chk, w = as.numeric(chk))]
  cy[, k := fifelse(treat == 1L, pmax(K_LO, pmin(K_HI, year - ev_year)), NA_integer_)]

  # event-time dummies; controls all zero; k = -1 is the reference
  ks <- setdiff(K_LO:K_HI, K_REF)
  for (kk in ks) {
    nm <- sprintf("k%s%d", ifelse(kk < 0, "m", "p"), abs(kk))
    cy[, (nm) := as.numeric(!is.na(k) & k == kk)]
  }
  xc <- sprintf("k%s%d", ifelse(ks < 0, "m", "p"), abs(ks))

  r <- twfe_cl(cy, xc)
  if (is.null(r)) { say("  singular -- skipped"); next }
  r[, k := as.integer(sub("p", "", sub("m", "-", sub("^k", "", term))))]
  r[, `:=`(lo = beta - 1.96 * se, hi = beta + 1.96 * se, species = pre)]
  base <- cy[treat == 1 & k %in% c(-3, -2, -1), sum(det) / sum(chk)]
  results[[pre]] <- r
  meta[[pre]] <- data.table(species = pre, label = SPP[[pre]],
                            treated = nrow(treated), controls = length(never),
                            base_rate = base)
  say(sprintf("  done (pre-fire baseline rate %.4f)", base))
}
R <- rbindlist(results); M <- rbindlist(meta)

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Event study: fire effect by years since the cell's own fire\n")
cat("===========================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Treated   : cells that burned >= ", sprintf("%.0f%%", 100 * DOSE_MIN),
    " of their area, in-window\n", sep = "")
cat("Controls  : never-burned cells (carry all event dummies = 0)\n")
cat("Model     : rate ~ cell FE + year FE + event-time dummies (k = -1 omitted)\n")
cat("Inference : cluster-robust by cell; 95% CI = beta +/- 1.96 SE\n")
cat("Outcome   : det_rate = detections / checklists (effort-normalised)\n\n")
print(as.data.frame(M[, .(species = label, treated_cells = treated,
                          control_cells = controls,
                          pre_fire_rate = round(base_rate, 4))]), row.names = FALSE)

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  k < 0  : years BEFORE the fire. These should sit near ZERO. If they\n")
cat("           drift, treated and control cells were already diverging and the\n")
cat("           design is not credible.\n")
cat("  k = 0  : the fire year itself (partly pre-, partly post-fire).\n")
cat("  k > 0  : years AFTER. The shape here is the recovery/decay curve.\n")
cat("  '*' marks coefficients whose 95% CI excludes zero.\n")

for (pre in names(SPP)) {
  if (is.null(results[[pre]])) next
  r <- results[[pre]][order(k)]
  b <- M[species == pre]$base_rate
  cat("\n=== ", SPP[[pre]], " ===\n", sep = "")
  cat(sprintf("%4s %10s %9s %20s %8s %s\n",
              "k", "beta", "se", "95% CI", "%base", ""))
  for (i in seq_len(nrow(r))) {
    sig <- if (r$lo[i] > 0 || r$hi[i] < 0) "*" else ""
    bar_n <- min(30, round(abs(r$beta[i]) / 0.004))
    bar <- if (r$beta[i] < 0) paste0(strrep("<", bar_n)) else paste0(strrep(">", bar_n))
    cat(sprintf("%4d %+10.4f %9.4f  [%+.4f, %+.4f] %+7.1f%% %s%s\n",
                r$k[i], r$beta[i], r$se[i], r$lo[i], r$hi[i],
                100 * r$beta[i] / b, sig, bar))
  }
  pre_k <- r[k < 0]
  cat(sprintf("  PRE-FIRE CHECK: %d of %d pre-fire coefficients differ from zero",
              sum(pre_k$lo > 0 | pre_k$hi < 0), nrow(pre_k)))
  cat(sprintf("  (max |beta| = %.4f)\n", max(abs(pre_k$beta))))
  post <- r[k >= 0]
  cat(sprintf("  POST-FIRE: mean beta years 0-2 = %+.4f | years 3-5 = %+.4f | years 6-8 = %+.4f\n",
              mean(post[k <= 2]$beta),
              if (nrow(post[k %in% 3:5])) mean(post[k %in% 3:5]$beta) else NA_real_,
              if (nrow(post[k %in% 6:8])) mean(post[k %in% 6:8]$beta) else NA_real_))
}

cat("\n-- SIDE BY SIDE (beta, as % of each species' own pre-fire rate) ------\n")
w <- dcast(merge(R, M[, .(species, base_rate)], by = "species")[
             , .(species, k, pct = round(100 * beta / base_rate, 1))],
           k ~ species, value.var = "pct")
print(as.data.frame(w), row.names = FALSE)
cat("\n  Negative = fewer detections than the year before the fire, net of the\n")
cat("  statewide trend. If the two habitat-loss species and the snag specialist\n")
cat("  diverge AFTER k=0 but not before, that is the pattern real fire response\n")
cat("  would produce.\n")

cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat("  * Later k have fewer cells: a cell burning in 2023 has no k=+8. Estimates\n")
cat("    at the far right are thin and their CIs widen.\n")
cat("  * k=0 straddles the fire, so it is not a clean post-fire year.\n")
cat("  * Two-way FE with staggered timing can be biased when effects vary by\n")
cat("    cohort; a Callaway-Sant'Anna estimator would be the robust follow-up.\n")
cat("  * MTBS under-maps 2023-2025, so some controls may have burned unrecorded.\n")
cat("  * Detection is not occupancy.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
