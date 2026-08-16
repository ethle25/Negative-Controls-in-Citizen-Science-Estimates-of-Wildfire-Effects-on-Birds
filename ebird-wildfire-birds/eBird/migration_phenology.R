# =============================================================================
# GAP 2: does wildfire shift bird MIGRATION timing?
#
# The research question asks about occurrence AND migration. Nothing in the
# project has touched migration. This closes that gap.
#
# WHAT "MIGRATION" CAN MEAN IN eBird DATA
#   Not individual tracking -- these are checklists, not tagged birds. What IS
#   measurable is PHENOLOGY: when in the year a species is detectable in a
#   place. For a migrant that is arrival, departure and season length. If fire
#   changes when a species is present, that shows up as a shift in these.
#
# THE THREE SPECIES ARE A NATURAL DESIGN
#   Olive-sided Flycatcher  long-distance migrant (winters in South America)
#   Brown Creeper           short-distance / altitudinal movements
#   Wrentit                 SEDENTARY -- among the most resident birds in North
#                           America; individuals may spend their whole lives in
#                           a few hundred metres.
#   Wrentit therefore acts as a PLACEBO: it cannot migrate, so any apparent
#   post-fire "timing shift" in Wrentit is detectability or observer behaviour,
#   not migration. It calibrates how much apparent shift is artefact.
#
# MEASURES, per cell-year (16-day resolution, so only shifts of ~2+ weeks are
# detectable):
#   peak_doy   effort-weighted mean day-of-year of detection, using the
#              DETECTION RATE per 16-day bin rather than raw counts, so uneven
#              seasonal birding effort does not drive the answer
#   first_doy  first bin of the year with a detection   (arrival proxy)
#   last_doy   last bin with a detection                (departure proxy)
#   season_len last_doy - first_doy
#
# Then the same within-cell design as everywhere else: does the measure shift
# after the cell's own fire, relative to never-burned controls, with cell and
# year fixed effects and cluster-robust SEs?
#
# WRITES NEW FILES ONLY. Reads existing frames read-only.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
out_sum <- file.path(dir_, "migration_phenology_summary.txt")

# Species and their movement descriptions come from species_registry.csv.
#
# ORDER MATTERS HERE, unlike the other scripts: migrants are printed first and
# the SEDENTARY species last, because the sedentary row is the placebo the reader
# is meant to discount the migrant rows against. Sorting by `migratory`
# descending then common_name puts residents after migrants and is deterministic,
# so adding a species cannot silently reshuffle the report. The placebo is
# identified by migratory == 0, not by naming Wrentit.
sp_reg <- cfg_species(active_only = TRUE)
sp_reg <- sp_reg[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]
if (!nrow(sp_reg)) stop("no active species in species_registry.csv has a built frame")
if (!any(sp_reg$migratory == 0))
  say("WARNING: no sedentary species active -- the phenology placebo is MISSING, ",
      "and any apparent timing shift is uncalibrated")
setorder(sp_reg, -migratory, common_name)
SPP <- setNames(lapply(seq_len(nrow(sp_reg)), function(i)
                  list(lab = sp_reg$common_name[i], mig = sp_reg$movement[i])),
                sp_reg$prefix)

MIN_DET  <- as.integer(Sys.getenv("PHEN_MIN_DET", "5"))   # detections per cell-year
MIN_BINS <- as.integer(Sys.getenv("PHEN_MIN_BINS", "3"))  # distinct bins with detections
DOSE_MIN <- as.numeric(Sys.getenv("PHEN_DOSE_MIN", "0.50"))
DOSE     <- Sys.getenv("PHEN_DOSE", "frac_max")   # any treatment_registry.csv name

twfe_cl <- function(dt, xcols, iters = 40) {
  z <- copy(dt)[, c("cell_id","year","yv","w", xcols), with = FALSE]
  cl <- z$cell_id
  for (cc in c("yv", xcols)) set(z, NULL, cc, as.numeric(z[[cc]]))
  for (i in seq_len(iters)) for (cc in c("yv", xcols)) {
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
  }
  X <- as.matrix(z[, ..xcols]); y <- z$yv; w <- z$w
  ok <- apply(X, 2, function(v) sum(abs(v)) > 1e-12); X <- X[, ok, drop = FALSE]
  XtWX <- t(X * w) %*% X
  if (abs(det(XtWX)) < 1e-26) return(NULL)
  Xi <- solve(XtWX); b <- as.vector(Xi %*% (t(X * w) %*% y))
  e <- as.vector(y - X %*% b); U <- X * (w * e)
  V <- Xi %*% (t(rowsum(U, cl)) %*% rowsum(U, cl)) %*% Xi
  data.table(term = colnames(X), beta = b, se = sqrt(diag(V)))
}

phen <- list(); res <- list(); meta <- list()
for (pre in names(SPP)) {
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  if (!file.exists(f)) { say("missing ", basename(f)); next }
  say("=== ", SPP[[pre]]$lab, " ===")
  d <- fread(f, showProgress = FALSE,
             select = c("cell_id","year","comp_doy","ever_burned","burn_frac_max",
                        "n_checklists","n_detections"))

  # per-cell-year phenology, computed from RATE per bin (effort-normalised)
  d[, rate_bin := n_detections / n_checklists]
  ph <- d[, {
    hit <- n_detections > 0
    if (sum(hit) >= MIN_BINS && sum(n_detections) >= MIN_DET) {
      wsum <- sum(rate_bin)
      # every element must be double: data.table requires consistent column
      # types across groups, and min() on an integer column returns integer
      # while the else-branch returns NA_real_
      .(peak_doy   = as.numeric(sum(rate_bin * comp_doy) / wsum),
        first_doy  = as.numeric(min(comp_doy[hit])),
        last_doy   = as.numeric(max(comp_doy[hit])),
        season_len = as.numeric(max(comp_doy[hit]) - min(comp_doy[hit])),
        ndet = as.numeric(sum(n_detections)), nbin = as.numeric(sum(hit)))
    } else .(peak_doy = NA_real_, first_doy = NA_real_, last_doy = NA_real_,
             season_len = NA_real_, ndet = as.numeric(sum(n_detections)),
             nbin = as.numeric(sum(hit)))
  }, by = .(cell_id, year)]
  ph <- ph[!is.na(peak_doy)]
  say(sprintf("  cell-years with usable phenology: %s (>=%d detections, >=%d bins)",
              format(nrow(ph), big.mark=","), MIN_DET, MIN_BINS))
  if (nrow(ph) < 200) { say("  too few -- skipping this species"); next }

  # Treated/control/dose from the shared definition. cell_treatment() rather than
  # build_panel() because the phenology measures are computed from 16-day bins and
  # build_panel() would have already collapsed those to cell-year.
  ct <- cell_treatment(pre, dose = DOSE)
  if (length(attr(ct, "dose_terms")) > 1L)
    stop("dose '", DOSE, "' declares ", length(attr(ct, "dose_terms")),
         " terms; this script thresholds on a SINGLE dose. Using only the first\n",
         "  would silently threshold on an arbitrary component (e.g. s1 of\n",
         "  class_doses). Use a single-term measure, or analyze()/sweep() which\n",
         "  handle multi-term doses properly.")
  DCOL <- attr(ct, "dose_terms")[1]
  never <- ct[treat == 0]$cell_id
  tr <- ct[treat == 1 & get(DCOL) >= DOSE_MIN, .(cell_id, ev_year)]

  pn <- ph[cell_id %in% c(tr$cell_id, never)]
  pn <- merge(pn, tr[, .(cell_id, ev_year)], by = "cell_id", all.x = TRUE)
  pn[, post := as.numeric(!is.na(ev_year) & year >= ev_year)]
  pn[, w := as.numeric(ndet)]                       # weight by detections
  nt <- uniqueN(pn[!is.na(ev_year)]$cell_id); nc <- uniqueN(pn[is.na(ev_year)]$cell_id)
  say(sprintf("  treated cells with phenology: %d | control cells: %d", nt, nc))
  if (nt < 15) { say("  too few treated cells -- estimates unreliable, still reported") }

  out <- list()
  for (m in c("peak_doy","first_doy","last_doy","season_len")) {
    z <- copy(pn); setnames(z, m, "yv")
    r <- twfe_cl(z, "post")
    if (!is.null(r)) out[[m]] <- data.table(measure = m, beta = r$beta[1], se = r$se[1])
  }
  res[[pre]] <- rbindlist(out)
  phen[[pre]] <- pn
  meta[[pre]] <- data.table(species = pre, lab = SPP[[pre]]$lab, mig = SPP[[pre]]$mig,
                            cellyears = nrow(pn), treated = nt, controls = nc,
                            mean_peak = mean(pn$peak_doy),
                            mean_first = mean(pn$first_doy),
                            mean_last = mean(pn$last_doy))
}
M <- rbindlist(meta)

# =============================================================================
sink(out_sum)
cat("Wildfire and bird MIGRATION timing (phenology)\n")
cat("==============================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Question  : does fire shift WHEN in the year a species is detectable?\n")
cat("Resolution: 16-day composite bins -- only shifts of ~2+ weeks are visible\n")
cat("Treated   : cells that burned >= ", sprintf("%.0f%%", 100*DOSE_MIN), " of their area\n", sep="")
cat("Design    : cell + year fixed effects, cluster-robust SEs, weighted by\n")
cat("            detections; measures computed from per-bin detection RATE so\n")
cat("            seasonal birding effort does not drive them\n")
cat("Filter    : cell-years with >= ", MIN_DET, " detections in >= ", MIN_BINS,
    " distinct bins\n\n", sep = "")

cat("-- SAMPLE AND BASELINE PHENOLOGY -------------------------------------\n")
print(as.data.frame(M[, .(species = lab, movement = mig,
                          cell_years = cellyears, treated = treated, controls = controls,
                          peak_doy = round(mean_peak), first_doy = round(mean_first),
                          last_doy = round(mean_last))]), row.names = FALSE)
cat("\n  Baseline phenology is a sanity check on the measures themselves:\n")
cat("  the long-distance migrant should show a LATER first_doy and a SHORTER\n")
cat("  season than the resident. If it does not, the measures are not capturing\n")
cat("  migration and nothing below should be believed.\n")

cat("\n-- POST-FIRE SHIFT IN TIMING (days) ----------------------------------\n")
cat("  Positive = later in the year. '*' = 95% CI excludes zero.\n")
for (pre in names(SPP)) {
  if (is.null(res[[pre]])) next
  cat("\n  === ", SPP[[pre]]$lab, "  [", SPP[[pre]]$mig, "] ===\n", sep = "")
  r <- res[[pre]]
  cat(sprintf("  %-12s %9s %8s %22s\n", "measure", "shift", "se", "95% CI"))
  for (i in seq_len(nrow(r))) {
    lo <- r$beta[i] - 1.96*r$se[i]; hi <- r$beta[i] + 1.96*r$se[i]
    cat(sprintf("  %-12s %+9.2f %8.2f  [%+.2f, %+.2f] %s\n",
                r$measure[i], r$beta[i], r$se[i], lo, hi,
                if (lo > 0 || hi < 0) "*" else ""))
  }
}

cat("\n-- HOW TO READ -------------------------------------------------------\n")
cat("  first_doy later  + last_doy earlier  -> compressed season\n")
cat("  first_doy earlier                    -> earlier arrival\n")
cat("  last_doy later                       -> delayed departure / longer stay\n")
cat("  peak_doy shift                       -> the season moved as a whole\n")
cat("\n  THE WRENTIT ROW IS THE CONTROL. It does not migrate, so a large shift\n")
cat("  there means the measure is picking up detectability or observer change,\n")
cat("  and the migrant results should be discounted by roughly that much.\n")

cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat("  * 16-day bins: sub-fortnight shifts are invisible. A real 5-day advance\n")
cat("    in arrival would not be detected here.\n")
cat("  * first_doy / last_doy are extremes and therefore sensitive to how much\n")
cat("    birding happened early and late in the year; peak_doy is more robust.\n")
cat("  * Requires enough detections per cell-year, which filters toward\n")
cat("    well-birded cells -- not a random sample of burned landscape.\n")
cat("  * This measures when a species is DETECTED, not when individuals move.\n")
cat("    A shift is consistent with changed phenology, changed abundance during\n")
cat("    part of the season, or changed detectability in altered vegetation.\n")
cat("  * MTBS under-maps 2023-2025; detection is not occupancy.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
