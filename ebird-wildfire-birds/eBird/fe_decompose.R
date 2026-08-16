# =============================================================================
# WHICH adjustment flips the sign -- cell FE, year FE, or both?
#
# WHY THIS EXISTS
#   burn_dose_summary.txt reports two numbers that disagree in SIGN for the
#   >= 90% burned cells:
#
#     raw pooled before/after, minus control trend : +0.0197
#     cell FE + year FE model                      : -0.0248
#
#   CLAUDE.md used to attribute the flip to year effects ("eBird detection rose
#   steeply over 2015-2026, so 'rose less than it should have' becomes
#   negative"). That was WRONG. det_rate -- the outcome this model actually uses
#   -- is nearly FLAT by year (never-burned cells: 0.0789 in 2015 -> 0.0836 in
#   2026). The steep climb belongs to any_detection, a different column driven
#   by checklist volume.
#
#   This script settles it by fitting the SAME estimator four ways:
#
#     none  : rate ~ 1 + treat + post + treat:post      <- textbook 2x2 DiD
#     cell  : rate ~ cellFE + post + treat:post         <- each cell its own baseline
#     year  : rate ~ yearFE + treat + post + treat:post <- each year its own baseline
#     both  : rate ~ cellFE + yearFE + post + treat:post <- what the project reports
#
#   ANSWER (all three species): it is the CELL fixed effect, every time. Year
#   effects move essentially nothing. That makes the result MORE defensible, not
#   less -- cell FE IS the within-cell design, the same fix Finding 1 rests on,
#   not a trend adjustment bolted on afterwards. The raw pooled number is a
#   BETWEEN-cell contrast carrying exactly the habitat confounding the project
#   rejects; it is not a competing estimate.
#
# REFACTORED 2026-07-30 onto panel_utils.R. The FE-mode switch now lives in the
# shared twfe(); this script only chooses regressors and reports.
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   FE_DOSE        burn measure from treatment_registry.csv (default frac_max)
#   FE_BOOT=0      point estimates only (skip the >= 90% bootstrap)
#   FE_BOOT_B      bootstrap replications (default 200)
# =============================================================================
suppressMessages({ library(data.table) })
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
DOSE    <- Sys.getenv("FE_DOSE", "frac_max")
DO_BOOT <- Sys.getenv("FE_BOOT", "1") != "0"
B_REPS  <- as.integer(Sys.getenv("FE_BOOT_B", "200"))
out_sum <- file.path(dir_, sprintf("fe_decompose_summary_%s.txt", PREFIX))

# --- panel: one line, was ~40 ------------------------------------------------
say("Building panel (", PREFIX, ", dose = ", DOSE, ")")
cy <- build_panel(PREFIX, dose = DOSE)
terms <- attr(cy, "dose_terms")
stopifnot(length(terms) == 1L)          # the threshold table needs a single dose
DCOL <- terms[1]
cy[, `:=`(one = 1, treat_d = as.numeric(treat),
          treatpost = as.numeric(treat) * as.numeric(post))]
say("panel ", fmt(nrow(cy)), " | treated ", fmt(uniqueN(cy[treat==1]$cell_id)),
    " | control ", fmt(uniqueN(cy[treat==0]$cell_id)))

MODES <- c("none","cell","year","both")

# regressors per FE mode: drop whatever the absorbed FE makes collinear.
# cell FE absorbs treat and the dose (both cell-constant); year FE absorbs the
# intercept.
xset <- function(mode, kind) {
  if (kind == "subset")
    switch(mode, none = c("one","treat_d","post_d","treatpost"),
                 cell = c("post_d","treatpost"),
                 year = c("treat_d","post_d","treatpost"),
                 both = c("post_d","treatpost"))
  else
    switch(mode, none = c("one", DCOL, "post_d", paste0("post_", DCOL)),
                 cell = c("post_d", paste0("post_", DCOL)),
                 year = c(DCOL, "post_d", paste0("post_", DCOL)),
                 both = c("post_d", paste0("post_", DCOL)))
}
KEY <- c(subset = "treatpost", dose = paste0("post_", DCOL))

# --- 1. raw pooled before/after, for reference -------------------------------
raw_did <- function(thr) {
  tb <- cy[treat == 1 & get(DCOL) >= thr]
  (tb[post == 1, sum(det)/sum(w)] - tb[post == 0, sum(det)/sum(w)]) -
  (cy[treat == 0 & post == 1, sum(det)/sum(w)] - cy[treat == 0 & post == 0, sum(det)/sum(w)])
}

# --- 2. point estimates: every FE mode x every burn threshold ----------------
say("Point estimates by FE mode")
THR <- c(0, 0.25, 0.50, 0.70, 0.90)
pts <- rbindlist(lapply(THR, function(thr) {
  keep <- cy[treat == 0 | get(DCOL) >= thr]
  r <- data.table(threshold = sprintf(">= %d%%", round(100*thr)),
                  cells = uniqueN(keep[treat == 1]$cell_id),
                  base = keep[treat == 1 & post == 0, sum(det)/sum(w)],
                  raw = raw_did(thr))
  for (m in MODES) r[[m]] <- twfe(keep, xset(m, "subset"), mode = m)[[KEY[["subset"]]]]
  r
}))

say("Point estimates, continuous dose")
dose_pts <- data.table(spec = "continuous dose (effect at 100% burned)")
for (m in MODES) dose_pts[[m]] <- twfe(cy, xset(m, "dose"), mode = m)[[KEY[["dose"]]]]

# --- 3. bootstrap CIs at >= 90%, where the sign conflict lives ---------------
boots <- NULL
if (DO_BOOT) {
  say("Bootstrapping >= 90% subset, ", B_REPS, " reps x 4 FE modes")
  k90 <- cy[treat == 0 | get(DCOL) >= 0.90]
  boot_mode <- function(m, B) {
    xc <- xset(m, "subset"); cells <- unique(k90$cell_id)
    out <- numeric(B)
    for (b in seq_len(B)) {
      s <- data.table(cell_id = sample(cells, length(cells), replace = TRUE))[, newid := .I]
      bb <- merge(s, k90, by = "cell_id", allow.cartesian = TRUE)[, cell_id := newid]
      out[b] <- tryCatch(twfe(bb, xc, iters = 15, mode = m)[["treatpost"]],
                         error = function(e) NA_real_)
    }
    out
  }
  boots <- rbindlist(lapply(MODES, function(m) {
    say("  mode: ", m); bs <- boot_mode(m, B_REPS)
    q <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)
    data.table(mode = m, est = twfe(k90, xset(m, "subset"), mode = m)[["treatpost"]],
               lo = q[1], hi = q[2], n_ok = sum(!is.na(bs)))
  }))
}

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Which adjustment flips the sign?  [", SPECIES, "]\n", sep = "")
cat("==================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Burn measure: ", DOSE, " (", attr(cy, "dose_label"), ")\n", sep = "")
cat("Panel     : shared build_panel(); ONLY the absorbed fixed effects change\n")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  burn_dose_summary reports +0.0197 raw and -0.0248 adjusted for the same\n")
cat("  >= 90% burned cells. CLAUDE.md used to blame year effects. But det_rate is\n")
cat("  nearly flat by year (never-burned: 0.0789 in 2015 -> 0.0836 in 2026);\n")
cat("  the steep climb belongs to any_detection, driven by checklist volume.\n")

cat("\n-- 1. EFFECT BY FE SPECIFICATION -------------------------------------\n")
cat("  'raw' = pooled before/after minus control trend (no adjustment at all)\n")
cat("  none / cell / year / both = which fixed effects are absorbed\n\n")
print(as.data.frame(pts[, .(threshold, cells, base = round(base, 4),
                            raw = round(raw, 4), none = round(none, 4),
                            cell = round(cell, 4), year = round(year, 4),
                            both = round(both, 4))]), row.names = FALSE)
cat("\n  READ ACROSS EACH ROW. The column where the sign changes is the\n")
cat("  adjustment doing the work. 'none' should reproduce 'raw' exactly -- if it\n")
cat("  does not, the panel is not the one the raw figure came from.\n")

cat("\n-- 2. CONTINUOUS DOSE-RESPONSE ---------------------------------------\n\n")
print(as.data.frame(dose_pts[, lapply(.SD, function(x)
        if (is.numeric(x)) round(x, 4) else x)]), row.names = FALSE)

cat("\n-- 3. BOOTSTRAP CIs AT >= 90% BURNED ---------------------------------\n")
if (is.null(boots)) cat("  skipped (FE_BOOT=0)\n") else {
  cat("  cell-level bootstrap, ", B_REPS, " reps\n\n", sep = "")
  print(as.data.frame(boots[, .(mode, est = round(est, 4), ci_low = round(lo, 4),
                                ci_high = round(hi, 4), n_ok)]), row.names = FALSE)
  cat("\n  NOTE: these cluster by CELL and are PERCENTILE intervals -- wrong unit\n")
  cat("  and wrong instrument. See wild_cluster_boot.R for the real verdict.\n")
  cat("  They are kept here only to show the four modes on equal footing.\n")
}

cat("\n-- HOW TO READ THE VERDICT -------------------------------------------\n")
cat("  * flip under 'year' but not 'cell' -> a trend-adjustment story\n")
cat("  * flip under 'cell' but not 'year' -> COMPOSITION: the mix of cells\n")
cat("    contributing checklists changed between before and after, and the raw\n")
cat("    pooled number was never a within-cell comparison at all\n")
cat("  * flip needs BOTH                  -> the two interact; report as such\n")
cat("\n  For all three species the answer is CELL. Year effects do nothing.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
