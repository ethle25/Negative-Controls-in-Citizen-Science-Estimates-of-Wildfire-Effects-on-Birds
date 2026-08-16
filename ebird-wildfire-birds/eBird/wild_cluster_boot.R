# =============================================================================
# Wild cluster bootstrap for the dose-response, clustered by FIRE
#
# WHY THIS EXISTS
#   fire_cluster_se.R resamples fires from a hat (percentile bootstrap). With
#   only ~225 fire clusters that is the wrong instrument, and this project has
#   the empirical proof: OSFL's dose interval FLIPPED VERDICT on rep count alone
#   -- B=200 gave [-0.0008, +0.0523] (fails), B=2000 gave [+0.0004, +0.0510]
#   (clears). A margin of +0.0004 on a +0.0283 estimate is 1.4%, well inside
#   Monte-Carlo error of a percentile endpoint.
#
#   Worse, the percentile intervals got the WIDTH right and the POSITION wrong:
#   Wrentit's implied SE was 0.0174 against an analytic 0.0171, but the interval
#   centred on -0.0367 rather than the actual -0.0298 because the resampled
#   distribution is skewed. The symmetric interval -0.0298 +/- 1.96*0.0171 =
#   [-0.0633, +0.0037] includes zero, and always did.
#
#   Cameron-Gelbach-Miller does not resample units. Every cluster is present in
#   every replication; what varies is a random +/-1 (Rademacher) sign flip on
#   each cluster's residuals. Standard fix for few-cluster inference.
#
# WHY G_eff IS REALLY ~225, NOT ~8,000
#   The dose regressor is exactly 0 for every never-burned control, so controls
#   contribute nothing to the score for that coefficient and the cluster-robust
#   meat matrix is built from the treated fire clusters alone. The arithmetic
#   handles this on its own -- no need to drop controls, which would change the
#   estimate.
#
# REFACTORED 2026-07-30 onto panel_utils.R + infer.R. The panel construction and
# the bootstrap itself now live in one place instead of being copied into six
# scripts. Behaviour is unchanged -- verify_published.R pins the outputs.
#
# ENV
#   EBIRD_PREFIX / EBIRD_SPECIES   which species
#   WCB_B          replications for the test          (default 9999)
#   WCB_DOSE       burn measure from treatment_registry.csv (default frac_max)
#   WCB_SUBSET=1   restrict to cells >= WCB_THRESH burned, test treat x post
#   WCB_GRIDN      grid points for the inverted CI    (default 61, 0 to skip)
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
DOSE    <- Sys.getenv("WCB_DOSE", "frac_max")
B_REPS  <- as.integer(Sys.getenv("WCB_B", "9999"))
USE_SUB <- Sys.getenv("WCB_SUBSET", "0") != "0"
THRESH  <- as.numeric(Sys.getenv("WCB_THRESH", "0.50"))
GRIDN   <- as.integer(Sys.getenv("WCB_GRIDN", "61"))
out_sum <- file.path(dir_, sprintf("wild_cluster_summary_%s.txt", PREFIX))

# --- panel + clusters: two lines, was ~45 -----------------------------------
say("Building panel (", PREFIX, ", dose = ", DOSE, ")")
p <- attach_fire(build_panel(PREFIX, dose = DOSE), PREFIX)
terms <- attr(p, "dose_terms")

if (USE_SUB) {
  stopifnot(length(terms) == 1L)                 # subsets need a single dose
  p <- p[treat == 0 | get(terms[1]) >= THRESH]
  p[, treatpost := as.numeric(treat) * as.numeric(post)]
  XCOLS <- c("post_d", "treatpost"); INTEREST <- "treatpost"
} else {
  XCOLS <- c("post_d", paste0("post_", terms)); INTEREST <- paste0("post_", terms[1])
}
n_fire_cl <- uniqueN(p[treat == 1]$fire)
say("panel ", fmt(nrow(p)), " | treated fire clusters ", fmt(n_fire_cl),
    " | total clusters ", fmt(uniqueN(p$cl)))

# --- fit + inference: shared code -------------------------------------------
say("Absorbing FE")
P  <- .prep(p, XCOLS)
ki <- match(INTEREST, XCOLS)
b_hat  <- .fit(P, P$y)
V_hat  <- crve(P, as.vector(P$y - P$X %*% b_hat))
se_hat <- sqrt(V_hat[ki, ki])
say("point estimate ", sprintf("%+.4f", b_hat[ki]),
    " | fire-clustered SE ", sprintf("%.4f", se_hat))

say("WCB test of H0: effect = 0  (B=", B_REPS, ")")
p0 <- wcb_p(P, ki, b_hat[ki], se_hat, B = B_REPS, b0 = 0)
say("  p = ", sprintf("%.4f", p0))

# --- CI by test inversion ----------------------------------------------------
# No percentile shortcut exists for WCB: the CI is the set of b0 the test does
# not reject, so it costs a full bootstrap per grid point.
ci_lo <- ci_hi <- NA_real_; grid_step <- NA_real_
if (GRIDN > 0) {
  B_GRID <- min(B_REPS, 1999L)
  say("Inverting the test over ", GRIDN, " grid points (B=", B_GRID, " each)")
  grid <- seq(b_hat[ki] - 5*se_hat, b_hat[ki] + 5*se_hat, length.out = GRIDN)
  grid_step <- diff(grid)[1]
  pg <- vapply(grid, function(g) wcb_p(P, ki, b_hat[ki], se_hat, B = B_GRID, b0 = g),
               numeric(1))
  keep <- which(pg >= 0.05)
  if (length(keep)) { ci_lo <- grid[min(keep)]; ci_hi <- grid[max(keep)] }
}

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("Wild cluster bootstrap, clustered by FIRE  [", SPECIES, "]\n", sep = "")
cat("==========================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Burn measure: ", DOSE, " (", attr(p, "dose_label"), ")\n", sep = "")
cat("Target    : ", INTEREST,
    if (USE_SUB) sprintf("  (subset >= %d%% burned)", round(100*THRESH))
    else "  (continuous dose at 100% burned)", "\n", sep = "")
cat("Clusters  : ", fmt(n_fire_cl), " fires (treated) + ",
    fmt(uniqueN(p[treat == 0]$fire)), " singleton controls\n", sep = "")
cat("Reps      : ", B_REPS, " for the test",
    if (GRIDN > 0) sprintf(", %d per grid point", min(B_REPS, 1999L)) else "", "\n", sep = "")

cat("\n-- WHY THIS EXISTS ---------------------------------------------------\n")
cat("  The percentile bootstrap in fire_cluster_se.R cannot resolve this: OSFL's\n")
cat("  dose interval flipped verdict between B=200 [-0.0008,+0.0523] and B=2000\n")
cat("  [+0.0004,+0.0510] on Monte-Carlo error alone. With ~225 clusters the wild\n")
cat("  cluster bootstrap is the correct instrument -- it holds every cluster in\n")
cat("  every replication and randomises the SIGN of each cluster's residuals.\n")
cat("\n  Only treated clusters carry information about this coefficient: the\n")
cat("  regressor is exactly 0 for every never-burned control, so the effective\n")
cat("  cluster count is ", fmt(n_fire_cl), ", not the nominal total.\n", sep = "")

cat("\n-- RESULT ------------------------------------------------------------\n\n")
cat("  point estimate                : ", sprintf("%+.4f", b_hat[ki]), "\n", sep = "")
cat("  fire-clustered SE             : ", sprintf("%.4f", se_hat), "\n", sep = "")
cat("  observed t                    : ", sprintf("%+.3f", b_hat[ki]/se_hat), "\n", sep = "")
cat("  WCB p-value, H0: effect = 0   : ", sprintf("%.4f", p0), "\n", sep = "")
cat("  verdict at 5%                 : ",
    ifelse(p0 < 0.05, "REJECTS zero", "does NOT reject zero"), "\n", sep = "")
if (GRIDN > 0) {
  cat("\n  WCB 95% CI (test inversion)   : [", sprintf("%+.4f", ci_lo), ", ",
      sprintf("%+.4f", ci_hi), "]\n", sep = "")
  cat("  excludes zero                 : ",
      ifelse(!is.na(ci_lo) && (ci_lo > 0 | ci_hi < 0), "YES", "no"), "\n", sep = "")
  cat("\n  KNOWN ISSUE: the inverted CI uses fewer reps per grid point than the\n")
  cat("  direct test and is granular to ", sprintf("%.5f", grid_step), ". It has\n", sep="")
  cat("  contradicted the direct p-value before (OSFL: p=0.0799 does not reject,\n")
  cat("  while the CI excluded zero). TRUST THE DIRECT P-VALUE.\n")
}

cat("\n-- HOW TO READ -------------------------------------------------------\n")
cat("  This p-value supersedes any percentile-bootstrap verdict on the same\n")
cat("  quantity. Where they disagree, the percentile verdict was an artifact of\n")
cat("  taking raw quantiles of an unstudentised coefficient at few clusters.\n")

cat("\n-- CAVEATS -----------------------------------------------------------\n")
cat("  * Rademacher weights need G >= ~12 for enough distinct draws; ",
    fmt(n_fire_cl), " treated\n    clusters is comfortably above that.\n", sep = "")
cat("  * The null is imposed (restricted residuals) -- the recommended variant.\n")
cat("    Unrestricted WCB under-rejects with few clusters.\n")
cat("  * FE absorbed once, replications refit in the partialled space\n")
cat("    (Frisch-Waugh). Standard for fast WCB; identical coefficients.\n")
cat("  * Controls are singleton clusters, so spatial correlation among\n")
cat("    never-burned cells is still unmodelled.\n")
cat("  * A null here is NOT proof of no effect: with ", fmt(n_fire_cl),
    " clusters this test\n    has limited power. Report the point estimate too.\n", sep = "")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
