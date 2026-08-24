# =============================================================================
# Regression test: does the code still produce the numbers already reported?
#
# WHY THIS EXISTS
#   Refactoring analysis code is dangerous in a way that refactoring most code is
#   not: a mistake does not crash, it silently returns a different number. Every
#   figure in CLAUDE.md, in the results bundle, and in anything written up comes
#   from a specific estimator on a specific panel. If a refactor moves one of
#   them, nobody finds out until a reviewer does.
#
#   So: pin the numbers FIRST, then refactor. Run this before touching anything,
#   run it after each change. A failure means the change altered a result, not
#   that the test is wrong.
#
# WHAT IT PINS
#   Structural facts (panel sizes, cluster counts) exactly, deterministic point
#   estimates tightly, and bootstrap-based quantities loosely -- those vary run
#   to run by design and a tight tolerance would produce false alarms.
#
# TIERS
#   default        structure + all point estimates.  ~1 min
#   VERIFY_FULL=1  adds Callaway-Sant'Anna and the wild cluster bootstrap. ~5 min
#
# EXIT CODE
#   0 = all pass. 1 = at least one mismatch. Usable in a shell chain:
#     Rscript verify_published.R && echo "safe to refactor"
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
setwd(dir_)
source(file.path(dir_, "panel_utils.R"))
FULL <- Sys.getenv("VERIFY_FULL", "0") != "0"

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

.results <- list()
chk <- function(name, actual, expected, tol, source_file = "") {
  pass <- !is.na(actual) && abs(actual - expected) <= tol
  .results[[length(.results) + 1L]] <<- data.table(
    check = name, expected = expected, actual = actual,
    tol = tol, pass = pass, source = source_file)
  message(sprintf("  %-46s %+10.4f vs %+10.4f  %s",
                  name, actual, expected, if (pass) "ok" else "*** FAIL ***"))
  invisible(pass)
}

# =============================================================================
# 1. panel structure -- exact. If these move, the panel definition changed.
# =============================================================================
say("1. panel structure")
p <- build_panel("wrentit", dose = "frac_max")
chk("panel: cell-year rows",        nrow(p),                          61745, 0, "dose_burn_summary")
chk("panel: treated cells",         uniqueN(p[treat == 1]$cell_id),    1009, 0, "dose_burn_summary")
chk("panel: control cells",         uniqueN(p[treat == 0]$cell_id),    8057, 0, "dose_burn_summary")
chk("panel: pre-fire treated rate", p[treat == 1 & post == 0, sum(det)/sum(w)],
                                                                     0.1526, 5e-4, "burn_dose_summary")

pf <- attach_fire(p, "wrentit")
chk("clusters: fires behind treated cells", uniqueN(pf[treat == 1]$fire), 225, 0, "fire_cluster_summary")

lk <- fread("cell_fire_lookup.csv", showProgress = FALSE)
chk("lookup: burned cells matched",  uniqueN(lk$cell_id),               3125, 0, "cell_fire_lookup_summary")
chk("lookup: distinct fires",        uniqueN(lk$fire_event_id),          480, 0, "cell_fire_lookup_summary")
chk("lookup: max cells in one fire", max(lk[, .N, by = fire_event_id]$N), 231, 0, "cell_fire_lookup_summary")

# =============================================================================
# 2. continuous dose, all three species -- deterministic given the seed
# =============================================================================
say("2. continuous dose (the headline)")
EXP_DOSE <- c(wrentit = -0.0298, brncre = -0.0149, olsfly = +0.0283)
for (s in names(EXP_DOSE)) {
  b <- twfe(build_panel(s, dose = "frac_max"), c("post_d","post_dose"))[["post_dose"]]
  chk(sprintf("dose @100%% burned: %s", s), b, EXP_DOSE[[s]], 5e-4, "burn_dose_summary")
}

# =============================================================================
# 3. heavy-burn subsets (Wrentit) -- the undiluted estimates
# =============================================================================
say("3. heavy-burn subsets, Wrentit")
EXP_SUB <- c("0" = -0.0065, "25" = -0.0183, "50" = -0.0151, "70" = -0.0177, "90" = -0.0248)
pw <- build_panel("wrentit", dose = "frac_max")
for (thr in names(EXP_SUB)) {
  keep <- copy(pw[treat == 0 | dose >= as.numeric(thr)/100])
  keep[, treatpost := as.numeric(treat) * as.numeric(post)]
  b <- twfe(keep, c("post_d","treatpost"))[["treatpost"]]
  chk(sprintf("subset >= %s%% burned", thr), b, EXP_SUB[[thr]], 5e-4, "burn_dose_summary")
}

# =============================================================================
# 4. which fixed effect flips the sign (fe_decompose)
# =============================================================================
say("4. FE decomposition, Wrentit >= 90%")
k90 <- copy(pw[treat == 0 | dose >= 0.90])
k90[, `:=`(one = 1, treat_d = as.numeric(treat),
           treatpost = as.numeric(treat) * as.numeric(post))]
fe_fit <- function(dt, xc, mode) {
  z <- copy(dt)[, c("cell_id","year","rate","w", xc), with = FALSE]
  for (cc in c("rate", xc)) set(z, NULL, cc, as.numeric(z[[cc]]))
  if (mode != "none") for (i in seq_len(if (mode == "both") 40 else 1))
    for (cc in c("rate", xc)) {
      if (mode %in% c("cell","both")) z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
      if (mode %in% c("year","both")) z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
    }
  X <- as.matrix(z[, ..xc]); XtWX <- t(X * z$w) %*% X
  as.vector(solve(XtWX, t(X * z$w) %*% z$rate))[match("treatpost", xc)]
}
chk("FE none (= raw pooled DiD)", fe_fit(k90, c("one","treat_d","post_d","treatpost"), "none"),
    +0.0197, 5e-4, "fe_decompose_summary")
chk("FE year only",               fe_fit(k90, c("treat_d","post_d","treatpost"), "year"),
    +0.0210, 5e-4, "fe_decompose_summary")
chk("FE cell only (flips sign)",  fe_fit(k90, c("post_d","treatpost"), "cell"),
    -0.0253, 5e-4, "fe_decompose_summary")
chk("FE both",                    fe_fit(k90, c("post_d","treatpost"), "both"),
    -0.0248, 5e-4, "fe_decompose_summary")

# =============================================================================
# 5. fire-clustered standard errors -- analytic, deterministic
# =============================================================================
say("5. fire-clustered SE")
EXP_SE <- list(wrentit = c(se = 0.0171, t = -1.744),
               brncre  = c(se = 0.0104, t = -1.427),
               olsfly  = c(se = 0.0140, t = +2.018))
for (s in names(EXP_SE)) {
  q  <- attach_fire(build_panel(s, dose = "frac_max"), s)
  xc <- c("post_d","post_dose"); b <- twfe(q, xc)
  z <- copy(q)[, c("cell_id","year","rate","w", xc), with = FALSE]
  for (cc in c("rate", xc)) set(z, NULL, cc, as.numeric(z[[cc]]))
  for (i in 1:40) for (cc in c("rate", xc)) {
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year] }
  X <- as.matrix(z[, ..xc]); e <- as.vector(z$rate - X %*% b)
  Xi <- solve(t(X * z$w) %*% X); S <- rowsum(X * z$w * e, q$cl, reorder = FALSE)
  se <- sqrt(diag(Xi %*% crossprod(S) %*% Xi))[2]
  chk(sprintf("fire-clustered SE: %s", s), se, EXP_SE[[s]][["se"]], 1e-3, "wild_cluster_summary")
  chk(sprintf("fire-clustered t:  %s", s), b[[2]]/se, EXP_SE[[s]][["t"]], 0.05, "wild_cluster_summary")
}

# =============================================================================
# 6. cross-numbering safety
#
# Everything above pins ONE numbering (wrentit/brncre/olsfly share one) and one
# burn source (`frame`). That is exactly the blind spot two production bugs used:
# a re-keyed .regrid_cache/burn.rds silently mis-joined every `cache` measure,
# and .perim_attrs() joined the fire lookup on cell_id so every `perimeter`
# measure was silently wrong for any separately built species. Neither moved a
# single number the checks above look at.
#
# So: exercise both other sources, and assert the cross-species guarantee
# directly -- the same GROUND must get the same dose under either numbering.
# =============================================================================
say("6. cross-numbering safety (cache + perimeter sources)")

chk("cache dose: extent_x_meansev (wrentit)",
    twfe(build_panel("wrentit", "extent_x_meansev"),
         c("post_d","post_extent","post_extent_sev"))[["post_extent"]],
    -0.0331, 5e-4, "panel_utils::.cache_sev")

# The cache guard must FAIL LOUDLY rather than hand back some other species'
# numbering. Tested with a prefix that cannot exist, so this stays a test of the
# GUARANTEE rather than of which species happen to have a cache today -- the
# earlier version used "calthr" and started failing the moment calthr's cache was
# built, which said nothing about whether the guard still worked.
#
# This is not hypothetical: until 2026-07-31 a missing frame left the fingerprint
# NULL and .cache_sev() returned the FIRST cache on disk, so a typo'd prefix
# silently got wrentit's severity composition.
chk("cache guard errors on unknown species",
    as.numeric(inherits(tryCatch(.cache_sev("__no_such_species__"), error = function(e) e),
                        "error")),
    1, 0, "panel_utils::.cache_sev")

# ...and every BUILT species must resolve to a cache of its own numbering.
for (s in c("wrentit","brncre","olsfly","calthr")) {
  if (!file.exists(sprintf("%s_cells_5km_16d.csv", s))) next
  chk(sprintf("cache resolves for built species: %s", s),
      as.numeric(!inherits(tryCatch(.cache_sev(s), error = function(e) e), "error")),
      1, 0, "panel_utils::.cache_sev")
}

if (file.exists("calthr_cells_5km_16d.csv")) {
  gxy <- function(s) unique(fread(sprintf("%s_cells_5km_16d.csv", s),
                                  select = c("cell_id","gx","gy"), showProgress = FALSE))
  d1 <- merge(unique(build_panel("wrentit", "fire_size")[treat == 1, .(cell_id, dose)]),
              gxy("wrentit"), by = "cell_id")
  d2 <- merge(unique(build_panel("calthr",  "fire_size")[treat == 1, .(cell_id, dose)]),
              gxy("calthr"),  by = "cell_id")
  m  <- merge(d1, d2, by = c("gx","gy"))
  chk("perimeter dose: cells shared across numberings", nrow(m), 1009, 0,
      "panel_utils::.perim_attrs")
  chk("perimeter dose: max |wrentit - calthr|", max(abs(m$dose.x - m$dose.y)), 0, 1e-9,
      "panel_utils::.perim_attrs")
} else message("  calthr frame absent -- skipping the cross-numbering comparison")

# =============================================================================
# 6b. model_burn_effect.R feature blocks, derived vs the literals they replaced
#
# The predictive AUCs (LightGBM 0.8913, delta-AUC +0.0460, the permutation
# ranking) are not pinned anywhere -- re-fitting four models x two splits costs
# ~9 min, far too slow for this guard. But every one of those numbers is a
# function of WHICH COLUMNS, IN WHICH ORDER, go into the design matrix. So pin
# that instead: if the registry yields the identical vectors, the models cannot
# have seen anything different.
#
# ORDER, not just membership. lightgbm and xgboost sample features
# (feature_fraction / colsample_bytree = 0.8) off a seeded RNG, so permuting the
# columns changes which ones each tree sees and moves the AUC in the 4th decimal.
# =============================================================================
say("6b. model_burn_effect feature blocks (registry vs the replaced literals)")
{
  fr <- names(fread("wrentit_cells_5km_16d.csv", nrows = 1, showProgress = FALSE))
  LIT <- list(
    effort = c("n_checklists","effort_hours","dur_median","dist_median",
               "obs_median","frac_traveling","n_locations"),
    season = c("comp_doy","year"),
    veg    = c("ndvi_l0","ndvi_l1","ndvi_l2","ndvi_l3","evi_l0","evi_l1","evi_l2",
               "evi_l3","ndvi_prevyr","ndvi_delta_yr","ndvi_trend"),
    burn   = c("ever_burned","burn_frac_cum","burn_frac_max","burn_sev_wt","burn_decay",
               "yrs_since_burn","max_sev_prior","n_fire_epochs",
               grep("^burn_frac_y", fr, value = TRUE),
               grep("^burn_sev_y",  fr, value = TRUE),
               "burn_frac_nbr","burn_decay_nbr"))
  for (b in names(LIT))
    chk(sprintf("block '%s' identical to literal (incl. order)", b),
        as.numeric(identical(cfg_features(b), LIT[[b]])), 1, 0, "feature_registry")
  chk("block 'geo' identical to literal",
      as.numeric(identical(cfg_features("geo", active_only = FALSE), c("gx","gy"))),
      1, 0, "feature_registry")
  chk("active modelling features (+1 runtime = 43)",
      length(cfg_features(c("effort","season","veg","burn"))), 42, 0, "feature_registry")
}

# =============================================================================
# 6c. every script that CALLS a registry function must LOAD one
#
# The checks above run analysis code, so they cannot catch a script that fails at
# line 46 before doing anything. That is a real gap and it bit immediately: the
# 2026-07-31 sweep that replaced 23 hardcoded "wrentit" defaults with
# cfg_default_prefix() left `cell_fire_lookup.R` calling a function it had never
# sourced. All 50 checks still passed -- verify_published.R never EXECUTES
# cell_fire_lookup.R or dnbr_validate.R -- and it only surfaced when
# ./run_species.sh actually ran it.
#
# Static, cheap, and permanent: if a file mentions cfg_*(), it must either be a
# registry library itself or source one.
# =============================================================================
say("6c. registry-function callers can resolve them")
{
  libs <- c("project_config.R", "panel_utils.R", "infer.R")
  bad <- character()
  for (f in setdiff(list.files(".", "\\.R$"), c(libs, "verify_published.R"))) {
    src <- readLines(f, warn = FALSE)
    code <- grep("^\\s*#", src, invert = TRUE, value = TRUE)
    if (!any(grepl("\\bcfg_[a-z_]+\\(", code))) next
    if (!any(grepl("source\\(.*(project_config|panel_utils)\\.R", code))) bad <- c(bad, f)
  }
  if (length(bad)) message("  !! call cfg_*() without sourcing a registry file: ",
                           paste(bad, collapse = ", "))
  chk("scripts calling cfg_*() that never source it", length(bad), 0, 0, "static scan")
}

# =============================================================================
# 7. FULL tier: Callaway-Sant'Anna and the wild cluster bootstrap
#    Both are stochastic; tolerances are deliberately loose.
# =============================================================================
if (FULL) {
  say("7. Callaway-Sant'Anna (VERIFY_FULL=1)")
  if (requireNamespace("did", quietly = TRUE)) {
    EXP_CS <- c(wrentit = -0.0208, brncre = -0.0022, olsfly = +0.0174)
    for (s in names(EXP_CS)) {
      out <- system2("Rscript", "cs_estimator.R", stdout = TRUE, stderr = FALSE,
                     env = c(sprintf("EBIRD_PREFIX=%s", s), "CS_THRESH=0.50", "CS_BITERS=200"))
      f <- sprintf("cs_summary_%s_50.txt", s)
      v <- as.numeric(sub(".*simple ATT\\) *: *([-+0-9.]+).*", "\\1",
                          grep("simple ATT", readLines(f), value = TRUE)[1]))
      chk(sprintf("CS simple ATT >=50%%: %s", s), v, EXP_CS[[s]], 2e-3, "cs_summary")
    }
  } else message("  did package not installed, skipping CS checks")

  say("8. wild cluster bootstrap p-values (stochastic, loose tolerance)")
  EXP_P <- c(wrentit = 0.110, brncre = 0.251, olsfly = 0.080)
  for (s in names(EXP_P)) {
    f <- sprintf("wild_cluster_summary_%s.txt", s)
    if (!file.exists(f)) { message("  missing ", f, ", skipping"); next }
    v <- as.numeric(sub(".*: *([0-9.]+) *$", "\\1",
                        grep("WCB p-value", readLines(f), value = TRUE)[1]))
    chk(sprintf("WCB p (from file): %s", s), v, EXP_P[[s]], 0.03, "wild_cluster_summary")
  }
  # --- 9. the four scripts retrofitted 2026-07-31 ---------------------------
  # These build their own panel shapes (event time, 16-day phenology, severity
  # composition), so they cannot be recomputed from build_panel() alone without
  # duplicating their logic here. Instead RE-RUN each and read its report back.
  # Slower, but it verifies the whole script rather than a fragment of it.
  say("9. re-running the four retrofitted scripts (~5 min)")
  # `section` anchors the search to ONE SPECIES' block.
  #
  # These reports list every active species, and their ORDER is registry-driven --
  # so "the first matching line" silently means "whichever species happens to sort
  # first". Activating Golden-crowned Kinglet on 2026-07-31 did exactly that: it
  # has migratory = 1 and sorts before Olive-sided Flycatcher, so the migration
  # check read gockin's -17.13 against OSFL's pinned +15.63 and reported a moved
  # published number. OSFL's value was untouched. A check that pins a per-species
  # figure must NAME THE SPECIES.
  runchk <- function(script, env, file, pattern, field, expected, tol, label,
                     section = NULL) {
    system2("Rscript", script, stdout = FALSE, stderr = FALSE, env = env)
    if (!file.exists(file)) { chk(label, NA_real_, expected, tol, script); return(invisible()) }
    txt <- readLines(file)
    # `section` may be several anchors applied IN SEQUENCE, so a check can say
    # "inside MODEL B, then inside the Wrentit block" rather than trusting that
    # Wrentit happens to sort first.
    for (s in section) {
      i <- grep(s, txt)[1]
      if (is.na(i)) { chk(label, NA_real_, expected, tol, basename(file)); return(invisible()) }
      txt <- txt[i:length(txt)]
    }
    ln <- grep(pattern, txt, value = TRUE)[1]
    v  <- if (is.na(ln)) NA_real_ else
            as.numeric(regmatches(ln, gregexpr("[-+][0-9]*\\.[0-9]+", ln))[[1]][field])
    chk(label, v, expected, tol, basename(file))
  }
  runchk("event_study.R", character(), "event_study_summary.txt",
         "^   1 ", 1, -0.0312, 5e-4, "event_study: Wrentit k=1",
         section = "^=== Wrentit")
  # filename gained a species suffix 2026-07-31 (it used to be a single shared
  # file that a second species would overwrite); the pinned NUMBER is unchanged
  runchk("did_burn_effect.R", character(), "burn_did_summary_wrentit.txt",
         "beta\\(post\\)", 1, -0.0024, 5e-4, "did_burn_effect: TWFE beta")
  runchk("migration_phenology.R", character(), "migration_phenology_summary.txt",
         "^  last_doy ", 1, 15.63, 0.05, "migration: OSFL last_doy shift",
         section = "=== Olive-sided Flycatcher")
  # NB pin a real number. An earlier version passed expected = NA with tol = Inf
  # as a "did it run" check -- but abs(NA - NA) <= Inf is NA, which chk() treats
  # as a failure, so the check could never pass.
  # The pattern must require a SIGNED NUMBER. "post x extent  " alone also
  # matched the explanatory header line above the table, which has no numbers, so
  # the check extracted NA and failed for the wrong reason.
  runchk("severity_effect.R", character(), "severity_effect_summary.txt",
         "post x extent +[-+][0-9]", 1, -0.0458, 5e-4,
         "severity_effect: Wrentit extent (Model B)",
         section = c("MODEL B", "=== Wrentit"))

} else {
  say("7-9. skipped (set VERIFY_FULL=1 for CS, wild cluster, and the four")
  say("     scripts that build their own panel shapes)")
}

# =============================================================================
# report
# =============================================================================
r <- rbindlist(.results)
cat("\n===============================================================\n")
cat(sprintf("VERIFY: %d checks, %d passed, %d FAILED\n",
            nrow(r), sum(r$pass), sum(!r$pass)))
cat("===============================================================\n")
if (any(!r$pass)) {
  cat("\nFAILURES -- a refactor changed a published number:\n\n")
  print(as.data.frame(r[pass == FALSE, .(check, expected, actual,
                                         diff = round(actual - expected, 5), source)]),
        row.names = FALSE)
  cat("\nDo NOT commit this state. Either the change is wrong, or the number in\n")
  cat("CLAUDE.md / the results bundle needs updating -- decide which, deliberately.\n")
}
fwrite(r, "verify_published_last_run.csv")
cat("\nfull results -> verify_published_last_run.csv\n")
quit(save = "no", status = if (any(!r$pass)) 1L else 0L)
