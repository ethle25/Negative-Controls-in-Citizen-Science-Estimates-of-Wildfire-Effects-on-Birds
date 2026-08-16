# =============================================================================
# Does burn severity move the prediction? -- four models on the cell x time frame
#
# QUESTION (from the prediction setting): given burn history and the vegetation
# run-up for a 5 km cell at a 16-day step, will Wrentit be observed? And then:
# HOW MUCH does burn severity actually contribute to that prediction?
#
# MODELS
#   lightgbm  gradient-boosted trees, leaf-wise
#   xgboost   gradient-boosted trees, level-wise
#   gbm       Friedman's original GBDT (the reference implementation)
#   glm       logistic regression, standardised features (the linear baseline)
#   NOTE lightgbm and xgboost ARE GBDTs; `gbm` is included as a third, distinct
#   implementation rather than a fourth model family.
#
# HOW THE BURN CONTRIBUTION IS MEASURED -- four independent ways, because
# feature importance alone is not an effect estimate:
#   1. DELTA-AUC     fit each model with and without the burn block. The gain is
#                    the honest "what does burn history buy you" number.
#   2. PERMUTATION   shuffle the burn block in the test set, measure AUC loss.
#   3. COUNTERFACTUAL set every burn feature to its unburned state and measure
#                    the shift in predicted probability.
#   4. PARTIAL DEPENDENCE over burn severity and years-since-burn.
#
# AND THE CONFOUND. Burned cells show a much higher detection rate than
# never-burned cells, but MTBS fires concentrate in chaparral, which is Wrentit
# habitat -- so the raw contrast is mostly habitat, not fire. Two guards:
#   * a WITH-GEOGRAPHY variant, where the model can use cell coordinates. Burn
#     features then have to beat "where is this cell", not stand in for it.
#   * a WITHIN-CELL before/after difference-in-differences, which removes cell
#     identity entirely and is the only estimate here that is not confounded by
#     which cells burn.
#
# VALIDATION SPLITS -- never random; cells repeat across time steps.
#   temporal : train year <= 2023, test 2024-2026  (forecasting)
#   spatial  : cells split 75/25, no cell in both  (generalising to new places)
# =============================================================================
suppressMessages({ library(data.table); library(lightgbm); library(xgboost); library(gbm) })
setDTthreads(0)
set.seed(1)

t0 <- Sys.time()
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
# NB: do NOT name this `log` -- it would mask base::log and silently break any
# log() arithmetic (logloss did exactly that: log(p) returned the logger's NULL).
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
# panel_utils.R pulls in project_config.R too, so this one line provides both the
# feature registry (cfg_features) and the treatment registry (cfg_treatments,
# .parse_terms, .terms_vars). This script previously loaded neither -- it was the
# only analysis script with no registry access at all.
source(file.path(dir_, "panel_utils.R"))
PREFIX  <- cfg_default_prefix()
SPECIES <- cfg_default_species()

# EFFORT NORMALISATION. any_detection is NOT effort-neutral: it climbs from
# 0.096 at one checklist to 0.489 at twenty, purely because more people looked.
# Two ways to handle that, both reported:
#   "covariate" - keep every cell-time, give the model the effort columns and
#                 let it condition on them. Best raw accuracy, but the target
#                 still partly measures observer volume.
#   "fixed1"    - keep only cell-times with EXACTLY ONE checklist. Then
#                 any_detection means "did one birding visit find the species",
#                 which is effort-constant by construction. Smaller and harder,
#                 but it is the honest binary question.
EFFORT_MODE <- Sys.getenv("EFFORT_MODE", "covariate")
stopifnot(EFFORT_MODE %in% c("covariate", "fixed1"))
tag <- sprintf("%s_%s", PREFIX, EFFORT_MODE)
# A non-default burn block must land in its own file. BURN_MEASURE did not enter
# the tag, so running the same species twice with different burn blocks silently
# overwrote the first result -- and the two are not comparable, which is the
# whole point of running both. The default path is unchanged, so every existing
# burn_effect_results_<prefix>_<mode>.csv keeps its name.
if (nzchar(Sys.getenv("BURN_MEASURE")) && Sys.getenv("BURN_MEASURE") != "block")
  tag <- sprintf("%s_%s", tag, Sys.getenv("BURN_MEASURE"))

frame_csv <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX))
out_sum   <- file.path(dir_, sprintf("burn_effect_summary_%s.txt", tag))
out_res   <- file.path(dir_, sprintf("burn_effect_results_%s.csv", tag))
stopifnot(file.exists(frame_csv))

# --- metrics ------------------------------------------------------------------
auc <- function(y, p) {                      # Mann-Whitney U, no extra deps
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)   # n1,n0 numeric: 1.4e5*4.8e5 overflows integer
}
logloss <- function(y, p) { p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
                            -mean(y * log(p) + (1 - y) * log(1 - p)) }
brier <- function(y, p) mean((y - p)^2)

say("Reading frame")
d <- fread(frame_csv, showProgress = FALSE)
say("  ", fmt(nrow(d)), " cell-time units")

# =============================================================================
# 1. feature blocks
# =============================================================================
# Blocks come from feature_registry.csv, not from five literal vectors. The
# registry rows are ordered to reproduce the previous hardcoded vectors EXACTLY,
# order included -- lightgbm and xgboost sample features off a seeded RNG
# (feature_fraction / colsample_bytree = 0.8), so a permuted column set is a
# different model. verify_published.R section 6b pins that equality, which is the
# only cheap way to guard AUCs that cost 9 minutes to recompute.
EFFORT <- cfg_features("effort")
SEASON <- cfg_features("season")
VEG    <- cfg_features("veg")
BURN   <- cfg_features("burn")
GEO    <- cfg_features("geo", active_only = FALSE)
TARGET <- "any_detection"

# BURN_MEASURE swaps the burn block for a named row of treatment_registry.csv,
# so the predictive side can ask the same "which burn representation?" question
# the causal side has always been able to ask. Default "block" = the full 22
# registry burn columns, i.e. the published behaviour.
#
# A registry measure contributes the columns its `terms` expressions REFERENCE
# (resolved with all.vars(), same as panel_utils.R), not a constructed dose:
# the models take raw columns and build their own interactions, so handing them
# a single collapsed dose would throw information away rather than swap it.
BURN_MEASURE <- Sys.getenv("BURN_MEASURE", "block")
if (BURN_MEASURE != "block") {
  tr <- cfg_treatments()[name == BURN_MEASURE]
  if (!nrow(tr)) stop("BURN_MEASURE '", BURN_MEASURE, "' is not in treatment_registry.csv. ",
                      "Available: block, ", paste(cfg_treatments()$name, collapse = ", "))
  if (tr$source[1] != "frame")
    stop("BURN_MEASURE '", BURN_MEASURE, "' is sourced from '", tr$source[1],
         "'. The predictive frame carries only `frame`-sourced burn columns; ",
         "cache/perimeter measures are not in *_cells_5km_16d.csv.")
  BURN <- .terms_vars(.parse_terms(tr$terms[1]))
  miss <- setdiff(BURN, names(d))
  if (length(miss)) stop("BURN_MEASURE '", BURN_MEASURE, "' references columns not in the frame: ",
                         paste(miss, collapse = ", "))
  say("burn block <- '", BURN_MEASURE, "' (", tr$label[1], "): ",
      length(BURN), " column(s) -- ", paste(BURN, collapse = ", "))
}
if (!length(BURN)) stop("burn block is empty; check feature_registry.csv / BURN_MEASURE")

# every registry column must actually exist in the frame, or the design matrix is
# silently short and the delta-AUC no longer means what the report says
miss <- setdiff(c(EFFORT, SEASON, VEG, BURN, GEO), names(d))
if (length(miss)) stop("feature_registry.csv names columns absent from ", basename(frame_csv),
                       ": ", paste(miss, collapse = ", "), "\n  Run cfg_check() to compare them.")

# NA handling, applied once so all four models see an identical matrix:
#   yrs_since_burn NA == "never burned" -> a sentinel far outside the observed
#     range, with ever_burned carrying the indicator (glm and gbm would
#     otherwise drop 66% of rows via complete cases)
#   dist_median NA == every checklist in the unit was Stationary -> 0 + flag
#   ndvi NA (0.4%) -> drop, so the comparison is complete-case for everyone
d[is.na(yrs_since_burn), yrs_since_burn := 100]
d[, dist_median_na := as.integer(is.na(dist_median))]
d[is.na(dist_median), dist_median := 0]
EFFORT <- c(EFFORT, "dist_median_na")
n_before <- nrow(d)
d <- d[complete.cases(d[, c(VEG, EFFORT, BURN, SEASON), with = FALSE])]
say("  dropped ", fmt(n_before - nrow(d)), " rows with missing vegetation (",
    sprintf("%.2f%%", 100 * (n_before - nrow(d)) / n_before), ")")

if (EFFORT_MODE == "fixed1") {
  n0 <- nrow(d)
  d <- d[n_checklists == 1L]
  say("  fixed-effort mode: ", fmt(n0), " -> ", fmt(nrow(d)),
      " cell-times with exactly one checklist")
}

FEATS_FULL <- c(EFFORT, SEASON, VEG, BURN)
FEATS_NOBURN <- c(EFFORT, SEASON, VEG)
# drop zero-variance columns (n_checklists and n_locations become constants once
# the frame is restricted to single-checklist units, and a constant column makes
# glm rank-deficient)
const <- names(which(sapply(d[, ..FEATS_FULL], function(v) {
  v <- v[!is.na(v)]; length(v) == 0 || max(v) == min(v) })))
if (length(const)) {
  say("  dropping constant features: ", paste(const, collapse = ", "))
  FEATS_FULL   <- setdiff(FEATS_FULL, const)
  FEATS_NOBURN <- setdiff(FEATS_NOBURN, const)
  EFFORT <- setdiff(EFFORT, const); BURN <- setdiff(BURN, const)
  VEG <- setdiff(VEG, const); SEASON <- setdiff(SEASON, const)
}
say("features: ", length(FEATS_FULL), " full | ", length(FEATS_NOBURN), " without burn | burn block = ", length(BURN))

# =============================================================================
# 2. splits
# =============================================================================
splits <- list()
splits$temporal <- list(train = d$year <= 2023, test = d$year >= 2024)
cellids <- unique(d$cell_id)
te_cells <- sample(cellids, floor(0.25 * length(cellids)))
splits$spatial <- list(train = !(d$cell_id %in% te_cells), test = d$cell_id %in% te_cells)
for (s in names(splits))
  say(sprintf("split %-9s train %s / test %s  (test prevalence %.3f)", s,
      fmt(sum(splits[[s]]$train)), fmt(sum(splits[[s]]$test)),
      mean(d[[TARGET]][splits[[s]]$test])))

# =============================================================================
# 3. fitting
# =============================================================================
fit_predict <- function(algo, feats, tr, te) {
  X <- as.matrix(d[, ..feats]); y <- d[[TARGET]]
  Xtr <- X[tr, , drop = FALSE]; ytr <- y[tr]
  Xte <- X[te, , drop = FALSE]; yte <- y[te]
  p <- switch(algo,
    lightgbm = {
      m <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                          num_leaves = 63, min_data_in_leaf = 100,
                          feature_fraction = 0.8, bagging_fraction = 0.8,
                          bagging_freq = 1, verbosity = -1, num_threads = 0),
                     lgb.Dataset(Xtr, label = ytr), nrounds = 400)
      predict(m, Xte)
    },
    xgboost = {
      m <- xgb.train(list(objective = "binary:logistic", eval_metric = "auc",
                          eta = 0.05, max_depth = 7, min_child_weight = 20,
                          subsample = 0.8, colsample_bytree = 0.8, nthread = 0),
                     xgb.DMatrix(Xtr, label = ytr), nrounds = 400, verbose = 0)
      predict(m, xgb.DMatrix(Xte))
    },
    gbm = {
      # Friedman's GBDT is single-threaded and far slower; a 150k-row training
      # subsample keeps it tractable. Reported as such -- it is a reference
      # implementation here, not a tuned competitor.
      idx <- which(tr); if (length(idx) > 150000) idx <- sample(idx, 150000)
      df <- as.data.frame(d[idx, ..feats]); df$.y <- y[idx]
      m <- gbm(.y ~ ., data = df, distribution = "bernoulli", n.trees = 300,
               interaction.depth = 5, shrinkage = 0.05, n.minobsinnode = 50,
               bag.fraction = 0.8, verbose = FALSE, n.cores = 1)
      predict(m, as.data.frame(d[te, ..feats]), n.trees = 300, type = "response")
    },
    glm = {
      mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv == 0] <- 1
      Ztr <- scale(Xtr, mu, sdv); Zte <- scale(Xte, mu, sdv)
      df <- as.data.frame(Ztr); df$.y <- ytr
      m <- glm(.y ~ ., data = df, family = binomial())
      as.numeric(predict(m, as.data.frame(Zte), type = "response"))
    })
  list(p = p, y = yte)
}

ALGOS <- c("lightgbm", "xgboost", "gbm", "glm")
res <- list(); preds <- list()
for (sp in names(splits)) {
  tr <- splits[[sp]]$train; te <- splits[[sp]]$test
  for (a in ALGOS) {
    for (variant in c("full", "noburn")) {
      feats <- if (variant == "full") FEATS_FULL else FEATS_NOBURN
      t1 <- Sys.time()
      o <- fit_predict(a, feats, tr, te)
      el <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
      res[[length(res) + 1L]] <- data.table(
        split = sp, algo = a, variant = variant, n_feat = length(feats),
        auc = auc(o$y, o$p), logloss = logloss(o$y, o$p), brier = brier(o$y, o$p),
        secs = round(el, 1))
      if (variant == "full") preds[[paste(sp, a)]] <- o
      say(sprintf("  %-9s %-9s %-7s AUC %.4f  (%.0fs)", sp, a, variant,
                  tail(res, 1)[[1]]$auc, el))
    }
  }
  # with-geography probe: can burn features survive the model knowing WHERE it is?
  for (a in c("lightgbm")) {
    o <- fit_predict(a, c(FEATS_FULL, GEO), tr, te)
    res[[length(res) + 1L]] <- data.table(split = sp, algo = paste0(a, "+geo"),
      variant = "full", n_feat = length(FEATS_FULL) + 2,
      auc = auc(o$y, o$p), logloss = logloss(o$y, o$p), brier = brier(o$y, o$p), secs = NA)
    o2 <- fit_predict(a, c(FEATS_NOBURN, GEO), tr, te)
    res[[length(res) + 1L]] <- data.table(split = sp, algo = paste0(a, "+geo"),
      variant = "noburn", n_feat = length(FEATS_NOBURN) + 2,
      auc = auc(o2$y, o2$p), logloss = logloss(o2$y, o2$p), brier = brier(o2$y, o2$p), secs = NA)
    say(sprintf("  %-9s %-9s geo probe: full %.4f / noburn %.4f", sp, a,
                tail(res, 2)[[1]]$auc, tail(res, 1)[[1]]$auc))
  }
}
R <- rbindlist(res)
fwrite(R, out_res)

# =============================================================================
# 4. permutation importance + counterfactual, on the lightgbm temporal model
# =============================================================================
say("Permutation importance (burn block, temporal split)")
tr <- splits$temporal$train; te <- splits$temporal$test
X <- as.matrix(d[, ..FEATS_FULL]); y <- d[[TARGET]]
mfull <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                        num_leaves = 63, min_data_in_leaf = 100,
                        feature_fraction = 0.8, bagging_fraction = 0.8,
                        bagging_freq = 1, verbosity = -1, num_threads = 0),
                   lgb.Dataset(X[tr, ], label = y[tr]), nrounds = 400)
Xte <- X[te, , drop = FALSE]; yte <- y[te]
base_auc <- auc(yte, predict(mfull, Xte))

perm_block <- function(cols) {
  Z <- Xte; ii <- sample(nrow(Z))
  for (cc in cols) Z[, cc] <- Z[ii, cc]
  base_auc - auc(yte, predict(mfull, Z))
}
imp <- data.table(block = c("effort","season","vegetation","burn"),
                  auc_drop = c(perm_block(EFFORT), perm_block(SEASON),
                               perm_block(VEG), perm_block(BURN)))
imp_single <- data.table(feature = BURN,
                         auc_drop = sapply(BURN, function(f) perm_block(f)))
setorder(imp_single, -auc_drop)

say("Counterfactual: set the burn block to its unburned state")
Zcf <- Xte
for (cc in BURN) Zcf[, cc] <- 0
Zcf[, "yrs_since_burn"] <- 100          # sentinel for "never burned"
p_obs <- predict(mfull, Xte); p_cf <- predict(mfull, Zcf)
burned_te <- d$ever_burned[te] == 1
cf <- data.table(
  group = c("all test units", "units whose cell had burned"),
  n = c(length(p_obs), sum(burned_te)),
  mean_p_observed = c(mean(p_obs), mean(p_obs[burned_te])),
  mean_p_unburned = c(mean(p_cf),  mean(p_cf[burned_te])),
  mean_shift = c(mean(p_obs - p_cf), mean((p_obs - p_cf)[burned_te])))

say("Partial dependence over burn severity and years-since-burn")
pdp <- function(col, grid) {
  s <- sample(nrow(Xte), min(20000, nrow(Xte)))
  Z <- Xte[s, , drop = FALSE]
  data.table(feature = col, value = grid,
             mean_pred = sapply(grid, function(g) { Z[, col] <- g; mean(predict(mfull, Z)) }))
}
pd <- rbind(
  pdp("burn_sev_wt",    c(0, 0.1, 0.25, 0.5, 1, 1.5, 2, 3)),
  pdp("yrs_since_burn", c(0.5, 1, 2, 3, 5, 8, 12, 20, 30, 100)),
  pdp("burn_frac_cum",  c(0, 0.1, 0.25, 0.5, 0.75, 1, 1.5, 2)))

# =============================================================================
# 5. within-cell before/after -- the estimate that is NOT confounded by which
#    cells burn, because every comparison is inside a single cell
# =============================================================================
say("Within-cell before/after (difference-in-differences)")
dd <- d[, .(cell_id, year, comp_doy, ever_burned, any_detection, n_detections,
            n_checklists)]
chg <- dd[, .(has_pre = any(ever_burned == 0), has_post = any(ever_burned == 1)),
          by = cell_id][has_pre & has_post]
did_cells <- chg$cell_id
treat <- dd[cell_id %in% did_cells]
tr_stat <- treat[, .(units = .N, det = sum(n_detections), chk = sum(n_checklists),
                     any_det = mean(any_detection)),
                 by = .(period = ifelse(ever_burned == 1, "after burn", "before burn"))]
tr_stat[, det_rate := det / chk]

# controls: cells that never burned at any point in the study, over the same years
ctrl_cells <- dd[, .(ever = any(ever_burned == 1)), by = cell_id][ever == FALSE]$cell_id
yr_split <- treat[ever_burned == 1, min(year)]
ctrl <- dd[cell_id %in% ctrl_cells]
ct_stat <- ctrl[, .(units = .N, det = sum(n_detections), chk = sum(n_checklists),
                    any_det = mean(any_detection)),
                by = .(period = ifelse(year >= yr_split, "after burn", "before burn"))]
ct_stat[, det_rate := det / chk]

# =============================================================================
# 6. report
# =============================================================================
sink(out_sum)
cat("Does burn severity move the prediction?   [", SPECIES, " | effort: ", EFFORT_MODE, "]\n", sep="")
cat("=======================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Frame     : ", basename(frame_csv), "  (", fmt(nrow(d)), " cell-time units)\n", sep = "")
cat("Effort    : ", EFFORT_MODE,
    if (EFFORT_MODE == "fixed1") "  (exactly one checklist per unit -- effort-constant)" else "  (all units; effort as covariates)", "\n", sep = "")
cat("Target    : ", TARGET, "  (prevalence ", sprintf("%.3f", mean(d[[TARGET]])), ")\n", sep = "")
cat("Features  : ", length(FEATS_FULL), " (effort ", length(EFFORT), ", season ",
    length(SEASON), ", vegetation ", length(VEG), ", burn ", length(BURN), ")\n", sep = "")
cat("Runtime   : ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min\n", sep = "")

cat("\n-- 1. MODEL PERFORMANCE ----------------------------------------------\n")
print(as.data.frame(R[order(split, algo, variant)]), row.names = FALSE)

cat("\n-- 2. WHAT THE BURN BLOCK BUYS (delta-AUC, full minus no-burn) -------\n")
w <- dcast(R, split + algo ~ variant, value.var = "auc")
w[, delta_auc := round(full - noburn, 4)]
print(as.data.frame(w[order(split, algo)]), row.names = FALSE)
cat("\n  This is the honest headline number: how much predictive power the whole\n")
cat("  burn history adds once effort, season and vegetation are already in the\n")
cat("  model. Vegetation is downstream of fire, so it already carries part of\n")
cat("  the fire signal -- this delta is the burn-specific remainder.\n")

cat("\n-- 3. PERMUTATION IMPORTANCE BY BLOCK (lightgbm, temporal) -----------\n")
cat("  baseline test AUC: ", sprintf("%.4f", base_auc), "\n", sep = "")
print(as.data.frame(imp[order(-auc_drop)]), row.names = FALSE)
cat("\n  top burn features individually:\n")
print(as.data.frame(head(imp_single, 10)), row.names = FALSE)

cat("\n-- 4. COUNTERFACTUAL: what if nothing had burned? ---------------------\n")
print(as.data.frame(cf), row.names = FALSE)
cat("\n  Every burn feature is set to its unburned state and the same fitted\n")
cat("  model is re-scored. The shift is the model's own belief about how much\n")
cat("  fire history changes the chance of observing a Wrentit.\n")

cat("\n-- 5. PARTIAL DEPENDENCE ---------------------------------------------\n")
print(as.data.frame(pd), row.names = FALSE)

cat("\n-- 6. THE CONFOUND, AND THE ONLY CLEAN ESTIMATE ----------------------\n")
cat("  Raw contrast across cells (from the frame summary): burned cells show a\n")
cat("  much higher detection rate than never-burned cells. That is mostly\n")
cat("  HABITAT -- MTBS fires concentrate in chaparral, which is where Wrentits\n")
cat("  live -- not a fire effect.\n")
cat("\n  WITHIN-CELL before/after, cells that burned during the study window\n")
cat("  (n = ", fmt(length(did_cells)), " cells; every comparison is inside one cell, so\n", sep = "")
cat("  cell identity and habitat cancel):\n")
print(as.data.frame(tr_stat[, .(period, units, chk, det, det_rate = round(det_rate, 4),
                                any_det = round(any_det, 4))]), row.names = FALSE)
cat("\n  NEVER-BURNED control cells over the same period split (year >= ", yr_split,
    "):\n", sep = "")
print(as.data.frame(ct_stat[, .(period, units, chk, det, det_rate = round(det_rate, 4),
                                any_det = round(any_det, 4))]), row.names = FALSE)
t_af <- tr_stat[period == "after burn"]$det_rate; t_bf <- tr_stat[period == "before burn"]$det_rate
c_af <- ct_stat[period == "after burn"]$det_rate; c_bf <- ct_stat[period == "before burn"]$det_rate
cat("\n  treated change : ", sprintf("%+.4f", t_af - t_bf),
    "   control change : ", sprintf("%+.4f", c_af - c_bf), "\n", sep = "")
cat("  DIFFERENCE-IN-DIFFERENCES : ", sprintf("%+.4f", (t_af - t_bf) - (c_af - c_bf)),
    " detection rate\n", sep = "")
cat("\n  CAVEATS on this DiD, which is descriptive and not a causal estimate:\n")
cat("   * no standard errors here; the effective n is CELLS (", fmt(length(did_cells)),
    "), not units\n", sep = "")
cat("   * burned cells may differ in effort composition before vs after a fire\n")
cat("   * time-since-burn is pooled -- a 1-year-old burn and a 10-year-old burn\n")
cat("     are averaged together, and the partial dependence in section 5 shows\n")
cat("     those are not the same thing\n")

cat("\n-- 7. READING THIS ---------------------------------------------------\n")
cat("  Use delta-AUC (section 2) for 'how much does burn history improve\n")
cat("  prediction'. Use the within-cell DiD (section 6) for 'what does fire do\n")
cat("  to Wrentit detection'. They answer different questions and should not be\n")
cat("  quoted interchangeably. Permutation importance (section 3) ranks features\n")
cat("  inside one fitted model and is NOT an effect size.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE in ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
