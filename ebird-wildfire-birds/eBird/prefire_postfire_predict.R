# =============================================================================
# BEFORE vs AFTER, IN THE SAME PLACES -- is fire relevant to the PREDICTION?
#
# THE PROBLEM THIS SOLVES
#   In the main predictive model the burn block is worth +0.046 AUC alone but
#   only +0.002 once coordinates are present. The reason is that burn history
#   works as a CHAPARRAL DETECTOR: fires burn chaparral, Wrentits live in
#   chaparral, so "has burned" is largely a statement about WHERE a cell is.
#
#   Restricting the sample to cells that burn AT SOME POINT removes that
#   channel. Every cell in the sample is burnable chaparral, so burn features
#   can no longer earn their keep by identifying habitat. Whatever predictive
#   value survives has to come from WHEN the fire happened relative to the
#   observation -- which is the question actually being asked.
#
# THE TAG
#   Each cell-time is labelled pre-fire or post-fire from its own cell's first
#   mapped fire. Cells must be observed BOTH before AND after, or they carry no
#   within-place contrast and would instead smuggle back a between-place one.
#
# THREE MODELS, IDENTICAL ROWS AND PARAMETERS
#   base        effort + season + vegetation, no fire information at all
#   + post      base plus a single 0/1 flag: has this cell's fire happened yet
#   + burn      base plus the full 22-column burn block
#
#   base -> +post isolates the value of knowing the fire has occurred.
#   base -> +burn is the same comparison for everything the project measures
#   about fire: extent, severity, recency, neighbourhood.
#
# A SECOND SAMPLE FOR CONTRAST
#   The same three fits are run on ALL cells, which reproduces the published
#   setting. The gap between the two samples is the point: it shows how much of
#   fire's apparent predictive value is location rather than timing.
#
# THIS IS PREDICTION, NOT CAUSATION. Delta-AUC is not an effect size, and this
# does not estimate what fire did to birds -- see Finding 5 and the negative
# controls in RESULTS_WRITEUP.md section 2.
#
# ENV  PF_PREFIX (default cfg_default_prefix())
# =============================================================================
suppressMessages({ library(data.table); library(lightgbm) })
setDTthreads(0)

.pf_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  p <- file.path(.pf_dir, "project_config.R"); if (file.exists(p)) sys.source(p, envir = globalenv())
})
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

PREFIX <- Sys.getenv("PF_PREFIX", ""); if (!nzchar(PREFIX)) PREFIX <- cfg_default_prefix()
TARGET <- "any_detection"

auc <- function(y, p) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fr <- fread(file.path(.pf_dir, sprintf("%s_cells_5km_16d.csv", PREFIX)), showProgress = FALSE)
say("frame ", format(nrow(fr), big.mark = ","), " rows")

# --- tag every cell-time pre / post its own cell's first fire ----------------
# joined on (gx, gy) -- NEVER cell_id, which is not stable across species
L <- fread(file.path(.pf_dir, "cell_fire_lookup.csv"), showProgress = FALSE)
setorder(L, gx, gy, fire_year)
first <- L[, .(fire_year = fire_year[1]), by = .(gx, gy)]
n0 <- nrow(fr)
fr <- merge(fr, first, by = c("gx","gy"), all.x = TRUE)
stopifnot(nrow(fr) == n0)
fr[, post_fire := fifelse(is.na(fire_year), 0L, as.integer(year > fire_year))]

# NA handling must match the pipeline: yrs_since_burn NA means NEVER BURNED,
# and dist_median NA means distance not reported. Dropping either would delete
# most never-burned cells -- see list_length_predict.R.
fr[is.na(yrs_since_burn), yrs_since_burn := 100]
fr[, dist_median_na := as.integer(is.na(dist_median))]
fr[is.na(dist_median), dist_median := 0]

BASE <- unlist(lapply(c("effort","season","veg"), cfg_features))
BURN <- cfg_features("burn")
BASE <- c(intersect(BASE, names(fr)), "dist_median_na")
BURN <- intersect(BURN, names(fr))
fr <- fr[complete.cases(fr[, c(BASE, BURN), with = FALSE])]
say("modelling rows ", format(nrow(fr), big.mark = ","),
    " | base ", length(BASE), " | burn ", length(BURN))

# --- the restricted sample: observed BOTH before and after its own fire ------
seen <- fr[!is.na(fire_year), .(pre = sum(post_fire == 0), post = sum(post_fire == 1)),
           by = .(gx, gy)]
both <- seen[pre > 0 & post > 0]
say(sprintf("cells with a fire: %s | observed BOTH sides: %s",
            format(nrow(seen), big.mark=","), format(nrow(both), big.mark=",")))
sub <- merge(fr, both[, .(gx, gy)], by = c("gx","gy"))
say("restricted sample ", format(nrow(sub), big.mark = ","), " cell-times")

fit <- function(d, feats, seed = 1) {
  tr <- d$year <= 2023; te <- d$year >= 2024
  if (sum(te) < 500 || length(unique(d[[TARGET]][te])) < 2) return(NA_real_)
  set.seed(seed)
  X <- as.matrix(d[, ..feats]); y <- d[[TARGET]]
  m <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                      num_leaves = 63, min_data_in_leaf = 100,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      bagging_freq = 1, verbose = -1, num_threads = 0),
                 lgb.Dataset(X[tr, ], label = y[tr]), nrounds = 400)
  auc(y[te], predict(m, X[te, ]))
}

run <- function(d, label) {
  say("fitting ", label)
  a_base <- fit(d, BASE)
  a_post <- fit(d, c(BASE, "post_fire"))
  a_burn <- fit(d, c(BASE, BURN))
  data.table(sample = label, n = nrow(d),
             cells = uniqueN(d[, paste(gx, gy)]),
             model = c("base (no fire)", "+ post_fire flag", "+ full burn block"),
             n_feat = c(length(BASE), length(BASE) + 1L, length(BASE) + length(BURN)),
             auc = c(a_base, a_post, a_burn),
             delta = c(0, a_post - a_base, a_burn - a_base))
}

res <- rbind(run(sub, "burned cells only (before vs after)"),
             run(fr,  "all cells (published setting)"))
fwrite(res, file.path(.pf_dir, sprintf("prefire_postfire_predict_%s.csv", PREFIX)))

# what the fitted model believes about the same cells pre vs post
mb <- local({
  tr <- sub$year <= 2023
  set.seed(1); X <- as.matrix(sub[, c(BASE, BURN), with = FALSE])
  m <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                      num_leaves = 63, min_data_in_leaf = 100,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      bagging_freq = 1, verbose = -1, num_threads = 0),
                 lgb.Dataset(X[tr, ], label = sub[[TARGET]][tr]), nrounds = 400)
  # compute the prediction as a COLUMN first. Calling predict() inside a `by=`
  # does not subset X to the group, so it silently returns the grand mean for
  # every group -- which is what an earlier version of this file reported.
  s <- copy(sub)
  s[, phat := predict(m, X)]
  s[, .(n = .N, obs = mean(get(TARGET)), p = mean(phat)), by = post_fire][order(post_fire)]
})

# --- does fire's predictive value DECAY with time since the fire? ------------
# One pair of models fitted once on the restricted train set, then scored
# SEPARATELY within each years-since-fire band of the test set. The band-level
# DELTA (burn - base) is the comparable quantity; raw AUC differs across bands
# because base rates and difficulty differ, so those are not comparable to each
# other. Expectation from veg_event_study.R: NDVI falls 27% at year 1 and is
# back to baseline by year 8, so if detection tracks the environment the delta
# should be largest early and fade.
by_age <- local({
  tr <- sub$year <= 2023; te <- which(sub$year >= 2024)
  y  <- sub[[TARGET]]
  mk <- function(feats) {
    set.seed(1); X <- as.matrix(sub[, ..feats])
    m <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                        num_leaves = 63, min_data_in_leaf = 100,
                        feature_fraction = 0.8, bagging_fraction = 0.8,
                        bagging_freq = 1, verbose = -1, num_threads = 0),
                   lgb.Dataset(X[tr, ], label = y[tr]), nrounds = 400)
    predict(m, X)
  }
  p_base <- mk(BASE); p_burn <- mk(c(BASE, BURN))
  d <- data.table(y = y[te], pb = p_base[te], pf = p_burn[te],
                  evt = sub$year[te] - sub$fire_year[te])
  d[, band := cut(evt, c(-100, -1, 1, 3, 5, 10, 100),
                  labels = c("before fire","0-1 yr","2-3 yr","4-5 yr","6-10 yr","11+ yr"))]
  d[!is.na(band), .(n = .N, prevalence = round(mean(y), 3),
                    auc_base = round(auc(y, pb), 4),
                    auc_burn = round(auc(y, pf), 4),
                    delta = round(auc(y, pf) - auc(y, pb), 4)), by = band][order(band)]
})

out <- file.path(.pf_dir, sprintf("prefire_postfire_predict_%s.txt", PREFIX))
sink(out)
cat("BEFORE vs AFTER A BURN -- IS FIRE RELEVANT TO THE PREDICTION?\n")
cat("=============================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, "   target ", TARGET, "\n", sep = "")
cat("Split     : temporal (train <= 2023, test 2024-2026)\n")
cat("Model     : LightGBM, identical parameters to model_burn_effect.R\n\n")
cat("  Restricting to cells that eventually burn removes the chaparral-detector\n")
cat("  channel: every cell in that sample is burnable habitat, so fire features\n")
cat("  can only earn AUC through TIMING, not through location.\n\n")
print(as.data.frame(res[, .(sample, cells, n, model, n_feat,
                            auc = round(auc, 4),
                            delta = sprintf("%+.4f", delta))]), row.names = FALSE)

cat("\n-- THE SAME CELLS, PRE vs POST (full burn model) ---------------------\n\n")
print(as.data.frame(mb[, .(post_fire, n, observed = round(obs, 4),
                           predicted = round(p, 4))]), row.names = FALSE)

cat("\n-- DOES FIRE'S PREDICTIVE VALUE DECAY WITH TIME SINCE THE FIRE? ------\n")
cat("   One model pair fitted once, scored within each band of the test set.\n")
cat("   Compare the DELTA column across bands. Raw AUCs are not comparable\n")
cat("   between bands -- base rates and difficulty differ.\n\n")
print(as.data.frame(by_age), row.names = FALSE)
cat("\n   For reference, NDVI in these cells (veg_event_study.R): -27% at year 1,\n")
cat("   recovering to pre-fire level by year 8.\n")

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  Compare the two samples. If the burn block's delta collapses in the\n")
cat("  restricted sample, its value in the published setting was LOCATION --\n")
cat("  it was telling the model which cells are chaparral, not what fire did.\n")
cat("  If the delta survives, fire timing carries real predictive information\n")
cat("  beyond habitat.\n")
cat("\n  For scale: in the published setting the burn block is worth about\n")
cat("  +0.046 AUC alone, and +0.002 once gx/gy are added.\n")
cat("\n  PREDICTION ONLY. Delta-AUC is not an effect size. The negative controls\n")
cat("  (RESULTS_WRITEUP.md section 2) show per-species fire effects on detection\n")
cat("  cannot currently be quoted, and nothing here changes that.\n")
sink()
cat(readLines(out), sep = "\n")
say("DONE")
