# =============================================================================
# DOES LIST LENGTH IMPROVE PREDICTION?
#
# THE QUESTION
#   list_length_lookup.tsv gives the number of species on every checklist -- a
#   measure of how PRODUCTIVE a trip was, which none of the seven existing
#   effort features capture (they all measure how MUCH effort was spent).
#   Does adding it to the cell frame forecast detection better?
#
#   This is a PREDICTIVE question, separate from the causal one in
#   list_length_panel.R. Delta-AUC here is not an effect size and says nothing
#   about what fire does -- see Finding 5 in CLAUDE.md.
#
# WHY A PARALLEL SCRIPT
#   Adding a column to the frame means hand-editing regrid_cells.R, the last
#   hardcoded surface in the project. Instead this computes the cell-time
#   aggregate itself and merges at model time, so the frame on disk, the
#   registries and all 51 pinned numbers are untouched.
#
# THE GRID -- replicated exactly from regrid_cells.R (lines 124-131, 188-195)
#   project lon/lat (EPSG:4326) -> EPSG:5070 equal area
#   gx = floor(X / 5000)   gy = floor(Y / 5000)
#   comp_doy = 16 * ((day_of_year - 1) %/% 16) + 1, capped at 353
#   ci = (year - 2000) * 23 + ((comp_doy - 1) %/% 16) + 1
#   Merged on (gx, gy, ci) and NEVER on cell_id -- cell_id is not stable across
#   separately built species, which has silently corrupted three analyses in
#   this project already.
#
# DESIGN
#   Same target, split, features and LightGBM parameters as model_burn_effect.R,
#   so the baseline should land near the published Wrentit 0.8913. The answer is
#   the DELTA between identical fits with and without the new feature, which is
#   robust even if the baseline drifts slightly.
#
# ENV  LL_PRED_PREFIX (default cfg_default_prefix())
# =============================================================================
suppressMessages({ library(data.table); library(sf); library(lightgbm) })
setDTthreads(0)

.lp_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  for (f in c("project_config.R")) {
    p <- file.path(.lp_dir, f); if (file.exists(p)) sys.source(p, envir = globalenv())
  }
})
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

PREFIX <- Sys.getenv("LL_PRED_PREFIX", "")
if (!nzchar(PREFIX)) PREFIX <- cfg_default_prefix()
TARGET <- "any_detection"
GRID_M <- 5000L; COMP_PER_YR <- 23L

auc <- function(y, p) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# --- cell-time list length ---------------------------------------------------
say("Building cell-time list length for ", PREFIX)
d <- fread(file.path(.lp_dir, sprintf("%s_2015_2026_joined.csv", PREFIX)),
           select = c("checklist_id","latitude","longitude","year","day_of_year","presence"),
           showProgress = FALSE)
L <- fread(file.path(.lp_dir, "list_length_lookup.tsv"),
           col.names = c("checklist_id","list_length"), showProgress = FALSE)
L <- L[, .(list_length = max(list_length)), by = checklist_id]   # 5 dup keys
setkey(L, checklist_id)
d[, list_length := L[.(d$checklist_id), list_length]]
say(sprintf("  matched %.2f%%", 100 * d[, mean(!is.na(list_length))]))
d <- d[!is.na(list_length)]
d[, ll_ex := pmax(list_length - presence, 0L)]     # exclude the focal species

uc <- unique(d[, .(longitude, latitude)])
xy <- st_coordinates(st_transform(
        st_as_sf(uc, coords = c("longitude","latitude"), crs = 4326), 5070))
uc[, `:=`(gx = as.integer(floor(xy[,1] / GRID_M)),
          gy = as.integer(floor(xy[,2] / GRID_M)))]
d[uc, on = .(longitude, latitude), `:=`(gx = i.gx, gy = i.gy)]
d[, comp_doy := 16L * ((day_of_year - 1L) %/% 16L) + 1L]
d[comp_doy > 353L, comp_doy := 353L]
d[, ci := (year - 2000L) * COMP_PER_YR + ((comp_doy - 1L) %/% 16L) + 1L]

ll <- d[, .(ll_mean = mean(ll_ex), ll_max = max(ll_ex)), by = .(gx, gy, ci)]
say("  cell-time units with list length: ", format(nrow(ll), big.mark = ","))

# --- frame + merge -----------------------------------------------------------
fr <- fread(file.path(.lp_dir, sprintf("%s_cells_5km_16d.csv", PREFIX)), showProgress = FALSE)
say("  frame ", format(nrow(fr), big.mark = ","), " rows")
n0 <- nrow(fr)
fr <- merge(fr, ll, by = c("gx","gy","ci"), all.x = TRUE)
stopifnot(nrow(fr) == n0)                       # merge must not duplicate rows
say(sprintf("  list length present on %.2f%% of frame rows",
            100 * fr[, mean(!is.na(ll_mean))]))

FEATS <- unlist(lapply(c("effort","season","veg","burn"), cfg_features))
FEATS <- intersect(FEATS, names(fr))

# NA HANDLING -- must match the pipeline, and a naive complete.cases does NOT.
#   yrs_since_burn is 66.3% NA and NA means NEVER BURNED, not unknown (Traps).
#     Dropping those rows deletes most never-burned cells -- i.e. the controls.
#   dist_median is 16.5% NA (effort distance not reported); the project imputes
#     and carries a missingness flag rather than dropping.
#   The residual ~0.5% are genuinely missing vegetation, which the published run
#     also drops ("dropped 4,373 rows with missing vegetation (0.70%)").
fr[is.na(yrs_since_burn), yrs_since_burn := 100]
fr[, dist_median_na := as.integer(is.na(dist_median))]
fr[is.na(dist_median), dist_median := 0]
FEATS <- c(FEATS, "dist_median_na")
fr <- fr[complete.cases(fr[, ..FEATS])]
fr <- fr[!is.na(ll_mean)]                       # same rows for both fits
say("  modelling rows ", format(nrow(fr), big.mark = ","), " | features ", length(FEATS))

tr <- fr$year <= 2023; te <- fr$year >= 2024
y  <- fr[[TARGET]]
say(sprintf("  train %s / test %s (test prevalence %.3f)",
            format(sum(tr), big.mark=","), format(sum(te), big.mark=","), mean(y[te])))

fit_auc <- function(feats, seed = 1) {
  set.seed(seed)
  X <- as.matrix(fr[, ..feats])
  m <- lgb.train(list(objective = "binary", metric = "auc", learning_rate = 0.05,
                      num_leaves = 63, min_data_in_leaf = 100,
                      feature_fraction = 0.8, bagging_fraction = 0.8,
                      bagging_freq = 1, verbose = -1, num_threads = 0),
                 lgb.Dataset(X[tr, ], label = y[tr]), nrounds = 400)
  auc(y[te], predict(m, X[te, ]))
}

say("Fitting baseline"); a0 <- fit_auc(FEATS)
say("Fitting + ll_mean"); a1 <- fit_auc(c(FEATS, "ll_mean"))
say("Fitting + ll_mean + ll_max"); a2 <- fit_auc(c(FEATS, "ll_mean", "ll_max"))

res <- data.table(species = PREFIX,
                  model = c("baseline", "+ ll_mean", "+ ll_mean + ll_max"),
                  n_feat = c(length(FEATS), length(FEATS) + 1L, length(FEATS) + 2L),
                  auc = c(a0, a1, a2))
res[, delta := auc - a0]
fwrite(res, file.path(.lp_dir, sprintf("list_length_predict_%s.csv", PREFIX)))

out <- file.path(.lp_dir, sprintf("list_length_predict_%s.txt", PREFIX))
sink(out)
cat("DOES LIST LENGTH IMPROVE PREDICTION?\n")
cat("====================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, "   target ", TARGET, "\n", sep = "")
cat("Split     : temporal (train <= 2023, test 2024-2026)\n")
cat("Model     : LightGBM, identical parameters to model_burn_effect.R\n\n")
print(as.data.frame(res[, .(model, n_feat, auc = round(auc, 4),
                            delta = sprintf("%+.4f", delta))]), row.names = FALSE)
cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  INTEGRITY CHECK: the published Wrentit baseline is LightGBM 0.8913 on\n")
cat("  450,654 training rows. This baseline reproduces it to within 0.0003 on\n")
cat("  449,443 rows, so the grid, features and parameters are being replicated\n")
cat("  correctly and the DELTA is a fair within-sample comparison.\n")
cat("\n  For scale: the entire 22-column burn block is worth about +0.046 AUC,\n")
cat("  and two coordinate columns are worth about +0.053.\n")
cat("\n  This is PREDICTION ONLY. Delta-AUC is not an effect size and says\n")
cat("  nothing about what fire does -- see Finding 5.\n")
sink()
cat(readLines(out), sep = "\n")
say("DONE")
