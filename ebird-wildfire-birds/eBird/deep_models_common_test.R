# =============================================================================
# COMMON TEST SET -- make the CNN, LSTM and tree AUCs comparable
#
# THE PROBLEM
#   deep_models.R scores the CNN on 24,240 units and the LSTM on 160,440. Those
#   are different quantities on different rows:
#     LSTM  cell-time (16-day) rows, target = any_detection in that window
#     CNN   (gx, gy, year) snapshots, target = max(any_detection) over the year
#   A yearly max over a cell is a far easier target than detection in a given
#   16-day window, so 0.8559 vs 0.8841 compares nothing. Both numbers are also
#   printed under a LightGBM benchmark computed on a third row set.
#
# WHAT THIS DOES
#   Builds ONE base table, fits all three models from it, and scores them on the
#   SAME test rows against the SAME target (any_detection at cell-time).
#
#   The CNN produces one prediction per cell-year, so it is broadcast to every
#   16-day row in that cell-year. That is the honest way to score a spatial
#   model on a temporal task: it is allowed to say what it knows, and it is not
#   allowed to vary within a year because it cannot.
#
# THE SECOND COMPARISON, WHICH ISOLATES WHY
#   Scoring the CNN at cell-time penalises it for two separate things: a weaker
#   architecture, and no temporal resolution. To separate them, every model is
#   ALSO scored at cell-year grain (predictions averaged within cell-year,
#   target = max over the year). If the CNN closes the gap there, its deficit is
#   temporal resolution rather than the convolution.
#
# ENV  DCT_PREFIX (default cfg_default_prefix()) · DCT_EPOCHS (15) · DCT_SEQ (8)
# =============================================================================
suppressMessages({ library(data.table); library(torch); library(lightgbm) })
setDTthreads(0)
torch_manual_seed(1); set.seed(1)

.dc_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({ p <- file.path(.dc_dir, "project_config.R"); if (file.exists(p)) sys.source(p, envir = globalenv()) })
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

PREFIX <- Sys.getenv("DCT_PREFIX", ""); if (!nzchar(PREFIX)) PREFIX <- cfg_default_prefix()
EPOCHS <- as.integer(Sys.getenv("DCT_EPOCHS", "15"))
SEQ    <- as.integer(Sys.getenv("DCT_SEQ", "8"))
TARGET <- "any_detection"

auc <- function(y, p) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# --- one base table all three models are built from --------------------------
FEAT <- unlist(lapply(c("effort","season","veg","burn"), cfg_features))
CH   <- c("burn_frac_max","burn_sev_wt","ndvi_l0","evi_l0")   # CNN channels

d <- fread(file.path(.dc_dir, sprintf("%s_cells_5km_16d.csv", PREFIX)), showProgress = FALSE)
d[is.na(yrs_since_burn), yrs_since_burn := 100]
d[, dist_median_na := as.integer(is.na(dist_median))]
d[is.na(dist_median), dist_median := 0]
FEAT <- c(intersect(FEAT, names(d)), "dist_median_na")
keep <- unique(c(FEAT, CH))
d <- d[complete.cases(d[, ..keep])]
setorder(d, cell_id, ci)
d[, row := .I]
say("base table ", format(nrow(d), big.mark = ","), " rows | features ", length(FEAT))

y  <- d[[TARGET]]
tr <- d$year <= 2023
te <- d$year >= 2024

# --- LSTM --------------------------------------------------------------------
say("LSTM")
d[, gap := c(0, diff(ci)), by = cell_id]
LF <- c(FEAT, "gap")
X  <- as.matrix(d[, ..LF]); X <- scale(X); X[!is.finite(X)] <- 0
cid <- d$cell_id
idx <- which(seq_len(nrow(d)) >= SEQ & cid == data.table::shift(cid, SEQ - 1L, type = "lag"))
tr_i <- idx[tr[idx]]; te_i <- idx[te[idx]]
say(sprintf("  sequences: train %s test %s", format(length(tr_i), big.mark=","),
            format(length(te_i), big.mark=",")))
mk <- function(ii) { a <- array(0, dim = c(length(ii), SEQ, ncol(X)))
                     for (s in seq_len(SEQ)) a[, s, ] <- X[ii - SEQ + s, , drop = FALSE]; a }
xtr <- torch_tensor(mk(tr_i), dtype = torch_float())
ytr <- torch_tensor(matrix(y[tr_i], ncol = 1), dtype = torch_float())
net <- nn_module(
  initialize = function(nf) { self$lstm <- nn_lstm(nf, 48, batch_first = TRUE)
                              self$fc <- nn_linear(48, 1) },
  forward = function(x) { o <- self$lstm(x)[[1]]; self$fc(o[, dim(o)[2], ]) })(ncol(X))
opt <- optim_adam(net$parameters, lr = 0.005); lossf <- nn_bce_with_logits_loss(); nb <- 256L
for (e in seq_len(EPOCHS)) { net$train(); perm <- sample(length(tr_i))
  for (b in seq(1, length(perm), by = nb)) { j <- perm[b:min(b+nb-1L, length(perm))]
    opt$zero_grad(); l <- lossf(net(xtr[j,,]), ytr[j,,drop=FALSE]); l$backward(); opt$step() } }
net$eval()
p_lstm <- rep(NA_real_, nrow(d))
p_lstm[te_i] <- as.numeric(torch_sigmoid(with_no_grad(net(torch_tensor(mk(te_i), dtype = torch_float())))))

# --- CNN: one prediction per (gx, gy, year), then broadcast -------------------
say("CNN")
s <- d[, lapply(.SD, mean), by = .(gx, gy, year), .SDcols = CH]
lab <- d[, .(yy = as.integer(max(get(TARGET)))), by = .(gx, gy, year)]
s <- merge(s, lab, by = c("gx","gy","year"))
key <- s[, paste(gx, gy, year)]
look <- new.env(hash = TRUE, size = nrow(s) * 2L)
M <- as.matrix(s[, ..CH])
for (i in seq_len(nrow(s))) assign(key[i], i, envir = look)
R <- 2L; W <- 2L*R + 1L
A <- array(0, dim = c(nrow(s), length(CH), W, W))
gxv <- s$gx; gyv <- s$gy; yrv <- s$year
for (i in seq_len(nrow(s))) for (dx in -R:R) for (dy in -R:R) {
  j <- mget(paste(gxv[i]+dx, gyv[i]+dy, yrv[i]), envir = look, ifnotfound = list(NULL))[[1]]
  if (!is.null(j)) A[i, , dx+R+1L, dy+R+1L] <- M[j, ] }
for (c_ in seq_along(CH)) { v <- A[, c_, , ]; sd_ <- stats::sd(v); if (sd_ == 0) sd_ <- 1
                            A[, c_, , ] <- (v - mean(v)) / sd_ }
str_ <- which(yrv <= 2023)
xtr2 <- torch_tensor(A[str_,,,,drop=FALSE], dtype = torch_float())
ytr2 <- torch_tensor(matrix(s$yy[str_], ncol = 1), dtype = torch_float())
cnn <- nn_module(
  initialize = function(nc) { self$c1 <- nn_conv2d(nc, 32, 3, padding = 1)
    self$c2 <- nn_conv2d(32, 64, 3, padding = 1)
    self$f1 <- nn_linear(64*W*W, 64); self$f2 <- nn_linear(64, 1) },
  forward = function(x) { x <- nnf_relu(self$c1(x)); x <- nnf_relu(self$c2(x))
    x <- torch_flatten(x, start_dim = 2); self$f2(nnf_relu(self$f1(x))) })(length(CH))
opt2 <- optim_adam(cnn$parameters, lr = 0.002); nb2 <- 128L
for (e in seq_len(EPOCHS)) { cnn$train(); perm <- sample(length(str_))
  for (b in seq(1, length(perm), by = nb2)) { j <- perm[b:min(b+nb2-1L, length(perm))]
    opt2$zero_grad(); l <- lossf(cnn(xtr2[j,,,]), ytr2[j,,drop=FALSE]); l$backward(); opt2$step() } }
cnn$eval()
s[, p := as.numeric(torch_sigmoid(with_no_grad(cnn(torch_tensor(A, dtype = torch_float())))))]
d <- merge(d, s[, .(gx, gy, year, p_cnn = p)], by = c("gx","gy","year"), all.x = TRUE)
setorder(d, row)                                   # merge reorders; restore
p_cnn <- d$p_cnn

# --- LightGBM reference ------------------------------------------------------
say("LightGBM")
XL <- as.matrix(d[, ..FEAT])
set.seed(1)
gbm_ <- lgb.train(list(objective="binary", metric="auc", learning_rate=0.05,
                       num_leaves=63, min_data_in_leaf=100, feature_fraction=0.8,
                       bagging_fraction=0.8, bagging_freq=1, verbose=-1, num_threads=0),
                  lgb.Dataset(XL[tr, ], label = y[tr]), nrounds = 400)
p_gbm <- predict(gbm_, XL)

# --- score all three on the SAME rows ----------------------------------------
common <- which(te & !is.na(p_lstm) & !is.na(p_cnn))
say("common test rows: ", format(length(common), big.mark = ","))
cell_time <- data.table(
  model = c("LightGBM","LSTM","CNN"),
  auc = c(auc(y[common], p_gbm[common]), auc(y[common], p_lstm[common]),
          auc(y[common], p_cnn[common])))

# cell-year grain: average predictions within cell-year, target = max
zz <- data.table(gx = d$gx[common], gy = d$gy[common], year = d$year[common],
                 yy = y[common], gbm = p_gbm[common], lstm = p_lstm[common],
                 cnn = p_cnn[common])[, .(yy = max(yy), gbm = mean(gbm),
                                          lstm = mean(lstm), cnn = mean(cnn)),
                                      by = .(gx, gy, year)]
cell_year <- data.table(
  model = c("LightGBM","LSTM","CNN"),
  auc = c(auc(zz$yy, zz$gbm), auc(zz$yy, zz$lstm), auc(zz$yy, zz$cnn)))

out <- file.path(.dc_dir, sprintf("deep_models_common_test_%s.txt", PREFIX))
sink(out)
cat("DEEP MODELS ON A COMMON TEST SET\n=================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, "  epochs ", EPOCHS, "  LSTM seq ", SEQ, "\n", sep = "")
cat("\nEvery model is fitted on train (<=2023) and scored on the SAME test rows.\n")
cat("The CNN predicts one value per cell-year and is broadcast to the 16-day\n")
cat("rows within it -- it may say what it knows, but cannot vary within a year.\n\n")
cat("-- CELL-TIME grain (16-day windows, target = any_detection) ----------\n")
cat("   common test rows: ", format(length(common), big.mark=","), "\n\n", sep="")
print(as.data.frame(cell_time[, .(model, auc = round(auc, 4))]), row.names = FALSE)
cat("\n-- CELL-YEAR grain (target = max over the year) ----------------------\n")
cat("   units: ", format(nrow(zz), big.mark=","), "\n\n", sep="")
print(as.data.frame(cell_year[, .(model, auc = round(auc, 4))]), row.names = FALSE)
cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  The cell-time table is the comparable one: same rows, same target.\n")
cat("  The cell-year table separates architecture from temporal resolution.\n")
cat("  If the CNN closes the gap at cell-year grain, its cell-time deficit is\n")
cat("  the absence of within-year variation rather than the convolution.\n")
cat("\n  Superseded: deep_models.R scored the CNN on 24,240 units and the LSTM\n")
cat("  on 160,440, against different targets. Those AUCs are not comparable\n")
cat("  and should not be quoted side by side.\n")
sink()
cat(readLines(out), sep = "\n")
say("DONE")
