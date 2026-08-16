# =============================================================================
# LSTM and CNN -- the two architectures from the spec, finally runnable.
#
# WHY THEY WERE NEVER BUILT
#   `torch` and `luz` were installed but libtorch (the compute backend) was not,
#   so nothing could run. Installed 2026-08-01; verified with a live tensor.
#
# WHAT EACH ONE IS FOR -- and it is NOT a better effect estimate
#   Neither of these is a causal estimator. The project has already shown that
#   two coordinate columns add +0.053 AUC while the entire 22-column burn block
#   adds +0.002 once geography is present. A more flexible model will mostly
#   learn GEOGRAPHY better. Expect these to improve forecasting and to say
#   nothing directly about what fire does. They are answering the predictive
#   question, and delta-AUC is still not an effect size.
#
#   The LSTM does have one genuinely new capability: the existing models treat
#   each cell-time as independent, with lags bolted on as columns. An LSTM reads
#   the actual SEQUENCE, so it can represent a recovery trajectory -- which is
#   the shape Finding 4 says matters (Wrentit crashes then rebounds). That is the
#   scientifically interesting part.
#
# THE CNN -- a deliberate deviation from the spec, stated plainly
#   The spec envisaged CNN patches cut from the raster imagery. That is not what
#   this does. Raster patch extraction on this machine is a serious build: the 42
#   MTBS rasters have 42 different origins with no common lattice, and the
#   project has already OOMed at 24 GB materialising their coordinates.
#
#   Instead this builds patches from the CELL GRID itself: for each focal cell, a
#   5x5 window of its neighbours' burn and vegetation values, assembled by (gx,gy)
#   offset from data already on disk. It is a real convolution over spatial
#   context -- "what is happening around this cell" -- at 5 km resolution rather
#   than 30 m. Cheap, safe, and honest about what it is.
#
# ENV
#   DL_SPECIES   prefixes to run (default: one per guild + the controls)
#   DL_EPOCHS    training epochs (default 15)
#   DL_SEQ       LSTM sequence length in 16-day steps (default 8 = ~4 months)
# =============================================================================
suppressMessages({ library(data.table); library(torch) })
setDTthreads(0)
torch_manual_seed(1); set.seed(1)
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))
EPOCHS <- as.integer(Sys.getenv("DL_EPOCHS", "15"))
SEQ    <- as.integer(Sys.getenv("DL_SEQ", "8"))
out_sum <- file.path(dir_, "deep_models_summary.txt")

auc <- function(y, p) {
  r <- rank(p); n1 <- as.numeric(sum(y == 1)); n0 <- as.numeric(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

FEAT <- c(cfg_features("effort"), cfg_features("season"),
          cfg_features("veg"), cfg_features("burn"))

# =============================================================================
# LSTM -- sequence of 16-day steps per cell
# =============================================================================
run_lstm <- function(pre) {
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  d <- fread(f, showProgress = FALSE)
  d[is.na(yrs_since_burn), yrs_since_burn := 100]
  d[, dist_median_na := as.integer(is.na(dist_median))]
  d[is.na(dist_median), dist_median := 0]
  cols <- intersect(c(FEAT, "dist_median_na"), names(d))
  d <- d[complete.cases(d[, ..cols])]
  setorder(d, cell_id, ci)

  # The panel is IRREGULAR -- a cell is birded sporadically, not every step. The
  # LSTM reads observed steps in order and is told how big each gap was, rather
  # than pretending the series is evenly spaced.
  d[, gap := c(0, diff(ci)), by = cell_id]
  cols <- c(cols, "gap")

  X <- as.matrix(d[, ..cols])
  X <- scale(X); X[!is.finite(X)] <- 0
  y <- d$any_detection
  cid <- d$cell_id; yr <- d$year

  # many-to-one: the previous SEQ observed steps of a cell predict detection now
  idx <- which(seq_len(nrow(d)) >= SEQ &
               cid == data.table::shift(cid, SEQ - 1L, type = "lag"))
  if (length(idx) < 2000) { say("  too few sequences for ", pre); return(NULL) }
  tr_i <- idx[yr[idx] <= 2023]; te_i <- idx[yr[idx] >= 2024]
  if (!length(te_i)) { say("  no test rows for ", pre); return(NULL) }

  mk <- function(ii) {
    a <- array(0, dim = c(length(ii), SEQ, ncol(X)))
    for (s in seq_len(SEQ)) a[, s, ] <- X[ii - SEQ + s, , drop = FALSE]
    a
  }
  say("  LSTM tensors: train ", length(tr_i), " test ", length(te_i))
  xtr <- torch_tensor(mk(tr_i), dtype = torch_float())
  ytr <- torch_tensor(matrix(y[tr_i], ncol = 1), dtype = torch_float())
  xte <- torch_tensor(mk(te_i), dtype = torch_float())

  net <- nn_module(
    initialize = function(nf) {
      self$lstm <- nn_lstm(nf, 48, batch_first = TRUE)
      self$fc   <- nn_linear(48, 1)
    },
    forward = function(x) {
      o <- self$lstm(x)[[1]]
      self$fc(o[ , dim(o)[2], ])
    })(ncol(X))

  opt <- optim_adam(net$parameters, lr = 0.005)
  lossf <- nn_bce_with_logits_loss()
  nb <- 256L
  for (e in seq_len(EPOCHS)) {
    net$train(); perm <- sample(length(tr_i))
    for (b in seq(1, length(perm), by = nb)) {
      j <- perm[b:min(b + nb - 1L, length(perm))]
      opt$zero_grad()
      l <- lossf(net(xtr[j, , ]), ytr[j, , drop = FALSE])
      l$backward(); opt$step()
    }
  }
  net$eval()
  p <- as.numeric(torch_sigmoid(with_no_grad(net(xte))))
  data.table(model = "LSTM", auc = auc(y[te_i], p), n_test = length(te_i))
}

# =============================================================================
# CNN -- 5x5 neighbourhood of CELLS around each focal cell
# =============================================================================
run_cnn <- function(pre) {
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  ch <- c("burn_frac_max","burn_sev_wt","ndvi_l0","evi_l0")     # 4 channels
  d <- fread(f, select = c("cell_id","gx","gy","year","ci","any_detection",
                           "n_checklists", ch), showProgress = FALSE)
  d <- d[complete.cases(d)]
  # one spatial snapshot per (gx, gy, year) keeps the tensor manageable
  s <- d[, lapply(.SD, mean), by = .(gx, gy, year), .SDcols = ch]
  lab <- d[, .(y = as.integer(max(any_detection))), by = .(gx, gy, year)]
  s <- merge(s, lab, by = c("gx","gy","year"))

  key <- s[, paste(gx, gy, year)]
  look <- new.env(hash = TRUE, size = nrow(s) * 2L)
  M <- as.matrix(s[, ..ch])
  for (i in seq_len(nrow(s))) assign(key[i], i, envir = look)

  R <- 2L; W <- 2L * R + 1L
  say("  CNN patches: ", nrow(s), " focal cells, ", W, "x", W, "x", length(ch))
  A <- array(0, dim = c(nrow(s), length(ch), W, W))
  gxv <- s$gx; gyv <- s$gy; yrv <- s$year
  for (i in seq_len(nrow(s))) {
    for (dx in -R:R) for (dy in -R:R) {
      k <- paste(gxv[i] + dx, gyv[i] + dy, yrv[i])
      j <- mget(k, envir = look, ifnotfound = list(NULL))[[1]]
      if (!is.null(j)) A[i, , dx + R + 1L, dy + R + 1L] <- M[j, ]
    }
  }
  for (c_ in seq_len(length(ch))) {
    v <- A[, c_, , ]; mu <- mean(v); sd_ <- stats::sd(v); if (sd_ == 0) sd_ <- 1
    A[, c_, , ] <- (v - mu) / sd_
  }

  tr <- which(yrv <= 2023); te <- which(yrv >= 2024)
  if (length(te) < 500) { say("  too few test cells for ", pre); return(NULL) }
  xtr <- torch_tensor(A[tr, , , , drop = FALSE], dtype = torch_float())
  ytr <- torch_tensor(matrix(s$y[tr], ncol = 1), dtype = torch_float())
  xte <- torch_tensor(A[te, , , , drop = FALSE], dtype = torch_float())

  net <- nn_module(
    initialize = function(nc) {
      self$c1 <- nn_conv2d(nc, 32, 3, padding = 1)
      self$c2 <- nn_conv2d(32, 64, 3, padding = 1)
      self$fc1 <- nn_linear(64 * W * W, 64); self$fc2 <- nn_linear(64, 1)
    },
    forward = function(x) {
      x <- nnf_relu(self$c1(x)); x <- nnf_relu(self$c2(x))
      x <- torch_flatten(x, start_dim = 2)
      self$fc2(nnf_relu(self$fc1(x)))
    })(length(ch))

  opt <- optim_adam(net$parameters, lr = 0.002)
  lossf <- nn_bce_with_logits_loss(); nb <- 128L
  for (e in seq_len(EPOCHS)) {
    net$train(); perm <- sample(length(tr))
    for (b in seq(1, length(perm), by = nb)) {
      j <- perm[b:min(b + nb - 1L, length(perm))]
      opt$zero_grad()
      l <- lossf(net(xtr[j, , , ]), ytr[j, , drop = FALSE])
      l$backward(); opt$step()
    }
  }
  net$eval()
  p <- as.numeric(torch_sigmoid(with_no_grad(net(xte))))
  data.table(model = "CNN", auc = auc(s$y[te], p), n_test = length(te))
}

# =============================================================================
sp_all <- cfg_species(active_only = FALSE)
sp_all <- sp_all[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]
want <- Sys.getenv("DL_SPECIES", "")
prefs <- if (nzchar(want)) strsplit(want, "[ ,]+")[[1]] else
  sp_all[, .SD[1], by = guild]$prefix          # one per guild by default
sp_all <- sp_all[prefix %in% prefs]

res <- rbindlist(lapply(seq_len(nrow(sp_all)), function(i) {
  pre <- sp_all$prefix[i]; say("=== ", sp_all$common_name[i], " (", pre, ")")
  out <- list()
  for (fn in list(run_lstm, run_cnn)) {
    r <- tryCatch(fn(pre), error = function(e) { say("  FAILED: ",
                    substr(conditionMessage(e), 1, 70)); NULL })
    if (!is.null(r)) { r[, `:=`(species = sp_all$common_name[i], prefix = pre,
                                guild = sp_all$guild[i])]; out[[length(out)+1]] <- r }
  }
  rbindlist(out)
}), fill = TRUE)

if (nrow(res)) fwrite(res, file.path(dir_, "deep_models_results.csv"))

sink(out_sum)
cat("LSTM and CNN -- the two architectures from the spec\n")
cat("==================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Split     : temporal (train <= 2023, test 2024-2026), same as the other models\n")
cat("Epochs    : ", EPOCHS, "   LSTM sequence length: ", SEQ, " x 16-day steps\n\n", sep = "")
if (!nrow(res)) cat("  no results\n") else
print(as.data.frame(res[order(model, -auc), .(species, guild, model,
      auc = round(auc, 4), n_test)]), row.names = FALSE)

cat("\n-- BENCHMARK --------------------------------------------------------\n")
cat("  Wrentit, same split, existing models:\n")
cat("    LightGBM 0.8913 | XGBoost 0.8880 | gbm 0.8680 | logistic 0.8311\n")
cat("    LightGBM + coordinates 0.9440\n")
cat("\n  Beating 0.8913 would be a genuine predictive gain. Not beating it is\n")
cat("  the more likely outcome and is itself worth reporting: gradient-boosted\n")
cat("  trees are hard to beat on tabular data of this size.\n")

cat("\n-- WHAT THESE DO NOT SHOW -------------------------------------------\n")
cat("  * Any causal effect. Neither architecture is an effect estimator, and\n")
cat("    delta-AUC remains not an effect size.\n")
cat("  * The CNN uses a 5x5 window of 5 km CELLS, not raster imagery patches.\n")
cat("    That is a deliberate deviation -- see the header for why.\n")
cat("  * The LSTM reads irregularly spaced observations with a gap feature; it\n")
cat("    is not a continuous-time model.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
