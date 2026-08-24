# =============================================================================
# Inference, in one place -- and a one-call API for running an analysis
#
# WHY THIS EXISTS
#   Two reasons.
#
#   1. DUPLICATION. Cluster-robust standard errors and the wild cluster
#      bootstrap were written out in four scripts (fire_cluster_se.R,
#      wild_cluster_boot.R, effort_placebo.R, and an inline test). Four copies
#      of the same maths is four chances to fix a bug in only three of them.
#
#   2. SPEED. Answering "what does this dose do to this species" took editing a
#      150-line script. It should take one line. analyze() and sweep() below are
#      the whole point of the modularity work: a question becomes a function
#      call, and a comparison becomes a table.
#
# WHICH INFERENCE TO USE -- this project learned the hard way
#   Use the WILD CLUSTER bootstrap. With ~225 fire clusters the percentile
#   bootstrap is the wrong instrument, and this project has the proof:
#   OSFL's dose interval FLIPPED VERDICT on rep count alone (B=200 said fail,
#   B=2000 said clear). Worse, percentile intervals got the WIDTH right but the
#   POSITION wrong -- Wrentit's centred on -0.0367 rather than its actual
#   -0.0298 because the resampled distribution is skewed, dragging the bound
#   past zero. boot_percentile() is kept only to reproduce old numbers.
#
# TYPICAL USE
#   source("panel_utils.R"); source("infer.R")
#
#   analyze("olsfly", "two_dose_sev")            # one species, one dose
#   sweep(dose = "frac_max")                     # all active species
#   sweep(species = "olsfly", dose = ALL_DOSES)  # one species, every dose
#   sweep()                                      # the full grid
# =============================================================================
suppressMessages(library(data.table))

.inf_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))

# analyze() calls cfg_default_prefix() and sweep() calls cfg_species(). In
# practice infer.R is always sourced after panel_utils.R, which provides both --
# but "in practice" is how cell_fire_lookup.R ended up calling a function it had
# never loaded. Source it here too, so infer.R stands on its own.
local({
  f <- file.path(.inf_dir, "project_config.R")
  if (file.exists(f) && !exists("cfg_default_prefix", mode = "function"))
    sys.source(f, envir = globalenv())
})

# --- core: demean once, then everything is small matrix algebra ---------------

#' Absorb cell + year FE and return the pieces every estimator needs.
#' Done once per panel; each bootstrap replication is then a tiny solve.
.prep <- function(panel, xcols, weights = TRUE) {
  z <- copy(panel)[, c("cell_id","year","rate","w", xcols), with = FALSE]
  for (cc in c("rate", xcols)) set(z, NULL, cc, as.numeric(z[[cc]]))
  if (!weights) set(z, NULL, "w", 1)
  for (i in seq_len(40)) for (cc in c("rate", xcols)) {
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
  }
  X <- as.matrix(z[, ..xcols])
  list(X = X, y = z$rate, w = z$w, Xw = X * z$w,
       Xi = solve(crossprod(X * z$w, X)), cl = panel$cl)
}

#' Cluster-robust variance matrix given residuals.
crve <- function(P, e) {
  S <- rowsum(P$Xw * e, P$cl, reorder = FALSE)
  P$Xi %*% crossprod(S) %*% P$Xi
}

.fit <- function(P, yy) as.vector(P$Xi %*% crossprod(P$Xw, yy))

#' Wild cluster bootstrap p-value for H0: beta = b0.
#' Null imposed (restricted residuals) + Rademacher signs per cluster --
#' the variant that does not under-reject with few clusters.
wcb_p <- function(P, ki, b_hat, se_hat, B = 999, b0 = 0) {
  y <- P$y; X <- P$X
  rest <- setdiff(seq_len(ncol(X)), ki)
  xi <- X[, ki]
  if (length(rest)) {
    Xr <- X[, rest, drop = FALSE]; Xrw <- Xr * P$w
    br <- as.vector(solve(crossprod(Xrw, Xr)) %*% crossprod(Xrw, y - b0 * xi))
    fitr <- as.vector(Xr %*% br) + b0 * xi
  } else fitr <- b0 * xi
  ur <- y - fitr
  t_obs <- (b_hat - b0) / se_hat
  ucl <- sort(unique(P$cl)); G <- length(ucl); cix <- match(P$cl, ucl)
  cnt <- 0L
  for (b in seq_len(B)) {
    ys <- fitr + sample(c(-1, 1), G, replace = TRUE)[cix] * ur
    bs <- .fit(P, ys)
    ts <- (bs[ki] - b0) / sqrt(crve(P, as.vector(ys - X %*% bs))[ki, ki])
    if (is.finite(ts) && abs(ts) >= abs(t_obs)) cnt <- cnt + 1L
  }
  (cnt + 1) / (B + 1)
}

#' Percentile bootstrap resampling CELLS. This is the project's ORIGINAL
#' inference and it is superseded -- kept because dose_burn_effect.R's published
#' intervals came from it and must stay reproducible. For a verdict use wcb_p().
#' New cell ids per draw, or the FE demeaning pools a cell with its own duplicate.
boot_cell <- function(dt, xcols, B = 200, iters = 15) {
  cells <- unique(dt$cell_id)
  out <- matrix(NA_real_, B, length(xcols), dimnames = list(NULL, xcols))
  for (b in seq_len(B)) {
    s <- data.table(cell_id = sample(cells, length(cells), replace = TRUE))[, newid := .I]
    bb <- merge(s, dt, by = "cell_id", allow.cartesian = TRUE)[, cell_id := newid]
    out[b, ] <- tryCatch(twfe(bb, xcols, iters), error = function(e) NA_real_)
  }
  out
}

#' Stratified percentile bootstrap: FIRES among treated, cells among controls.
#' Right unit, still the wrong instrument at ~225 clusters -- see header. The
#' informative output is the WIDTH RATIO against boot_cell(), not the verdict.
boot_fire <- function(dt, xcols, B = 200, iters = 15) {
  dtt <- dt[treat == 1]; dtc <- dt[treat == 0]
  fires <- unique(dtt$fire); ccells <- unique(dtc$cell_id)
  keep <- setdiff(names(dtt), "fire")        # same column set from both sides
  out <- matrix(NA_real_, B, length(xcols), dimnames = list(NULL, xcols))
  for (b in seq_len(B)) {
    sf <- data.table(fire = sample(fires, length(fires), replace = TRUE))[, draw := .I]
    tb <- merge(sf, dtt, by = "fire", allow.cartesian = TRUE)
    tb[, cell_id := as.integer(draw * 100000L + cell_id)]
    sc <- data.table(cell_id = sample(ccells, length(ccells), replace = TRUE))[, draw := .I]
    cb <- merge(sc, dtc, by = "cell_id", allow.cartesian = TRUE)
    cb[, cell_id := as.integer(1000000000L + draw)]
    out[b, ] <- tryCatch(twfe(rbind(tb[, ..keep], cb[, ..keep]), xcols, iters),
                         error = function(e) NA_real_)
  }
  out
}

boot_ci <- function(m, k) quantile(m[, k], c(0.025, 0.975), na.rm = TRUE)

# =============================================================================
# THE ONE-CALL API
# =============================================================================

#' Run one species x one burn measure and return a tidy result row per dose term.
#'
#' @param species  prefix from species_registry.csv, e.g. "olsfly"
#' @param dose     name from treatment_registry.csv, e.g. "two_dose_sev"
#' @param B        wild cluster replications. 999 for exploring, 9999 to quote.
#' @param weights  checklist-weighted (TRUE, matches the published spec)
#' @param quiet    suppress progress
#' @return data.table: species, dose, term, estimate, se, t, p, n_cells, n_fires
analyze <- function(species = NULL, dose = "frac_max", B = 999,
                    weights = TRUE, quiet = FALSE) {
  if (is.null(species)) species <- cfg_default_prefix()
  if (!quiet) message(sprintf("  analyze: %-8s x %-14s", species, dose))
  p  <- attach_fire(build_panel(species, dose = dose), species)
  tm <- attr(p, "dose_terms")
  xc <- c("post_d", paste0("post_", tm))
  P  <- .prep(p, xc, weights = weights)
  b  <- .fit(P, P$y)
  V  <- crve(P, as.vector(P$y - P$X %*% b))
  base <- p[treat == 1 & post == 0, sum(det)/sum(w)]
  out <- rbindlist(lapply(seq_along(tm), function(k) {
    ki <- match(paste0("post_", tm[k]), xc)
    se <- sqrt(V[ki, ki])
    data.table(species = species, dose = dose, term = tm[k],
               estimate = b[ki], se = se, t = b[ki]/se,
               p = wcb_p(P, ki, b[ki], se, B = B),
               pct_of_base = 100 * b[ki] / base,
               n_cells = uniqueN(p[treat == 1]$cell_id),
               n_fires = uniqueN(p[treat == 1]$fire))
  }))
  # added here, not in sweep(), so analyze() and sweep() return the same shape
  out[, sig := fifelse(p < 0.05, "*", "")]
  out[]
}

ALL_SPECIES <- function() {
  s <- cfg_species(active_only = TRUE)$prefix
  s[vapply(s, function(x) file.exists(file.path(.inf_dir,
       sprintf("%s_cells_5km_16d.csv", x))), logical(1))]
}
ALL_DOSES <- function() cfg_treatments()$name

#' Run a grid of species x burn measures and return one table.
#' Defaults to every species with a built frame and every ACTIVE burn measure.
sweep <- function(species = NULL, dose = NULL, B = 999, weights = TRUE) {
  if (is.null(species)) species <- ALL_SPECIES()
  if (is.null(dose))    dose    <- cfg_treatments(active_only = TRUE)$name
  grid <- CJ(species = species, dose = dose, sorted = FALSE)
  message(sprintf("sweep: %d species x %d dose(s) = %d fits, B=%d",
                  length(species), length(dose), nrow(grid), B))
  rbindlist(lapply(seq_len(nrow(grid)), function(i)
    tryCatch(analyze(grid$species[i], grid$dose[i], B = B, weights = weights),
             error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL })))
}

#' Print a sweep result compactly.
show_sweep <- function(r) {
  print(as.data.frame(r[, .(species, dose, term,
                            est = round(estimate, 4), se = round(se, 4),
                            t = round(t, 2), p = round(p, 4),
                            pct = round(pct_of_base, 1), fires = n_fires, sig)]),
        row.names = FALSE)
  invisible(r)
}
