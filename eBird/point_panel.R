# =============================================================================
# STEP 2 -- the SAME SPOT before and after. The fix for the effort placebo.
#
# WHY THIS EXISTS -- the one criticism that currently lands
#   effort_placebo.R found that after a fire, birders visit 43% FEWER DISTINCT
#   LOCATIONS per 5 km cell (p = 0.0001), and file 7% fewer checklists
#   (p = 0.015). A cell fixed effect compares a cell to itself, but it CANNOT
#   absorb a change in WHICH SPOTS INSIDE IT are sampled. If post-fire access
#   concentrates on roadsides, the habitat being sampled changes and every
#   species estimate in this project is contaminated by an unknown amount.
#
#   Comparing the SAME LOCATION before and after removes that channel entirely.
#   This is the only design in the project that does.
#
# SECOND BENEFIT
#   Distance stops being blurred. A 5 km cell centre is accurate to about half a
#   cell; a birding spot is accurate to the 500 m snap.
#
# THE SNAP
#   4.4M checklists sit on 498,566 distinct coordinates -- NOT the 36,519 that
#   CLAUDE.md claimed, which is distinct NAMED localities. point_fire_distance.R
#   snaps to a 500 m grid (px, py), collapsing them to ~156k. Two birders 400 m
#   apart are the same place relative to a fire kilometres away.
#
# UNIT
#   (spot x year). Location fixed effects + year fixed effects, clustered by the
#   causing FIRE -- spots inside one perimeter share one ignition, exactly as
#   cells do. Outcome is det_rate = detections / checklists at that spot-year.
#
# WHAT IT STILL CANNOT DO
#   Detection is not occupancy. And a spot observed only before or only after
#   its fire contributes nothing to a within-spot comparison, so the usable
#   sample is much smaller than the checklist count suggests -- see MIN_YEARS.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

.pp_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  f <- file.path(.pp_dir, "panel_utils.R")
  if (file.exists(f) && !exists("twfe", mode = "function")) sys.source(f, envir = globalenv())
})

#' Spot x year panel with proximity treatment.
#'
#' @param prefix     species prefix
#' @param R          radius km; a fire inside R treats the spot
#' @param MIN_YEARS  a spot must be birded in at least this many distinct years,
#'                   or it can carry no before/after contrast at all
#' @param decay      "linear" (1 - d/R) or "exp"
build_point_panel <- function(prefix = NULL, R = 10, MIN_YEARS = 2L,
                              decay = c("linear","exp")) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  decay <- match.arg(decay)

  jf  <- file.path(.pp_dir, sprintf("%s_2015_2026_joined.csv", prefix))
  mf  <- file.path(.pp_dir, "point_grid_map.csv")
  df  <- file.path(.pp_dir, "point_fire_distance.csv")
  for (f in c(jf, mf, df)) if (!file.exists(f)) stop("missing input: ", basename(f))

  message("  reading checklists")
  d <- fread(jf, select = c("latitude","longitude","year","presence"), showProgress = FALSE)
  m <- fread(mf, showProgress = FALSE)
  d <- merge(d, m, by = c("latitude","longitude"))
  d[, c("latitude","longitude") := NULL]

  # spot-year outcome. `presence` is the zero-filled detection flag; det_rate is
  # detections / checklists, so it is already effort-normalised for VOLUME --
  # which is the thing cell FE could absorb. What it could not absorb, and what
  # this design fixes, is which SPOT was visited.
  sy <- d[, .(det = sum(presence), chk = .N), by = .(px, py, year)]
  keep <- sy[, .(ny = uniqueN(year)), by = .(px, py)][ny >= MIN_YEARS]
  sy <- merge(sy, keep[, .(px, py)], by = c("px","py"))
  message(sprintf("  spot-years %s on %s spots (>= %d years birded)",
                  format(nrow(sy), big.mark=","),
                  format(nrow(keep), big.mark=","), MIN_YEARS))

  message("  attaching fire distances")
  D <- fread(df, showProgress = FALSE)
  D <- D[dist_km <= R]
  obs <- sy[, .(y_min = min(year), y_max = max(year)), by = .(px, py)]
  D <- merge(D[, .(px, py, fire_event_id, fire_year, dist_km)], obs, by = c("px","py"))
  D <- D[fire_year > y_min & fire_year <= y_max]     # must have a before AND after
  if (!nrow(D)) stop("no qualifying spot-fire pairs at R = ", R)

  setorder(D, px, py, fire_year, dist_km)
  ev <- D[, .(ev_year = fire_year[1], dist_km = dist_km[1],
              fire = fire_event_id[1]), by = .(px, py)]
  ev[, dose := if (decay == "linear") pmax(0, 1 - dist_km / R)
              else exp(-dist_km / (R / 3))]

  ctl <- fsetdiff(unique(sy[, .(px, py)]), unique(ev[, .(px, py)]))
  set.seed(1)
  ctl[, `:=`(ev_year = sample(ev$ev_year, .N, replace = TRUE),
             dist_km = NA_real_, fire = NA_character_, dose = 0)]

  p <- merge(sy, rbind(ev, ctl), by = c("px","py"))
  p[, `:=`(treat = as.integer(!is.na(fire)), post = as.integer(year >= ev_year))]
  p[, `:=`(rate = det / chk, w = as.numeric(chk),
           post_d = as.numeric(post), post_dose = as.numeric(post) * dose)]
  # twfe()/.prep() key the unit fixed effect off `cell_id`; here the unit is a
  # SPOT, so give it that name rather than duplicating the estimator
  p[, cell_id := .GRP, by = .(px, py)]
  p[is.na(fire), fire := paste0("CTRL_", cell_id)]
  p[, cl := as.integer(factor(fire))]

  message(sprintf("  treated spots %s / %s fires | control spots %s",
                  format(uniqueN(p[treat==1]$cell_id), big.mark=","),
                  format(uniqueN(p[treat==1]$fire), big.mark=","),
                  format(uniqueN(p[treat==0]$cell_id), big.mark=",")))
  setattr(p, "dose_terms", "dose"); setattr(p, "radius_km", R)
  p[]
}

#' Fit it. Same estimator and inference as everything else in the project.
analyze_point <- function(prefix = NULL, R = 10, MIN_YEARS = 2L, B = 999,
                          decay = "linear") {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  p  <- build_point_panel(prefix, R = R, MIN_YEARS = MIN_YEARS, decay = decay)
  xc <- c("post_d","post_dose")
  P  <- .prep(p, xc)
  b  <- .fit(P, P$y)
  se <- sqrt(crve(P, as.vector(P$y - P$X %*% b))[2,2])
  base <- p[treat == 1 & post == 0, sum(det)/sum(w)]
  data.table(species = prefix, grain = "point", radius_km = R,
             estimate = b[2], se = se, t = b[2]/se,
             p = wcb_p(P, 2L, b[2], se, B = B),
             pct_of_base = 100 * b[2] / base,
             n_spots = uniqueN(p[treat==1]$cell_id),
             n_fires = uniqueN(p[treat==1]$fire))[]
}
