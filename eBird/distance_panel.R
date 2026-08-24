# =============================================================================
# STEP 1 -- treatment by PROXIMITY to fire, not by being burned.
#
# WHY THIS IS A SEPARATE FILE
#   build_panel() defines "treated" as `ever_burned` flipping 0 -> 1 inside the
#   window. Every published number in this project flows through it. A distance
#   treatment changes WHO IS TREATED -- not just the dose value -- so bolting it
#   onto build_panel() would put the headline estimates one bug away from moving.
#   This is a parallel builder that reuses the shared ESTIMATOR (twfe, crve,
#   wcb_p) but derives treatment independently. panel_utils.R is untouched.
#
# THE QUESTION
#   Today a cell is treated only if fire burned it: 1,009 cells behind 225 fires.
#   Verified at checklist grain, 100% of the 206,379 fire-matched rows sit INSIDE
#   a perimeter, so a cell 500 m from a 100,000-acre scar is encoded identically
#   to one 200 km away. step1_gate.R measured what proximity buys: at 10 km,
#   405 distinct fires qualify (+80%), which clears the ~285 the Wrentit result
#   needs to reach t = 1.96.
#
#   That is a COUNT, not an effect. This file produces the effect.
#
# THE DOSE
#   dose = 1 - dist_km / R, clipped at 0. So dose = 1 for a fire at the cell and
#   decays linearly to 0 at the radius. Bounded 0-1 like frac_max, so the
#   coefficient is again "effect at full dose" and is directly comparable to the
#   published -0.0298. `decay = "exp"` gives exp(-dist/(R/3)) instead, for
#   checking that the answer is not an artifact of the linear shape.
#
# WHAT IT CANNOT SHOW
#   That distant fires matter. A larger radius admits more fires but weaker
#   treatment, trading bias for variance. Read the radius sweep as a whole.
# =============================================================================
suppressMessages({ library(data.table) })

.dp_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  f <- file.path(.dp_dir, "panel_utils.R")
  if (file.exists(f) && !exists("twfe", mode = "function")) sys.source(f, envir = globalenv())
})

#' Cell-year panel with treatment defined by proximity to fire.
#'
#' @param prefix species prefix
#' @param R      radius in km. Cells with a qualifying fire inside R are treated.
#' @param decay  "linear" (default) or "exp"
#' @param binary TRUE -> dose is 1/0 rather than distance-graded
#' @return data.table like build_panel()'s, plus `fire` (the qualifying fire id)
#'         and `dist_km`. attr "dose_terms" = "dose".
build_dist_panel <- function(prefix = NULL, R = 10, decay = c("linear","exp"),
                             binary = FALSE) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  decay <- match.arg(decay)

  fdf <- file.path(.dp_dir, "cell_fire_distance.csv")
  if (!file.exists(fdf))
    stop("cell_fire_distance.csv not found -- run fire_distance_lookup.R first")
  fr <- file.path(.dp_dir, sprintf("%s_cells_5km_16d.csv", prefix))
  if (!file.exists(fr)) stop("cell frame not found: ", fr)

  d <- fread(fr, select = c("cell_id","gx","gy","year","n_checklists","n_detections"),
             showProgress = FALSE)
  D <- fread(fdf, showProgress = FALSE)

  # JOIN ON (gx, gy). cell_id is not stable across separately built species and
  # the lookup carries whichever numbering it was generated under.
  map <- unique(d[, .(cell_id, gx, gy)])
  D <- merge(D[, .(gx, gy, fire_event_id, fire_year, dist_km)], map,
             by = c("gx","gy"))

  # A cell can only support a before/after comparison if it was birded both
  # before and after the fire. Without this the "fire count" is just a count of
  # fires near somewhere, which is not a design.
  obs <- d[, .(y_min = min(year), y_max = max(year)), by = cell_id]
  D <- merge(D, obs, by = "cell_id")
  D <- D[dist_km <= R & fire_year > y_min & fire_year <= y_max]

  if (!nrow(D)) stop("no qualifying (cell, fire) pairs at R = ", R, " km")

  # the EARLIEST qualifying fire defines the event, matching build_panel()'s
  # min(year) convention; ties broken by proximity
  setorder(D, cell_id, fire_year, dist_km)
  ev <- D[, .(ev_year = fire_year[1], dist_km = dist_km[1],
              fire = fire_event_id[1]), by = cell_id]

  treated <- ev$cell_id
  control <- setdiff(unique(d$cell_id), treated)

  ev[, dose := if (binary) 1
              else if (decay == "linear") pmax(0, 1 - dist_km / R)
              else exp(-dist_km / (R / 3))]

  # Controls need a pseudo event year or `post` absorbs the treatment. Same draw
  # and same seed position as build_panel(), so the two are comparable.
  set.seed(1)
  ctrl <- data.table(cell_id = control,
                     ev_year = sample(ev$ev_year, length(control), replace = TRUE),
                     dist_km = NA_real_, fire = NA_character_, dose = 0)

  panel <- rbind(d[cell_id %in% treated][, treat := 1L],
                 d[cell_id %in% control][, treat := 0L])
  panel <- merge(panel, rbind(ev, ctrl), by = "cell_id")
  panel[, post := as.integer(year >= ev_year)]

  cy <- panel[, .(det = sum(n_detections), chk = sum(n_checklists)),
              by = .(cell_id, year, treat, post, ev_year, dose, dist_km, fire)][chk > 0]
  cy[, `:=`(rate = det / chk, w = as.numeric(chk), post_d = as.numeric(post),
            post_dose = as.numeric(post) * as.numeric(dose))]
  # controls are singleton clusters, as everywhere else in this project
  cy[is.na(fire), fire := paste0("CTRL_", cell_id)]
  cy[, cl := as.integer(factor(fire))]

  setattr(cy, "dose_terms", "dose")
  setattr(cy, "radius_km", R)
  setattr(cy, "decay", decay)
  cy[]
}

#' Fit the distance panel and return one tidy row.
analyze_dist <- function(prefix = NULL, R = 10, decay = "linear", B = 999,
                         binary = FALSE, quiet = FALSE) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  p  <- build_dist_panel(prefix, R = R, decay = decay, binary = binary)
  xc <- c("post_d", "post_dose")
  P  <- .prep(p, xc)
  b  <- .fit(P, P$y)
  V  <- crve(P, as.vector(P$y - P$X %*% b))
  se <- sqrt(V[2, 2])
  base <- p[treat == 1 & post == 0, sum(det) / sum(w)]
  out <- data.table(
    species = prefix, radius_km = R, decay = decay,
    estimate = b[2], se = se, t = b[2] / se,
    p = wcb_p(P, 2L, b[2], se, B = B),
    pct_of_base = 100 * b[2] / base,
    n_cells = uniqueN(p[treat == 1]$cell_id),
    n_fires = uniqueN(p[treat == 1]$fire),
    n_ctrl  = uniqueN(p[treat == 0]$cell_id))
  if (!quiet) message(sprintf("  %-8s R=%2d  est %+.4f  t %+.2f  p %.3f  fires %d",
                              prefix, R, out$estimate, out$t, out$p, out$n_fires))
  out[]
}
