# =============================================================================
# Shared panel construction, swappable burn-impact measure, shared estimator
#
# WHY THIS EXISTS
#   TWO problems, one file.
#
#   1. DUPLICATION. The ~40 lines that define treated/control, assign pseudo
#      event years to controls, attach the dose and aggregate to cell-year were
#      copied into SIX scripts (dose_burn_effect, fe_decompose, fire_cluster_se,
#      cs_estimator, wild_cluster_boot, effort_placebo). Change the panel
#      definition and six files must change together or they silently diverge --
#      and divergent panels produce results that disagree for reasons nobody can
#      trace later.
#
#   2. THE DOSE WAS HARDCODED. Every causal script used burn_frac_max. But the
#      research question says SEVERITY, and severity lives in a different place
#      (the regrid cache's per-class fractions) from extent (the frame). Testing
#      a different burn measure meant editing each script.
#
#   build_panel() fixes both: one definition, and the burn measure is a name
#   from treatment_registry.csv.
#
# REPRODUCIBILITY
#   set.seed(1) is applied immediately before the control pseudo-event-year
#   draw, exactly where dose_burn_effect.R applies it, so a panel built here is
#   IDENTICAL to the one the published numbers came from. Verified by
#   panel_selftest() below -- run it after any change to this file.
#
# USAGE
#   source(file.path(dir_, "panel_utils.R"))
#   p  <- build_panel("wrentit", dose = "frac_max")
#   tm <- attr(p, "dose_terms")          # e.g. "dose", or c("extent","high_sev")
#   b  <- twfe(p, c("post_d", paste0("post_", tm)))
# =============================================================================
suppressMessages(library(data.table))

.pu_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))

# session cache for cell-frame reads; see build_panel(). clear with pu_clear_cache()
.pu_cache <- new.env(parent = emptyenv())
pu_clear_cache <- function() rm(list = ls(.pu_cache), envir = .pu_cache)

# cfg_species(), cfg_catalog(), cfg_add_species() etc. live in project_config.R.
# Source it here so ONE source("panel_utils.R") gives the whole config surface --
# otherwise cfg_treatments() works and cfg_species() mysteriously does not.
local({
  f <- file.path(.pu_dir, "project_config.R")
  if (file.exists(f) && !exists("cfg_species", mode = "function")) sys.source(f, envir = globalenv())
})

# --- treatment registry -------------------------------------------------------

cfg_treatments <- function(active_only = FALSE) {
  f <- file.path(.pu_dir, "treatment_registry.csv")
  if (!file.exists(f)) stop("treatment_registry.csv not found in ", .pu_dir)
  t <- fread(f, showProgress = FALSE)
  if (active_only) t <- t[active == 1]
  t[]
}

#' Parse a registry `terms` string into named expressions.
#' "extent=frac_tot;high_sev=s4" -> c(extent = "frac_tot", high_sev = "s4")
#'
#' Splits on the FIRST "=" only. Splitting on every "=" silently truncates any
#' expression containing >=, <= or == -- e.g. "dose=as.numeric(x>=2)" became
#' "as.numeric(x>" and failed to parse.
.parse_terms <- function(s) {
  chunks <- trimws(strsplit(s, ";", fixed = TRUE)[[1]])
  nm <- sub("=.*$", "", chunks)
  ex <- sub("^[^=]*=", "", chunks)
  setNames(trimws(ex), trimws(nm))
}

#' Column names a set of term expressions actually references.
#' Derived with all.vars() rather than a hardcoded list, so ANY column in the
#' frame can be used in a registry expression without editing this file.
.terms_vars <- function(terms) {
  unique(unlist(lapply(terms, function(e) all.vars(parse(text = e)))))
}

#' Per-cell FIRE ATTRIBUTES from the MTBS perimeter shapefile, joined through
#' cell_fire_lookup.R's dominant-fire mapping.
#'
#' These are properties of the FIRE, not of the cell's burn: how big it was, what
#' month it started, whether it was a wildfire or a prescribed burn, and the
#' per-fire dNBR calibration. Nothing in the project uses them yet, and they are
#' the only burn information in MTBS that the cell frame does not carry.
#'
#' Returns one row per cell: fire_ac, fire_month, fire_doy, is_wildfire,
#' dnbr_offst, dnbr_stddv, low_t, mod_t, high_t.
#'
#' KEYED TO THE REQUESTED SPECIES' cell_id, via (gx, gy).
#' cell_fire_lookup.csv is written from whichever frame cell_fire_lookup.R last
#' ran on (wrentit by default), so its cell_id column belongs to THAT numbering.
#' This function used to merge it onto the panel on cell_id directly -- the exact
#' unsafe join attach_fire() was fixed to avoid, twelve lines below. Measured:
#' the lookup matches wrentit on 3,125 of 3,125 cells and calthr on 0 of 3,125,
#' and build_panel("calthr","fire_size") returned silently with 760 of 1,009
#' treated cells at dose 0 and the other 249 attached to unrelated fires.
#' (gx, gy) is the geography and is stable across separately built species.
.perim_attrs <- function(prefix = NULL) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  ckey <- paste0("__perim|", prefix)        # cache per numbering, not globally
  if (!is.null(.pu_cache[[ckey]])) return(copy(.pu_cache[[ckey]]))
  suppressMessages(requireNamespace("sf", quietly = TRUE)) ||
    stop("sf needed for perimeter-sourced burn measures")
  lkf <- file.path(.pu_dir, "cell_fire_lookup.csv")
  if (!file.exists(lkf)) stop("cell_fire_lookup.csv not found -- run cell_fire_lookup.R first")
  lk <- fread(lkf, showProgress = FALSE)
  if (!all(c("gx","gy") %in% names(lk)))
    stop("cell_fire_lookup.csv has no gx/gy columns. Re-run cell_fire_lookup.R -- ",
         "joining it on cell_id alone is unsafe across species.")
  shp <- file.path(.pu_dir, "mtbs_perim", "mtbs_perims_DD.shp")
  p <- sf::st_read(shp, quiet = TRUE, query =
    "SELECT Event_ID, Ig_Date, BurnBndAc, Incid_Type, dNBR_offst, dNBR_stddv, Low_T, Mod_T, High_T
       FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
  a <- as.data.table(sf::st_drop_geometry(p))
  setnames(a, tolower(names(a)))
  a[, ig_date := as.Date(ig_date)]
  a[, `:=`(fire_ac = as.numeric(burnbndac),
           fire_month = as.integer(format(ig_date, "%m")),
           fire_doy = as.integer(format(ig_date, "%j")),
           is_wildfire = as.integer(grepl("Wildfire", incid_type, ignore.case = TRUE)))]
  # a cell may appear for several fire years; keep its most recent fire, matching
  # how fire_cluster_se.R and cs_estimator.R pick the causing fire. Grouped by
  # (gx, gy) rather than cell_id -- they are 1:1 within the lookup, and this keeps
  # the dominant-fire choice independent of any numbering.
  setorder(lk, gx, gy, -fire_year)
  dom <- lk[, .(fire_event_id = fire_event_id[1]), by = .(gx, gy)]

  # translate the lookup's geography into THIS species' cell_id
  fr <- file.path(.pu_dir, sprintf("%s_cells_5km_16d.csv", prefix))
  if (!file.exists(fr)) stop("cell frame not found: ", fr)
  map <- unique(fread(fr, select = c("cell_id","gx","gy"), showProgress = FALSE))
  dom <- merge(dom, map, by = c("gx","gy"))
  if (nrow(dom) < 0.75 * uniqueN(lk[, .(gx, gy)]))
    warning(sprintf("perim_attrs: only %d of %d lookup cells matched '%s' on (gx,gy). ",
                    nrow(dom), uniqueN(lk[, .(gx, gy)]), prefix),
            "Expect ~all. The frame and the lookup may cover different grids.")

  out <- merge(dom[, .(cell_id, fire_event_id)], a,
               by.x = "fire_event_id", by.y = "event_id", all.x = TRUE)
  out <- out[, .(cell_id, fire_ac, fire_month, fire_doy, is_wildfire,
                 dnbr_offst = as.numeric(dnbr_offst), dnbr_stddv = as.numeric(dnbr_stddv),
                 low_t = as.numeric(low_t), mod_t = as.numeric(mod_t),
                 high_t = as.numeric(high_t))]
  .pu_cache[[ckey]] <- out
  copy(out)
}

#' Per-cell severity composition from the regrid cache, restricted to in-window
#' fires. Same construction as severity_effect.R so the two agree.
.cache_sev <- function(prefix = NULL) {
  # THE CACHE IS KEYED BY cell_id, AND cell_id IS NOT STABLE ACROSS SPECIES.
  # regrid_cells.R now writes burn_<tag>.rds, one per numbering, so a new
  # species cannot silently re-key another's cache. Before that it wrote a single
  # burn.rds: running the calthr regrid on 2026-07-31 left it keyed to calthr
  # while wrentit/brncre/olsfly still expected their own, and every cache-sourced
  # measure would have mis-joined with no error. Pick the file matching THIS
  # species' numbering, and fail loudly if none does.
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  cdir <- file.path(.pu_dir, ".regrid_cache")
  fr   <- file.path(.pu_dir, sprintf("%s_cells_5km_16d.csv", prefix))
  # The frame is REQUIRED, because it is the only thing that says which numbering
  # this species uses. The old code left the fingerprint NULL when the frame was
  # missing and then accepted the FIRST cache it found -- so a typo'd or unbuilt
  # prefix silently received another species' severity composition, keyed to a
  # numbering that has nothing to do with it. That is precisely the failure this
  # function exists to prevent, so it must be an error, not a fallback.
  if (!file.exists(fr))
    stop("no cell frame for '", prefix, "': ", basename(fr), " not found.\n",
         "  A cache-sourced burn measure needs the frame to identify the cell\n",
         "  numbering; without it any cache would be a guess. Build the species\n",
         "  first (./run_species.sh ", prefix, "), or check the prefix spelling.")
  d  <- unique(fread(fr, select = c("cell_id","gx","gy"), showProgress = FALSE))
  fp <- c(5000, nrow(d), sum(as.numeric(d$gx) * d$cell_id),
                         sum(as.numeric(d$gy) * d$cell_id))
  cands <- c(list.files(cdir, "^burn_[0-9a-f]+\\.rds$", full.names = TRUE),
             file.path(cdir, "burn.rds"))
  cands <- cands[file.exists(cands)]
  if (!length(cands)) stop("no regrid burn cache in ", cdir,
    " -- rerun regrid_cells.R, or use a frame-sourced burn measure.")
  bl <- NULL
  for (cf in cands) {
    b <- readRDS(cf)
    # fp is never NULL now -- a missing frame errors above rather than matching
    # the first cache indiscriminately
    if (isTRUE(all.equal(as.numeric(b$grid), fp))) { bl <- b; break }
  }
  if (is.null(bl))
    stop("no regrid cache matches the cell numbering of '", prefix, "'.\n",
         "  Checked: ", paste(basename(cands), collapse = ", "), "\n",
         "  cell_id is not comparable across separately built species.\n",
         "  Fix: EBIRD_PREFIX=", prefix, " Rscript regrid_cells.R  (~6 min),\n",
         "  or use a frame-sourced measure (frac_max, frac_cum, sev_weighted, ...).")
  cb <- as.data.table(bl$cellburn)
  # classes 5 (increased greenness) and 6 (mask) are excluded from the dose:
  # they are not burn severity, and together are 0.4% of burned area
  w <- dcast(cb[fire_year >= 2015 & sev %in% 1:4, .(frac = sum(frac)), by = .(cell_id, sev)],
             cell_id ~ sev, value.var = "frac", fill = 0)
  setnames(w, c("cell_id", paste0("s", 1:4)))
  w[, frac_tot := s1 + s2 + s3 + s4]
  w[, mean_sev := (1*s1 + 2*s2 + 3*s3 + 4*s4) / pmax(frac_tot, 1e-9)]
  w[]
}

# --- the panel ----------------------------------------------------------------

#' Per-CELL treatment assignment: who is treated, when their fire was, and their
#' dose. This is the piece every causal script shares, independent of the time
#' granularity it then aggregates to.
#'
#' build_panel() aggregates to cell-YEAR, which suits most analyses. But
#' migration_phenology.R works at 16-day resolution and event_study.R needs event
#' time, so they need the assignment WITHOUT the aggregation. Exposing it here
#' keeps one definition of "treated" instead of five.
#'
#' @return data.table: cell_id, treat, ev_year, + one column per dose term.
#'         attr "dose_terms" names the dose columns.
cell_treatment <- function(prefix = NULL, dose = "frac_max") {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  tr <- cfg_treatments()[name == dose]
  if (!nrow(tr)) stop("dose '", dose, "' is not in treatment_registry.csv. Available: ",
                      paste(cfg_treatments()$name, collapse = ", "))
  terms <- .parse_terms(tr$terms[1])

  fp <- file.path(.pu_dir, sprintf("%s_cells_5km_16d.csv", prefix))
  if (!file.exists(fp)) stop("cell frame not found: ", fp)
  have <- strsplit(readLines(fp, n = 1), ",")[[1]]
  need <- unique(c("cell_id","year","ever_burned", .terms_vars(terms)))
  ckey <- paste("ct", prefix, paste(sort(intersect(need, have)), collapse = ","), sep = "|")
  if (is.null(.pu_cache[[ckey]]))
    .pu_cache[[ckey]] <- fread(fp, select = intersect(need, have), showProgress = FALSE)
  d <- copy(.pu_cache[[ckey]])

  st <- d[, .(pre = sum(ever_burned == 0), post = sum(ever_burned == 1)), by = cell_id]
  treated <- st[pre > 0 & post > 0]$cell_id
  control <- st[post == 0]$cell_id
  ev <- d[cell_id %in% treated & ever_burned == 1, .(ev_year = min(year)), by = cell_id]

  # controls need a pseudo event year or `post` absorbs the treatment; drawn from
  # the treated distribution. SEED HERE, matching dose_burn_effect.R.
  set.seed(1)
  ctrl_ev <- data.table(cell_id = control,
                        ev_year = sample(ev$ev_year, length(control), replace = TRUE))

  if (tr$source[1] == "cache") {
    src <- .cache_sev(prefix)[cell_id %in% treated]
  } else if (tr$source[1] == "perimeter") {
    src <- .perim_attrs(prefix)[cell_id %in% treated]
  } else {
    src <- d[cell_id %in% treated & ever_burned == 1]
  }
  dose_dt <- if (tr$source[1] == "frame") {
    src[, lapply(terms, function(e) eval(parse(text = e), envir = .SD)), by = cell_id]
  } else {
    x <- src[, lapply(terms, function(e) eval(parse(text = e), envir = .SD))]
    x[, cell_id := src$cell_id]; x
  }
  setcolorder(dose_dt, "cell_id")

  out <- rbind(data.table(cell_id = treated, treat = 1L),
               data.table(cell_id = control, treat = 0L))
  out <- merge(out, rbind(ev, ctrl_ev), by = "cell_id")
  out <- merge(out, dose_dt, by = "cell_id", all.x = TRUE)
  for (tn in names(terms)) set(out, which(!is.finite(out[[tn]])), tn, 0)
  setattr(out, "dose_terms", names(terms))
  setattr(out, "dose_label", tr$label[1])
  out[]
}

#' Build the cell-year DiD panel.
#'
#' @param prefix   species prefix (e.g. "wrentit"); defaults to EBIRD_PREFIX
#' @param dose     name from treatment_registry.csv
#' @param extra    extra frame columns to carry through, aggregated by mean
#' @return data.table with cell_id, year, treat, post, det, chk, rate, w,
#'         one column per dose term, post_d, and post_<term> interactions.
#'         attr(x, "dose_terms") lists the term names.
build_panel <- function(prefix = NULL, dose = "frac_max", extra = NULL,
                        extra_sum = NULL) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  tr <- cfg_treatments()[name == dose]
  if (!nrow(tr)) stop("dose '", dose, "' is not in treatment_registry.csv. Available: ",
                      paste(cfg_treatments()$name, collapse = ", "))
  terms <- .parse_terms(tr$terms[1])

  # Columns the registry expression references are derived, not hardcoded -- a
  # fixed list silently broke every measure using a column that was not on it
  # (burn_decay, n_fire_epochs, the six time-band columns...).
  need <- unique(c("cell_id","year","ever_burned","n_checklists","n_detections",
                   "burn_frac_max", .terms_vars(terms), extra, extra_sum))
  fp <- file.path(.pu_dir, sprintf("%s_cells_5km_16d.csv", prefix))
  if (!file.exists(fp)) stop("cell frame not found: ", fp)
  have <- strsplit(readLines(fp, n = 1), ",")[[1]]

  # Frame reads dominate a sweep: the CSV is ~200 MB and a 3-species x 10-dose
  # grid would read it 30 times. Cache per (prefix, column set) for the session.
  ckey <- paste(prefix, paste(sort(intersect(need, have)), collapse = ","), sep = "|")
  if (is.null(.pu_cache[[ckey]])) {
    .pu_cache[[ckey]] <- fread(fp, select = intersect(need, have), showProgress = FALSE)
  }
  d <- copy(.pu_cache[[ckey]])

  # treated = observed unburned at some point AND burned at some point, in window
  st <- d[, .(pre = sum(ever_burned == 0), post = sum(ever_burned == 1)), by = cell_id]
  treated <- st[pre > 0 & post > 0]$cell_id
  control <- st[post == 0]$cell_id
  ev <- d[cell_id %in% treated & ever_burned == 1, .(ev_year = min(year)), by = cell_id]

  # Controls need a pseudo "post" or the post term absorbs the treatment. Drawn
  # from the treated event-year distribution. SEED HERE, matching
  # dose_burn_effect.R, so the panel reproduces the published numbers exactly.
  set.seed(1)
  ctrl_ev <- data.table(cell_id = control,
                        ev_year = sample(ev$ev_year, length(control), replace = TRUE))

  # --- the dose, from whichever source the registry names --------------------
  # cache source: already one row per cell, so the expressions are vectorised.
  # frame source: many rows per cell, so the expressions AGGREGATE (max(), sum()).
  if (tr$source[1] == "cache") {
    src <- .cache_sev(prefix)[cell_id %in% treated]
    dose_dt <- src[, lapply(terms, function(e) eval(parse(text = e), envir = .SD))]
    dose_dt[, cell_id := src$cell_id]
  } else if (tr$source[1] == "perimeter") {
    src <- .perim_attrs(prefix)[cell_id %in% treated]
    dose_dt <- src[, lapply(terms, function(e) eval(parse(text = e), envir = .SD))]
    dose_dt[, cell_id := src$cell_id]
    for (tn in names(terms))                      # unmatched fires -> zero dose
      set(dose_dt, which(!is.finite(dose_dt[[tn]])), tn, 0)
  } else {
    src <- d[cell_id %in% treated & ever_burned == 1]
    dose_dt <- src[, lapply(terms, function(e) eval(parse(text = e), envir = .SD)),
                   by = cell_id]
  }
  setcolorder(dose_dt, "cell_id")
  stopifnot(ncol(dose_dt) == length(terms) + 1L, !anyDuplicated(dose_dt$cell_id))

  panel <- rbind(d[cell_id %in% treated][, treat := 1L],
                 d[cell_id %in% control][, treat := 0L])
  panel <- merge(panel, rbind(ev, ctrl_ev), by = "cell_id")
  panel <- merge(panel, dose_dt, by = "cell_id", all.x = TRUE)
  for (tn in names(terms)) set(panel, which(is.na(panel[[tn]])), tn, 0)  # controls: zero dose
  panel[, post := as.integer(year >= ev_year)]

  # `extra` columns are averaged over the 16-day steps in the year; `extra_sum`
  # columns are totalled. effort_hours must be SUMMED (it is a total, and
  # effort_placebo.R's published numbers use the sum) while the per-checklist
  # medians must be AVERAGED.
  agg <- c(list(det = quote(sum(n_detections)), chk = quote(sum(n_checklists))),
           setNames(lapply(extra, function(x) bquote(mean(.(as.name(x)), na.rm = TRUE))), extra),
           setNames(lapply(extra_sum, function(x) bquote(sum(.(as.name(x)), na.rm = TRUE))), extra_sum))
  # ev_year is cell-constant so adding it to `by` changes no groups; it is
  # carried because build_cs_panel() needs the cohort and would otherwise have
  # to recompute the treated/control split from scratch.
  cy <- panel[, eval(as.call(c(quote(list), agg))),
              by = c("cell_id","year","treat","post","ev_year", names(terms))][chk > 0]
  cy[, `:=`(rate = det / chk, w = as.numeric(chk), post_d = as.numeric(post))]
  for (tn in names(terms)) set(cy, NULL, paste0("post_", tn),
                               as.numeric(cy$post) * as.numeric(cy[[tn]]))
  setattr(cy, "dose_terms", names(terms))
  setattr(cy, "dose_name", dose)
  setattr(cy, "dose_label", tr$label[1])
  cy[]
}

#' Panel in the shape did::att_gt() wants: cell-year, plus `first_treat` = the
#' period a cell is first treated (0 = never treated).
#'
#' Differs from build_panel() in three ways that att_gt() forces:
#'   * BINARY treatment. att_gt takes staggered on/off, not a continuous dose,
#'     so treatment is "burned >= thresh" and the dose only picks the subset.
#'     This is a DIFFERENT ESTIMAND from the continuous-dose headline.
#'   * Controls are genuinely never-treated -- no pseudo event years needed.
#'   * Cohorts below `min_cohort` are dropped. The 2025 cohort has 9 cells at
#'     >= 50% and almost no post-period, and returned a degenerate ATT of
#'     -0.5001 with SE 0.0016 before this filter existed.
#' Cells with a single observation are dropped (att_gt requires >= 2).
build_cs_panel <- function(prefix = NULL, dose = "frac_max", thresh = 0.50,
                           min_cohort = 20L) {
  cy <- build_panel(prefix, dose = dose)
  tm <- attr(cy, "dose_terms"); stopifnot(length(tm) == 1L)
  trt <- unique(cy[treat == 1 & get(tm[1]) >= thresh, .(cell_id, ev_year)])
  csz <- trt[, .N, by = ev_year][order(ev_year)]
  drop_g <- csz[N < min_cohort]$ev_year
  if (length(drop_g)) {
    message("  dropping cohorts with < ", min_cohort, " cells: ",
            paste(sprintf("%d (n=%d)", csz[N < min_cohort]$ev_year,
                          csz[N < min_cohort]$N), collapse = ", "))
    trt <- trt[!ev_year %in% drop_g]
  }
  keep <- cy[cell_id %in% trt$cell_id | treat == 0,
             .(cell_id, year, det, chk, rate, w, treat)]
  keep <- merge(keep, trt[, .(cell_id, first_treat = as.integer(ev_year))],
                by = "cell_id", all.x = TRUE)
  keep[is.na(first_treat), first_treat := 0L]
  nobs <- keep[, .N, by = cell_id][N >= 2]
  keep <- keep[cell_id %in% nobs$cell_id]
  setattr(keep, "thresh", thresh); setattr(keep, "dose_name", dose)
  keep[]
}

# --- shared estimator ---------------------------------------------------------

#' Weighted two-way FE by iterative demeaning.
#'
#' @param yname   outcome column (effort_placebo.R swaps this for effort measures)
#' @param mode    which fixed effects to absorb: "both" (cell+year), "cell",
#'                "year", or "none". fe_decompose.R uses this to find which FE
#'                flips a sign. With "none" the caller must supply an intercept
#'                column; with "year" the intercept is absorbed.
#' @param weights FALSE gives every cell-year equal weight. Note: weighting by
#'                checklist count is circular when the OUTCOME is checklist count.
#'
#' Everything is cast to double first: data.table := into an INTEGER column
#' truncates the demeaned values to 0 and silently returns NaN.
#'
#' One FE dimension is exact in a single pass; two need iteration.
twfe <- function(dt, xcols, iters = 40, weights = TRUE, yname = "rate",
                 mode = c("both","cell","year","none")) {
  mode <- match.arg(mode)
  z <- copy(dt)[, c("cell_id","year",yname,"w", xcols), with = FALSE]
  for (cc in c(yname, xcols)) set(z, NULL, cc, as.numeric(z[[cc]]))
  if (!weights) set(z, NULL, "w", 1)
  # NB: no is.finite() filter here. It would run on every call, which inside a
  # 200-rep bootstrap means hundreds of extra scans of a 62k-row table for no
  # benefit -- `rate` is always finite (det/chk with chk > 0). Callers using a
  # different outcome that can be NA (effort_placebo.R) pre-filter instead.
  if (mode != "none") {
    n_it <- if (mode == "both") iters else 1L
    for (i in seq_len(n_it)) for (cc in c(yname, xcols)) {
      if (mode %in% c("cell","both"))
        z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
      if (mode %in% c("year","both"))
        z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
    }
  }
  X <- as.matrix(z[, ..xcols]); y <- z[[yname]]; wt <- z$w
  XtWX <- t(X * wt) %*% X
  if (abs(det(XtWX)) < 1e-18) return(setNames(rep(NA_real_, length(xcols)), xcols))
  as.vector(solve(XtWX, t(X * wt) %*% y)) |> setNames(xcols)
}

#' Attach the dominant fire per cell as a cluster label (from cell_fire_lookup.R).
#' Controls become singleton clusters.
#' CRITICAL: joins on (gx, gy), NOT cell_id.
#'
#' cell_id is NOT stable across separately-built species. regrid_cells.R
#' assigns it with `cell_id := .I` over unique (gx,gy) in ROW-APPEARANCE order,
#' which depends on the order of the source joined table. wrentit, brncre and
#' olsfly happen to share a numbering because they were built from the same
#' combined zero-fill; a species built from a per-species zero-fill gets a
#' DIFFERENT numbering for the same geography (verified: calthr matches wrentit
#' on 0 of 13,026 cells).
#'
#' Joining the fire lookup on cell_id therefore silently mis-assigns fires for
#' any newly built species -- it quietly turned 827 of 1,009 treated cells into
#' singleton clusters and inflated the fire count from 225 to 936. (gx, gy) is
#' the geography and is stable by construction.
attach_fire <- function(cy, prefix = NULL) {
  if (is.null(prefix)) prefix <- cfg_default_prefix()
  lk <- fread(file.path(.pu_dir, "cell_fire_lookup.csv"), showProgress = FALSE)
  if (!all(c("gx","gy") %in% names(lk)))
    stop("cell_fire_lookup.csv has no gx/gy columns. Re-run cell_fire_lookup.R -- ",
         "joining it on cell_id alone is unsafe across species.")
  d <- fread(file.path(.pu_dir, sprintf("%s_cells_5km_16d.csv", prefix)),
             select = c("cell_id","gx","gy","year","ever_burned"), showProgress = FALSE)
  map <- unique(d[, .(cell_id, gx, gy)])
  ev  <- d[ever_burned == 1, .(ev_year = min(year)), by = cell_id]
  ev  <- merge(ev, map, by = "cell_id")
  evt <- merge(ev, lk[, .(gx, gy, fire_year, fire_event_id)],
               by = c("gx","gy"), allow.cartesian = TRUE)[fire_year <= ev_year]
  setorder(evt, cell_id, -fire_year)
  cf <- evt[, .(fire = fire_event_id[1]), by = cell_id]

  out <- merge(cy, cf, by = "cell_id", all.x = TRUE)
  n_unmatched <- uniqueN(out[treat == 1 & is.na(fire)]$cell_id)
  n_treated   <- uniqueN(out[treat == 1]$cell_id)
  if (n_treated > 0 && n_unmatched > 0.25 * n_treated)
    warning(sprintf("attach_fire: %d of %d treated cells (%0.0f%%) matched no fire. ",
                    n_unmatched, n_treated, 100*n_unmatched/n_treated),
            "Expect ~0. This usually means the frame and the lookup disagree.")
  out[is.na(fire), fire := paste0("CTRL_", cell_id)]
  out[, cl := as.integer(factor(fire))]
  setattr(out, "dose_terms", attr(cy, "dose_terms"))
  out[]
}

# --- self test ----------------------------------------------------------------

#' Confirm build_panel("wrentit","frac_max") still reproduces the published
#' continuous dose of -0.0298. Run after ANY change to this file.
panel_selftest <- function(tol = 5e-4) {
  p <- build_panel("wrentit", dose = "frac_max")
  b <- twfe(p, c("post_d","post_dose"))[["post_dose"]]
  ok <- abs(b - (-0.0298)) < tol
  message(sprintf("selftest: continuous dose = %+.4f (published -0.0298) -> %s",
                  b, if (ok) "MATCH" else "*** MISMATCH ***"))
  invisible(ok)
}
