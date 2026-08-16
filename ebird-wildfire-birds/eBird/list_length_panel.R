# =============================================================================
# LIST LENGTH -- does normalising by recording activity remove the artifact?
#
# THE PROBLEM  (LIST_LENGTH_PRESPEC.md holds the pass/fail rule, fixed first)
#   Per-species occurrence effects are unusable: three negative controls
#   (Mallard, American Coot, Black Phoebe) show LARGER apparent fire effects
#   than the focal species. After a fire birders visit 43% fewer distinct
#   locations (p = 0.0001), which depresses recording for EVERY species.
#
# THE IDEA
#   If that shock is common across species, dividing a species out by total
#   recording activity cancels it -- the same cancellation that makes
#   guild_contrast.R work, applied per species instead of per guild.
#
# WHY LIST LENGTH EXCLUDING THE FOCAL SPECIES
#   Raw list length is a MEDIATOR, not just a nuisance. Fire depresses it via
#   observer behaviour (what we want to remove) but a genuine collapse in birds
#   would ALSO depress it (what we do not want to remove). Conditioning on it
#   raw would over-correct and eat part of the real effect.
#
#   So every measure below uses `ll_ex` = list length MINUS the focal species.
#   The denominator is then other species' detections: it still carries the
#   observer shock, but is not moved by the focal species' own abundance. This
#   is the same assumption the guild contrast already relies on -- that the
#   artifact is shared across species -- and no stronger.
#
# TWO FORMS, BOTH ON THE det_rate SCALE so they are comparable to the published
# point-grain control mean of 0.0188:
#   share      det / sum(ll_ex) * LLbar   -- rescaled by the global mean so that
#                                            if ll_ex were constant it equals
#                                            det_rate exactly
#   covariate  det_rate, with mean ll_ex entered as a regressor
#
# ARCHITECTURE
#   A PARALLEL builder. It does not modify point_panel.R or panel_utils.R, so
#   none of the 51 pinned numbers can move. It reuses the shared estimator
#   (.prep / .fit / crve / wcb_p) so the inference is identical.
#
# ENV  LL_R (radius km, default 10) · LL_MIN_YEARS (default 2) · LL_B (default 999)
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

.ll_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  for (f in c("project_config.R", "panel_utils.R", "infer.R")) {
    p <- file.path(.ll_dir, f)
    if (file.exists(p)) sys.source(p, envir = globalenv())
  }
})

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

#' Attach list length to a species' checklists.
#'
#' Returns the checklist table with `ll_ex` (list length excluding the focal
#' species) and hard-fails if the join is poor -- a silent NA rate here would
#' bias exactly the subset this is meant to fix, since group checklists are
#' longer and more productive than solo ones.
ll_checklists <- function(prefix) {
  jf <- file.path(.ll_dir, sprintf("%s_2015_2026_joined.csv", prefix))
  lf <- file.path(.ll_dir, "list_length_lookup.tsv")
  for (f in c(jf, lf)) if (!file.exists(f)) stop("missing input: ", basename(f))

  d <- fread(jf, select = c("checklist_id","latitude","longitude","year","presence"),
             showProgress = FALSE)
  L <- fread(lf, col.names = c("checklist_id","list_length"), showProgress = FALSE)
  # The lookup is built per source file and concatenated, and a handful of keys
  # (5 of 7,015,454 on the 2026-08-02 build -- four group ids, one event id)
  # appear in both date ranges. A keyed lookup returns BOTH rows and the
  # assignment then fails with a length mismatch, so collapse first. `max`
  # matches the group reconciliation rule in list_length_lookup.sh.
  nd <- nrow(L)
  L <- L[, .(list_length = max(list_length)), by = checklist_id]
  if (nrow(L) < nd) say(sprintf("  collapsed %d duplicate lookup key(s)", nd - nrow(L)))
  setkey(L, checklist_id)
  d[, list_length := L[.(d$checklist_id), list_length]]

  miss <- d[, mean(is.na(list_length))]
  say(sprintf("  list length matched %.2f%% of %s checklists",
              100 * (1 - miss), format(nrow(d), big.mark = ",")))
  if (miss > 0.05)
    stop(sprintf("list-length join failed on %.1f%% of checklists -- refusing to ",
                 100 * miss),
         "continue. Group ids ('G...') are 14.1% of checklist_id and must be ",
         "present in the lookup; check list_length_lookup.sh emitted both keys.")
  d <- d[!is.na(list_length)]

  # exclude the focal species from its own denominator -- see header
  d[, ll_ex := pmax(list_length - presence, 0L)]
  d[]
}

#' Spot x year panel with a chosen outcome. Mirrors build_point_panel(), which
#' it deliberately does not modify.
build_ll_panel <- function(prefix, outcome = c("det_rate","share","covariate"),
                           R = 10, MIN_YEARS = 2L) {
  outcome <- match.arg(outcome)
  mf <- file.path(.ll_dir, "point_grid_map.csv")
  df <- file.path(.ll_dir, "point_fire_distance.csv")
  for (f in c(mf, df)) if (!file.exists(f)) stop("missing input: ", basename(f))

  d <- ll_checklists(prefix)
  m <- fread(mf, showProgress = FALSE)
  d <- merge(d, m, by = c("latitude","longitude"))

  LLbar <- d[, mean(ll_ex)]                       # global scale anchor
  sy <- d[, .(det = sum(presence), chk = .N,
              sum_ll = sum(ll_ex), mean_ll = mean(ll_ex)), by = .(px, py, year)]
  keep <- sy[, .(ny = uniqueN(year)), by = .(px, py)][ny >= MIN_YEARS]
  sy <- merge(sy, keep[, .(px, py)], by = c("px","py"))

  # Both forms sit on the det_rate scale: with ll_ex constant at LLbar, `share`
  # reduces exactly to det/chk, so a coefficient here is directly comparable to
  # the published point-grain numbers.
  sy[, det_rate := det / chk]
  sy[, share    := fifelse(sum_ll > 0, det / sum_ll * LLbar, NA_real_)]
  sy <- sy[is.finite(share)]

  say(sprintf("  spot-years %s | mean ll_ex %.2f | outcome '%s'",
              format(nrow(sy), big.mark = ","), LLbar, outcome))

  D <- fread(df, showProgress = FALSE)[dist_km <= R]
  obs <- sy[, .(y_min = min(year), y_max = max(year)), by = .(px, py)]
  D <- merge(D[, .(px, py, fire_event_id, fire_year, dist_km)], obs, by = c("px","py"))
  D <- D[fire_year > y_min & fire_year <= y_max]
  if (!nrow(D)) stop("no qualifying spot-fire pairs at R = ", R)

  setorder(D, px, py, fire_year, dist_km)
  ev <- D[, .(ev_year = fire_year[1], dist_km = dist_km[1],
              fire = fire_event_id[1]), by = .(px, py)]
  ev[, dose := pmax(0, 1 - dist_km / R)]

  ctl <- fsetdiff(unique(sy[, .(px, py)]), unique(ev[, .(px, py)]))
  set.seed(1)                                     # same draw as point_panel.R
  ctl[, `:=`(ev_year = sample(ev$ev_year, .N, replace = TRUE),
             dist_km = NA_real_, fire = NA_character_, dose = 0)]

  p <- merge(sy, rbind(ev, ctl), by = c("px","py"))
  p[, `:=`(treat = as.integer(!is.na(fire)), post = as.integer(year >= ev_year))]
  p[, `:=`(rate = if (outcome == "share") share else det_rate,
           w = as.numeric(chk),
           post_d = as.numeric(post), post_dose = as.numeric(post) * dose,
           ll = mean_ll - LLbar)]
  p[, cell_id := .GRP, by = .(px, py)]
  p[is.na(fire), fire := paste0("CTRL_", cell_id)]
  p[, cl := as.integer(factor(fire))]
  setattr(p, "outcome", outcome)
  p[]
}

#' Fit. Identical estimator and inference to the rest of the project.
analyze_ll <- function(prefix, outcome = "det_rate", R = 10, MIN_YEARS = 2L,
                       B = 999) {
  p  <- build_ll_panel(prefix, outcome, R = R, MIN_YEARS = MIN_YEARS)
  xc <- c("post_d", "post_dose")
  if (outcome == "covariate") xc <- c(xc, "ll")   # list length as a regressor
  P  <- .prep(p, xc)
  b  <- .fit(P, P$y)
  se <- sqrt(crve(P, as.vector(P$y - P$X %*% b))[2, 2])
  data.table(species = prefix, outcome = outcome,
             estimate = b[2], se = se, t = b[2] / se,
             p = wcb_p(P, 2L, b[2], se, B = B),
             n_spots = uniqueN(p[treat == 1]$cell_id),
             n_fires = uniqueN(p[treat == 1]$fire))[]
}

# =============================================================================
if (sys.nframe() == 0L) {
  R  <- as.numeric(Sys.getenv("LL_R", "10"))
  MY <- as.integer(Sys.getenv("LL_MIN_YEARS", "2"))
  B  <- as.integer(Sys.getenv("LL_B", "999"))

  reg <- cfg_species(active_only = FALSE)
  # the three controls are the test; wrentit + olsfly ride along as the
  # published comparison, exactly as in point_controls_results.csv
  want <- c("mallar3", "y00475", "blkpho", "wrentit", "olsfly")
  want <- want[want %in% reg$prefix]

  res <- rbindlist(lapply(c("det_rate", "share", "covariate"), function(oc)
    rbindlist(lapply(want, function(pre) {
      say("=== ", pre, " / ", oc)
      tryCatch(analyze_ll(pre, oc, R = R, MIN_YEARS = MY, B = B),
               error = function(e) { say("  FAILED: ", conditionMessage(e)); NULL })
    }))), fill = TRUE)

  res <- merge(res, reg[, .(prefix, common_name, guild)],
               by.x = "species", by.y = "prefix", all.x = TRUE)
  fwrite(res, file.path(.ll_dir, "list_length_results.csv"))

  out <- file.path(.ll_dir, "list_length_summary.txt")
  sink(out)
  cat("DOES LIST-LENGTH NORMALISATION REMOVE THE OBSERVER ARTIFACT?\n")
  cat("============================================================\n")
  cat("Generated : ", format(Sys.time()), "\n", sep = "")
  cat("Verdict rule fixed in advance -- see LIST_LENGTH_PRESPEC.md\n")
  cat("Grain     : point (same birding spot before/after), radius ", R, " km\n", sep = "")
  cat("Inference : fire-clustered + wild cluster bootstrap, B = ", B, "\n\n", sep = "")
  cat("  det_rate  = detections / checklists          (the published outcome)\n")
  cat("  share     = detections / other-species detections, rescaled\n")
  cat("  covariate = det_rate with mean list length as a regressor\n\n")

  for (oc in c("det_rate", "share", "covariate")) {
    r <- res[outcome == oc]
    if (!nrow(r)) next
    cat("-- ", oc, " ", strrep("-", 60 - nchar(oc)), "\n", sep = "")
    print(as.data.frame(r[order(guild, species),
      .(common_name, guild, estimate = round(estimate, 4),
        se = round(se, 4), t = round(t, 2), p = round(p, 3), n_fires)]),
      row.names = FALSE)
    ct <- r[guild == "control"]
    if (nrow(ct)) {
      cat(sprintf("\n   CONTROLS: mean |effect| %.4f | %d of %d with p < 0.05\n\n",
                  mean(abs(ct$estimate)), sum(ct$p < 0.05), nrow(ct)))
    }
  }

  cat("\n-- VERDICT (rule fixed before the data existed) ----------------------\n")
  for (oc in c("share", "covariate")) {
    ct <- res[outcome == oc & guild == "control"]
    if (nrow(ct) < 3) { cat("  ", oc, ": incomplete, cannot judge\n", sep = ""); next }
    ok <- all(ct$p >= 0.05) && mean(abs(ct$estimate)) < 0.0100
    cat(sprintf("  %-10s %s  (all controls p>=0.05: %s | mean |effect| %.4f < 0.0100: %s)\n",
                oc, if (ok) "CLEAN" else "NOT CLEAN",
                all(ct$p >= 0.05), mean(abs(ct$estimate)),
                mean(abs(ct$estimate)) < 0.0100))
  }
  cat("\n  Published baseline for comparison: point-grain control mean |effect|\n")
  cat("  = 0.0188, with all three controls significant.\n")
  cat("\n  A CLEAN result restores per-species occurrence work. It does NOT\n")
  cat("  vindicate the guild contrast: this cancels effects COMMON to all\n")
  cat("  species, and differential detectability -- burning off dense chaparral\n")
  cat("  makes a skulking bird easier to see -- would survive it.\n")
  sink()
  cat(readLines(out), sep = "\n")
  say("DONE -> list_length_results.csv / list_length_summary.txt")
}
