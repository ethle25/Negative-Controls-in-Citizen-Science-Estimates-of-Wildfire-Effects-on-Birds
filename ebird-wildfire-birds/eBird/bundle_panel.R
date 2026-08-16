# =============================================================================
# STEP 3 -- BUNDLES of birds: species richness and guild aggregates.
#
# WHY
#   Francisco's first To Do item was "different measures of biodiversity ...
#   bundles of bird species". Every result so far is one species at a time, which
#   asks "does fire move THIS bird". A bundle asks two different questions:
#
#     GUILD RATE  do shrub-living birds AS A GROUP decline where fire burns?
#                 Pools information across species, so it can be sharper than any
#                 single member -- the shrub guild has 6 species, the snag guild 6.
#     RICHNESS    does fire change HOW MANY KINDS of bird are seen, not just
#                 which? This is a genuinely different outcome and the project
#                 cannot currently ask it at all.
#
# WHY IT NEEDS NO NEW COLUMNS
#   Adding a column to the frame is the last hardcoded surface in this project
#   (it means hand-editing regrid_cells.R). This sidesteps it: the per-species
#   frames already hold n_detections, so a bundle is a JOIN across frames, not a
#   new build. Nothing in the pipeline changes.
#
# THE JOIN KEY IS (gx, gy, year)
#   NOT cell_id. Species built from separate zero-fills have different cell
#   numberings for identical geography -- calthr matches wrentit on 0 of 13,026
#   cells. This trap has bitten the project three times.
#
# THE CONTROLS ARE THE POINT
#   The `control` guild -- Mallard, American Coot, Black Phoebe -- are wetland and
#   open-water birds. Their bundle should show NO fire effect. If it does, the
#   design is broken, and that is worth more than another confirming guild.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

.bp_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  f <- file.path(.bp_dir, "panel_utils.R")
  if (file.exists(f) && !exists("twfe", mode = "function")) sys.source(f, envir = globalenv())
})

#' Detections per (gx, gy, year) for every built species, in one table.
#' Cached in the session -- reading 19 frames is the expensive part.
.bundle_raw <- function() {
  if (!is.null(.pu_cache[["__bundle"]])) return(copy(.pu_cache[["__bundle"]]))
  s <- cfg_species(active_only = FALSE)
  s <- s[file.exists(file.path(.bp_dir, sprintf("%s_cells_5km_16d.csv", prefix)))]
  message("  reading ", nrow(s), " species frames")
  L <- lapply(seq_len(nrow(s)), function(i) {
    f <- fread(file.path(.bp_dir, sprintf("%s_cells_5km_16d.csv", s$prefix[i])),
               select = c("gx","gy","year","n_checklists","n_detections"),
               showProgress = FALSE)
    f <- f[, .(det = sum(n_detections), chk = sum(n_checklists)), by = .(gx, gy, year)]
    f[, `:=`(prefix = s$prefix[i], guild = s$guild[i])]
    f
  })
  out <- rbindlist(L)
  .pu_cache[["__bundle"]] <- out
  copy(out)
}

#' Build a cell-year panel whose OUTCOME is a bundle.
#'
#' @param guild  which guild to bundle ("shrub","forest","snag","control"), or
#'               "all" for every species
#' @param what   "rate"     = detections summed over the guild / checklists
#'               "richness" = how many of the guild's species were detected
#' @param dose   burn measure from treatment_registry.csv
build_bundle_panel <- function(guild = "shrub", what = c("rate","richness"),
                               dose = "frac_max") {
  what <- match.arg(what)
  b <- .bundle_raw()
  # NEVER name an argument after a data.table column. Inside b[...], `guild`
  # resolves to the COLUMN, so `b[guild == guild]` compares it with itself and
  # matches everything -- and `..guild` does not reach a function argument, which
  # is what produced "object 'guild' not found". CLAUDE.md logs this exact trap
  # biting twice before (`prefix`, `common_name`); this is the third time.
  g_ <- guild
  if (g_ != "all") b <- b[guild == g_]
  if (!nrow(b)) stop("no built species in guild '", g_, "'")
  nsp <- uniqueN(b$prefix)

  # checklists are a property of the CELL-TIME, identical across species, so take
  # one copy rather than summing it nsp times
  agg <- b[, .(det = sum(det),
               n_det_species = sum(det > 0),
               chk = chk[1]), by = .(gx, gy, year)]
  agg[, y := if (what == "rate") det / chk else n_det_species]

  # treatment comes from the SHARED definition so this is comparable to every
  # single-species result; the reference species only supplies the geography
  ref <- cfg_species()$prefix[1]
  ct  <- cell_treatment(ref, dose = dose)
  tm  <- attr(ct, "dose_terms")
  if (length(tm) > 1L) stop("bundle panels take a single-term dose; got ", length(tm))
  map <- unique(fread(file.path(.bp_dir, sprintf("%s_cells_5km_16d.csv", ref)),
                      select = c("cell_id","gx","gy"), showProgress = FALSE))
  ct  <- merge(ct, map, by = "cell_id")             # -> (gx, gy), the safe key

  p <- merge(agg, ct[, .(gx, gy, treat, ev_year, dose = get(tm[1]))], by = c("gx","gy"))
  p[, post := as.integer(year >= ev_year)]
  p[, `:=`(rate = y, w = as.numeric(chk),
           post_d = as.numeric(post), post_dose = as.numeric(post) * dose)]
  p[, cell_id := .GRP, by = .(gx, gy)]

  setattr(p, "dose_terms", "dose")
  setattr(p, "guild", guild); setattr(p, "what", what); setattr(p, "n_species", nsp)
  p[]
}

#' Fit a bundle and return one tidy row, clustered by fire.
analyze_bundle <- function(guild = "shrub", what = "rate", dose = "frac_max",
                           B = 999) {
  p   <- build_bundle_panel(guild, what, dose)
  ref <- cfg_species()$prefix[1]
  p   <- attach_fire(p, ref)
  xc  <- c("post_d","post_dose")
  P   <- .prep(p, xc)
  bh  <- .fit(P, P$y)
  se  <- sqrt(crve(P, as.vector(P$y - P$X %*% bh))[2,2])
  base <- p[treat == 1 & post == 0, weighted.mean(rate, w)]
  data.table(guild = guild, outcome = what, n_species = attr(p, "n_species"),
             estimate = bh[2], se = se, t = bh[2]/se,
             p = wcb_p(P, 2L, bh[2], se, B = B),
             pct_of_base = 100 * bh[2] / base,
             n_fires = uniqueN(p[treat==1]$fire))[]
}
