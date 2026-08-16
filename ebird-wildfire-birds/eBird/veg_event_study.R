# =============================================================================
# DID THE FIRE ACTUALLY CHANGE THE ENVIRONMENT? -- an event study on VEGETATION
#
# THE QUESTION
#   Physically, a place before and after a fire is obviously different. Yet the
#   burn block adds only +0.002 AUC once coordinates are known, and
#   yrs_since_burn ranks below every other burn feature. Is that because the
#   fire signal is genuinely weak AT THIS GRAIN -- 5 km cells, 16-day steps,
#   1 km MODIS -- or because it is redundant with information the model already
#   has?
#
#   This asks the prior question that the whole project has assumed rather than
#   checked: does fire produce a MEASURABLE environmental change in the data we
#   actually hold? If it does not, no amount of modelling will find a fire
#   effect on birds, and that is a measurement finding rather than an ecological
#   one.
#
# WHY THIS IS NOT A CLAIM ABOUT BIRDS
#   The outcome here is NDVI. No bird is involved. This is a validation of the
#   treatment, not a test of an outcome, so it does not add to the project's
#   multiple-comparison burden in the way another species x measure fit would.
#
# DESIGN
#   Same event-study logic as event_study.R, with greenness in place of
#   detection. Event time = year - fire_year, so negative values are the years
#   BEFORE that cell's own fire -- which is the "same coordinates, before"
#   comparison the frame supports through ndvi_prevyr / ndvi_delta_yr.
#   Never-burned cells provide the counterfactual trend.
#
#   Cell fixed effects: every cell is compared to ITSELF, so the chaparral
#   confound (burned cells average NDVI 0.44-0.49 vs 0.4306 never-burned --
#   burned country is greener to begin with) cannot drive the result.
#
# JOIN RULE
#   cell_fire_lookup.csv is merged on (gx, gy), NEVER on cell_id. cell_id is not
#   stable across separately built species and that bug has silently corrupted
#   three analyses in this project already.
#
# ENV  VE_PREFIX (default cfg_default_prefix())
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

.ve_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
local({
  p <- file.path(.ve_dir, "project_config.R"); if (file.exists(p)) sys.source(p, envir = globalenv())
})
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

PREFIX <- Sys.getenv("VE_PREFIX", ""); if (!nzchar(PREFIX)) PREFIX <- cfg_default_prefix()

fr <- fread(file.path(.ve_dir, sprintf("%s_cells_5km_16d.csv", PREFIX)),
            select = c("cell_id","gx","gy","year","ci","ndvi_l0","ndvi_prevyr",
                       "ndvi_delta_yr","n_checklists"), showProgress = FALSE)
fr <- fr[!is.na(ndvi_l0)]
say("frame rows with NDVI: ", format(nrow(fr), big.mark = ","))

# --- attach each cell's FIRST fire, on (gx, gy) ------------------------------
L <- fread(file.path(.ve_dir, "cell_fire_lookup.csv"), showProgress = FALSE)
setorder(L, gx, gy, fire_year)
first <- L[, .(fire_year = fire_year[1], overlap = overlap_frac[1]), by = .(gx, gy)]
n_before <- uniqueN(fr[, .(gx, gy)])
fr <- merge(fr, first, by = c("gx","gy"), all.x = TRUE)
matched <- uniqueN(fr[!is.na(fire_year), .(gx, gy)])
say(sprintf("cells %s | with a mapped fire %s (%.1f%%)",
            format(n_before, big.mark=","), format(matched, big.mark=","),
            100 * matched / n_before))
if (matched == 0) stop("no cells matched the fire lookup on (gx,gy) -- check the lookup")

fr[, evt := ifelse(is.na(fire_year), NA_integer_, year - fire_year)]

# --- raw picture: greenness by event time, heavily burned cells only ---------
HEAVY <- 0.5
say("Event-time NDVI, cells >= ", HEAVY * 100, "% burned")
b <- fr[!is.na(evt) & overlap >= HEAVY & evt >= -5L & evt <= 8L]
tab <- b[, .(n_cellyr = .N, cells = uniqueN(paste(gx, gy)),
             ndvi = mean(ndvi_l0), delta_yr = mean(ndvi_delta_yr, na.rm = TRUE)),
         by = evt][order(evt)]

ctrl <- fr[is.na(fire_year), .(ndvi_ctrl = mean(ndvi_l0),
                               delta_ctrl = mean(ndvi_delta_yr, na.rm = TRUE)), by = year]

# --- the estimate: within-cell, controlling for the year -----------------------
# NDVI ~ post + cell FE + year FE, on burned cells vs never-burned cells.
# Never-burned cells get a pseudo event year drawn from the treated distribution
# so that `post` is defined for them too -- otherwise `post` would estimate the
# common trend rather than the fire effect (a trap already logged in CLAUDE.md).
p <- fr[, .(gx, gy, year, ndvi = ndvi_l0, fire_year, overlap)]
p[, treat := as.integer(!is.na(fire_year))]
set.seed(1)
fy <- p[treat == 1, unique(fire_year)]
u  <- unique(p[treat == 0, .(gx, gy)])
u[, pseudo := sample(fy, .N, replace = TRUE)]
p <- merge(p, u, by = c("gx","gy"), all.x = TRUE)
p[, ev := fifelse(treat == 1, fire_year, pseudo)]
p[, post := as.integer(year > ev)]
p[, tp := treat * post]

fit_did <- function(dt) {
  z <- copy(dt)
  z[, uid := .GRP, by = .(gx, gy)]
  # CAST FIRST. `tp` is treat*post, i.e. INTEGER, and `:=` assigning a double
  # into an integer column TRUNCATES -- every demeaned fractional value becomes
  # 0, variance collapses to exactly zero and the fit is singular. This trap is
  # logged in CLAUDE.md and it bit again here. .prep() in infer.R does the same
  # cast for the same reason.
  for (cc in c("ndvi","tp")) set(z, NULL, cc, as.numeric(z[[cc]]))
  for (i in seq_len(40)) for (cc in c("ndvi","tp")) {
    z[, (cc) := get(cc) - mean(get(cc)), by = uid]
    z[, (cc) := get(cc) - mean(get(cc)), by = year]
  }
  X <- as.matrix(z[, .(tp)]); y <- z$ndvi
  b <- solve(crossprod(X), crossprod(X, y))
  e <- y - X %*% b
  # cluster by cell -- deliberately the loose version; the point here is the
  # magnitude and sign of an environmental change, not a significance claim
  S <- rowsum(as.vector(X) * as.vector(e), z$uid)
  se <- sqrt(solve(crossprod(X)) %*% crossprod(S) %*% solve(crossprod(X)))
  c(est = as.numeric(b), se = as.numeric(se))
}

res <- rbindlist(lapply(c(0, 0.25, 0.5, 0.9), function(th) {
  d <- p[treat == 0 | overlap >= th]
  r <- fit_did(d)
  data.table(threshold = th, est = r[["est"]], se = r[["se"]],
             t = r[["est"]] / r[["se"]],
             treated_cells = uniqueN(d[treat == 1, paste(gx, gy)]))
}))

base_pre <- p[treat == 1 & post == 0, mean(ndvi)]

out <- file.path(.ve_dir, "veg_event_study_summary.txt")
sink(out)
cat("DID THE FIRE ACTUALLY CHANGE THE ENVIRONMENT?\n")
cat("=============================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Species   : ", PREFIX, " (cells only -- the OUTCOME HERE IS NDVI, not birds)\n", sep = "")
cat("Design    : within-cell, cell + year fixed effects, treat x post\n")
cat("Purpose   : validate that fire produces a measurable environmental change\n")
cat("            at 5 km / 16-day / 1 km MODIS resolution.\n\n")

cat("-- 1. GREENNESS BY YEAR RELATIVE TO EACH CELL'S OWN FIRE (>= 50% burned) --\n")
cat("   evt < 0 is BEFORE that cell's fire, at the SAME coordinates.\n\n")
print(as.data.frame(tab[, .(evt, cells, cell_years = n_cellyr,
                            ndvi = round(ndvi, 4),
                            vs_prev_yr = round(delta_yr, 4))]), row.names = FALSE)

cat("\n-- 2. WITHIN-CELL DiD ON NDVI ---------------------------------------\n")
cat("   pre-fire mean NDVI in treated cells: ", round(base_pre, 4), "\n\n", sep = "")
print(as.data.frame(res[, .(threshold, estimate = round(est, 4),
                            se = round(se, 4), t = round(t, 2),
                            pct_of_pre = round(100 * est / base_pre, 1),
                            treated_cells)]), row.names = FALSE)

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  A clearly negative estimate means fire IS visible in the vegetation at\n")
cat("  this grain, and the weak burn-block importance is about REDUNDANCY --\n")
cat("  the model already reads the same change through NDVI.\n")
cat("  An estimate near zero means fire is NOT measurable at this resolution,\n")
cat("  which would be a measurement finding and would explain the weak fire\n")
cat("  signal throughout the project.\n")
cat("\n  Cell-clustered SEs here are the LOOSE version and are known to be too\n")
cat("  narrow in this project. This is a magnitude check, not a significance\n")
cat("  claim, and no bird is involved.\n")
sink()
cat(readLines(out), sep = "\n")
say("DONE -> veg_event_study_summary.txt")
