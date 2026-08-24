# =============================================================================
# GAP 1: does wildfire SEVERITY (not just extent) affect occurrence?
#
# The research question asks about satellite-derived wildfire SEVERITY. Every
# analysis so far has used burn EXTENT -- how much of a cell burned -- which is
# a different thing. A cell can be 100% burned at low severity or 30% burned to
# bare mineral soil, and dose_burn_effect.R cannot tell those apart.
#
# The one earlier attempt at a severity contrast FAILED: `max_sev_prior` records
# the worst patch anywhere in a 5 km cell, so in any heavily burned cell it
# saturates at class 4 (285 of 290 cells). Useless as a gradient.
#
# THE FIX. The regrid cache retains per-(cell, fire_year, severity) burned
# fractions -- exactly what is needed. Instead of one dose, this decomposes the
# burned area into FOUR doses, one per MTBS severity class:
#
#   rate_it = a_i + b_t + p*post_it + SUM_s g_s * (post_it * frac_s_i) + e_it
#
# g_s is the effect per unit of cell area burned AT SEVERITY s. Comparing g_1
# (unburned-low) with g_4 (high) is the severity gradient the question asks for,
# holding extent fixed: two cells equally burned but at different severities get
# different predictions.
#
# MTBS classes: 1 unburned-low, 2 low, 3 moderate, 4 high, 5 increased
# greenness, 6 mask. Only 1-4 form the severity ladder; 5 is regrowth and 6 is
# cloud/water, both excluded from the dose (they are 0.4% of burned area).
#
# Inference: cluster-robust by cell. Outcome det_rate (effort-normalised).
# Design otherwise identical to dose_burn_effect.R so results are comparable.
#
# WRITES NEW FILES ONLY. Reads the existing frames and regrid cache read-only.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
out_sum <- file.path(dir_, "severity_effect_summary.txt")

# Species from the registry, not a literal list -- see event_study.R for why.
sp_reg <- cfg_species(active_only = TRUE)
sp_reg <- sp_reg[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]
if (!nrow(sp_reg)) stop("no active species in species_registry.csv has a built frame")
SPP <- setNames(as.list(sprintf("%s (%s)", sp_reg$common_name, sp_reg$habitat)),
                sp_reg$prefix)
SEV_LAB <- c("1 unburned-low", "2 low", "3 moderate", "4 high")

# --- the two models, declared in treatment_registry.csv ------------------------
# Model A ("class_doses") is one dose per MTBS severity class; Model B
# ("extent_x_meansev") is extent plus extent x mean severity. Both were written
# out inline here, duplicating rows that already exist in the registry -- so the
# registry could be edited without this script noticing. Now the terms ARE the
# registry's. The defaults reproduce the published numbers exactly.
SEV_A <- Sys.getenv("SEV_DOSE_A", "class_doses")
SEV_B <- Sys.getenv("SEV_DOSE_B", "extent_x_meansev")
sev_terms <- function(nm) {
  tr <- cfg_treatments()[name == nm]
  if (!nrow(tr)) stop("burn measure '", nm, "' is not in treatment_registry.csv")
  if (tr$source[1] != "cache")
    stop("severity_effect.R needs a `cache`-sourced measure (per-class fractions); '",
         nm, "' is sourced from ", tr$source[1])
  .parse_terms(tr$terms[1])
}
TERMS_A <- sev_terms(SEV_A); TERMS_B <- sev_terms(SEV_B)

# display labels: the severity ladder keeps its words, anything else is named
# after its registry term so a new measure prints sensibly without an edit here
sev_label <- function(tn) {
  if (tn %in% paste0("s", 1:4)) return(SEV_LAB[as.integer(sub("s", "", tn))])
  switch(tn, extent = "post x extent", extent_sev = "post x extent x mean-sev",
         paste("post x", tn))
}

# --- within estimator, cluster-robust SEs -------------------------------------
twfe_cl <- function(dt, xcols, iters = 40) {
  z <- copy(dt)[, c("cell_id","year","rate","w", xcols), with = FALSE]
  cl <- z$cell_id
  for (cc in c("rate", xcols)) set(z, NULL, cc, as.numeric(z[[cc]]))
  for (i in seq_len(iters)) for (cc in c("rate", xcols)) {
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = cell_id]
    z[, (cc) := get(cc) - weighted.mean(get(cc), w), by = year]
  }
  X <- as.matrix(z[, ..xcols]); y <- z$rate; w <- z$w
  ok <- apply(X, 2, function(v) sum(abs(v)) > 1e-12)
  X <- X[, ok, drop = FALSE]
  XtWX <- t(X * w) %*% X
  if (abs(det(XtWX)) < 1e-26) return(NULL)
  Xi <- solve(XtWX); b <- as.vector(Xi %*% (t(X * w) %*% y))
  e <- as.vector(y - X %*% b); U <- X * (w * e)
  V <- Xi %*% (t(rowsum(U, cl)) %*% rowsum(U, cl)) %*% Xi
  data.table(term = colnames(X), beta = b, se = sqrt(diag(V)))
}

# --- per-cell severity composition of the in-window fire ----------------------
# Uses panel_utils::.cache_sev(), NOT a direct readRDS of a fixed filename.
# The cache is keyed by cell_id, which is NOT stable across separately built
# species, and the regrid now writes one cache per numbering (burn_<tag>.rds).
# Reading ".regrid_cache/burn.rds" directly would pick an arbitrary numbering and
# silently mis-join; .cache_sev(prefix) selects the matching one and errors if
# none does. It is also the identical construction, so the two cannot drift.
# Computed INSIDE the species loop for that reason -- it is per-numbering, and
# hoisting it out is what made it look species-independent.
results <- list(); meta <- list(); comp_tab <- NULL
for (pre in names(SPP)) {
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  if (!file.exists(f)) { say("missing ", basename(f)); next }
  say("=== ", SPP[[pre]], " ===")
  d <- fread(f, showProgress = FALSE,
             select = c("cell_id","year","ever_burned","n_checklists","n_detections"))
  # Treated/control from the shared definition (panel_utils.R) so this agrees
  # with dose_burn_effect.R by construction. The severity mix still comes from
  # sev_comp, which is the same construction as panel_utils::.cache_sev().
  ct <- cell_treatment(pre, dose = "frac_max")
  treated <- ct[treat == 1]$cell_id
  never   <- ct[treat == 0]$cell_id
  ev <- ct[treat == 1, .(cell_id, ev_year)]

  # Matched to THIS species' cell numbering. .cache_sev() errors rather than
  # silently using another species' cache -- catch it so one species without a
  # regrid cache skips instead of aborting the whole multi-species run.
  sev_comp <- tryCatch(.cache_sev(pre), error = function(e) e)
  if (inherits(sev_comp, "error")) {
    say("  no regrid burn cache for '", pre, "' -- SKIPPED. Fix with:")
    say("    EBIRD_PREFIX=", pre, " Rscript regrid_cells.R   (~6 min)")
    next
  }
  say("  cells with an in-window fire: ", nrow(sev_comp))

  tr <- merge(ev, sev_comp, by = "cell_id")          # treated cells with severity mix
  say(sprintf("  treated with severity data: %d | controls: %d", nrow(tr), length(never)))
  if (is.null(comp_tab)) comp_tab <- tr

  cy <- rbind(
    merge(d[cell_id %in% tr$cell_id], tr, by = "cell_id"
          )[, .(det = sum(n_detections), chk = sum(n_checklists)),
            by = .(cell_id, year, ev_year, s1, s2, s3, s4, frac_tot, mean_sev)],
    d[cell_id %in% never, .(det = sum(n_detections), chk = sum(n_checklists)),
      by = .(cell_id, year)][, `:=`(ev_year = NA_integer_, s1 = 0, s2 = 0, s3 = 0,
                                    s4 = 0, frac_tot = 0, mean_sev = 0)],
    use.names = TRUE)
  cy <- cy[chk > 0]
  cy[, `:=`(rate = det / chk, w = as.numeric(chk),
            post = as.numeric(!is.na(ev_year) & year >= ev_year))]
  # Evaluate each registry term against the cell-year table, then interact it
  # with post. Controls carry s1..s4 / frac_tot / mean_sev = 0, so every term is
  # 0 for them, exactly as the hardcoded version did.
  mk <- function(terms) {
    for (tn in names(terms))
      set(cy, NULL, paste0("p_", tn),
          as.numeric(cy$post) * as.numeric(eval(parse(text = terms[[tn]]), envir = cy)))
    paste0("p_", names(terms))
  }
  # mk() MUST run before twfe_cl(), not inside its argument list: R forces
  # arguments lazily, so `copy(dt)` inside twfe_cl was evaluated first and copied
  # cy before mk() had added the p_* columns to it -- "columns not found".
  xA <- mk(TERMS_A)   # Model A: one dose per severity class (default class_doses)
  xB <- mk(TERMS_B)   # Model B: extent + extent x mean severity (extent_x_meansev)
  rA <- twfe_cl(cy, c("post", xA))
  rB <- twfe_cl(cy, c("post", xB))
  base <- cy[post == 0 & !is.na(ev_year), sum(det) / sum(chk)]
  results[[pre]] <- list(A = rA, B = rB)
  meta[[pre]] <- data.table(species = pre, label = SPP[[pre]],
                            treated = nrow(tr), controls = length(never), base = base)
  say(sprintf("  done (pre-fire baseline %.4f)", base))
}
M <- rbindlist(meta)

# =============================================================================
sink(out_sum)
cat("Wildfire SEVERITY and bird occurrence\n")
cat("=====================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Question  : does severity matter, holding burned EXTENT fixed?\n")
cat("Outcome   : det_rate = detections / checklists (effort-normalised)\n")
cat("Design    : within-cell before/after vs never-burned controls,\n")
cat("            cell + year fixed effects, cluster-robust SEs by cell\n")
cat("Severity  : MTBS classes 1 unburned-low, 2 low, 3 moderate, 4 high\n")
cat("            (5 increased-greenness and 6 mask excluded from the dose)\n\n")
print(as.data.frame(M[, .(species = label, treated_cells = treated,
                          control_cells = controls, pre_fire_rate = round(base, 4))]),
      row.names = FALSE)

cat("\n-- BURNED-AREA COMPOSITION OF TREATED CELLS --------------------------\n")
cat("  (identical across species -- same cells, same fires)\n")
cc <- comp_tab
cat(sprintf("  mean fraction of cell burned at each severity, over %d treated cells:\n", nrow(cc)))
for (k in 1:4)
  cat(sprintf("    %-16s %.4f   (cells with any: %d)\n", SEV_LAB[k],
              mean(cc[[paste0("s", k)]]), sum(cc[[paste0("s", k)]] > 0)))
cat(sprintf("    %-16s %.4f\n", "TOTAL burned", mean(cc$frac_tot)))
cat(sprintf("    area-weighted mean severity: %.2f (range %.2f - %.2f)\n",
            mean(cc$mean_sev), min(cc$mean_sev), max(cc$mean_sev)))
cat("  NOTE this is why max_sev_prior failed: nearly every treated cell contains\n")
cat("  SOME class-4 patch. The fractions vary; the maximum does not.\n")

cat("\n-- MODEL A: one dose per severity class ------------------------------\n")
cat("  beta_s = change in detection rate per unit of cell area burned at s.\n")
cat("  If severity matters, |beta| should GROW from class 1 to class 4.\n")
for (pre in names(SPP)) {
  if (is.null(results[[pre]])) next
  r <- results[[pre]]$A; b <- M[species == pre]$base
  cat("\n  === ", SPP[[pre]], " (baseline ", sprintf("%.4f", b), ") ===\n", sep = "")
  cat(sprintf("  %-18s %10s %9s %22s %9s\n", "term", "beta", "se", "95% CI", "%base"))
  for (i in seq_len(nrow(r))) {
    lo <- r$beta[i] - 1.96*r$se[i]; hi <- r$beta[i] + 1.96*r$se[i]
    lab <- if (r$term[i] == "post") "post (extent 0)" else
           sev_label(sub("^p_", "", r$term[i]))
    cat(sprintf("  %-18s %+10.4f %9.4f  [%+.4f, %+.4f] %+8.1f%% %s\n",
                lab, r$beta[i], r$se[i], lo, hi, 100*r$beta[i]/b,
                if (lo > 0 || hi < 0) "*" else ""))
  }
}

cat("\n-- MODEL B: does severity add anything OVER extent? ------------------\n")
cat("  post x extent            = effect per unit burned area, severity ignored\n")
cat("  post x extent x mean-sev = EXTRA effect per unit of severity on top\n")
cat("  If the second row's CI excludes 0, severity carries information that\n")
cat("  extent alone does not.\n")
for (pre in names(SPP)) {
  if (is.null(results[[pre]])) next
  r <- results[[pre]]$B; b <- M[species == pre]$base
  cat("\n  === ", SPP[[pre]], " ===\n", sep = "")
  for (i in seq_len(nrow(r))) {
    lo <- r$beta[i] - 1.96*r$se[i]; hi <- r$beta[i] + 1.96*r$se[i]
    lab <- if (r$term[i] == "post") "post (common trend)" else
           sev_label(sub("^p_", "", r$term[i]))
    cat(sprintf("  %-28s %+9.4f  [%+.4f, %+.4f] %s\n", lab, r$beta[i], lo, hi,
                if (lo > 0 || hi < 0) "*" else ""))
  }
}

cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat("  * Severity classes are correlated within a cell (a big fire produces all\n")
cat("    four), so the four coefficients are not cleanly separable -- read the\n")
cat("    GRADIENT across them, not each one alone.\n")
cat("  * MTBS severity is thematic (dnbr6), not continuous dNBR; the class\n")
cat("    thresholds are per-fire calibrated, so class 3 in one fire is not\n")
cat("    exactly class 3 in another.\n")
cat("  * Same limitations as the dose analysis: 5 km aggregation, MTBS\n")
cat("    under-mapping 2023-2025, detection is not occupancy.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
