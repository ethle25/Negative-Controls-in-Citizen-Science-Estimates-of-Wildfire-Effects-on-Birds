# =============================================================================
# SEVERITY, REFIT WITH FIRE CLUSTERING -- the test the best result has not faced.
#
# WHY
#   severity_effect.R clusters by CELL. This project has already established
#   that cell clustering is too narrow: cells inside one perimeter share one
#   ignition, one weather system, one fuel type and one post-fire access regime
#   -- one treatment draw, not eight. When the dose-response was refit with FIRE
#   clustering and a wild cluster bootstrap, NO species survived.
#
#   The severity result is the project's strongest -- 4 of 6 snag species
#   significant, all positive, and the three negative controls silent -- and it
#   has never faced that test. This runs it.
#
# WHAT CHANGES AND WHAT MUST NOT
#   ONLY the standard errors. The panel and the point estimates are reproduced
#   exactly as severity_effect.R builds them, and the script asserts that Wrentit
#   still lands on its published -0.0458. If a point estimate moves, the panel
#   was rebuilt wrong and the run is void.
#
# INFERENCE
#   Cluster-robust by fire (crve) plus the wild cluster bootstrap (wcb_p),
#   Rademacher signs, null imposed -- the same instrument that overturned the
#   percentile verdicts for the dose-response.
#
# ENV
#   SFC_B   wild cluster replications (default 999; use 9999 to quote)
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
source(file.path(dir_, "infer.R"))
B_REPS  <- as.integer(Sys.getenv("SFC_B", "999"))
out_sum <- file.path(dir_, "severity_fire_cluster_summary.txt")

sp <- cfg_species(active_only = FALSE)
sp <- sp[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]

res <- rbindlist(lapply(seq_len(nrow(sp)), function(i) {
  pre <- sp$prefix[i]
  say("=== ", sp$common_name[i], " (", pre, ")")
  f <- file.path(dir_, sprintf("%s_cells_5km_16d.csv", pre))
  d <- fread(f, select = c("cell_id","year","ever_burned","n_checklists","n_detections"),
             showProgress = FALSE)

  sev <- tryCatch(.cache_sev(pre), error = function(e) e)
  if (inherits(sev, "error")) { say("  no regrid cache -- skipped"); return(NULL) }

  ct <- cell_treatment(pre, dose = "frac_max")
  never <- ct[treat == 0]$cell_id
  tr <- merge(ct[treat == 1, .(cell_id, ev_year)], sev, by = "cell_id")

  # EXACTLY severity_effect.R's panel: controls carry ev_year = NA and therefore
  # post = 0 throughout. That is a different convention from build_panel()'s
  # pseudo event year, and changing it would move the point estimates.
  cy <- rbind(
    merge(d[cell_id %in% tr$cell_id], tr, by = "cell_id"
          )[, .(det = sum(n_detections), chk = sum(n_checklists)),
            by = .(cell_id, year, ev_year, s1, s2, s3, s4, frac_tot, mean_sev)],
    d[cell_id %in% never, .(det = sum(n_detections), chk = sum(n_checklists)),
      by = .(cell_id, year)][, `:=`(ev_year = NA_integer_, s1 = 0, s2 = 0, s3 = 0,
                                    s4 = 0, frac_tot = 0, mean_sev = 0)],
    use.names = TRUE)[chk > 0]
  cy[, `:=`(rate = det / chk, w = as.numeric(chk),
            post = as.numeric(!is.na(ev_year) & year >= ev_year))]
  cy[, `:=`(p_frac = post * frac_tot, p_msev = post * frac_tot * mean_sev,
            treat = as.integer(!is.na(ev_year)))]

  cyf <- attach_fire(cy, pre)          # joins on (gx,gy); adds `fire` and `cl`
  xc  <- c("post","p_frac","p_msev")
  P   <- .prep(cyf, xc)
  b   <- .fit(P, P$y)
  V   <- crve(P, as.vector(P$y - P$X %*% b))

  rbindlist(lapply(2:3, function(k) data.table(
    species = sp$common_name[i], prefix = pre, guild = sp$guild[i],
    term    = c("extent","extent x mean-sev")[k - 1],
    estimate = b[k], se_fire = sqrt(V[k,k]), t_fire = b[k]/sqrt(V[k,k]),
    p_wcb   = wcb_p(P, k, b[k], sqrt(V[k,k]), B = B_REPS),
    n_fires = uniqueN(cyf[treat == 1]$fire))))
}))

# --- the integrity check: point estimates must not have moved ------------------
w <- res[prefix == "wrentit" & term == "extent"]$estimate
if (length(w) && abs(w - (-0.0458)) > 5e-4)
  stop("Wrentit extent is ", round(w,4), ", published is -0.0458. The panel was ",
       "rebuilt differently -- these standard errors describe a different model. VOID.")
say("integrity check passed: Wrentit extent reproduces -0.0458")

fwrite(res, file.path(dir_, "severity_fire_cluster_results.csv"))

sev <- res[term == "extent x mean-sev"]
sink(out_sum)
cat("SEVERITY REFIT WITH FIRE CLUSTERING\n===================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Change    : ONLY the standard errors. Point estimates are severity_effect.R's.\n")
cat("Inference : cluster-robust by FIRE + wild cluster bootstrap, B = ", B_REPS, "\n\n", sep="")
cat("-- SEVERITY TERM (post x extent x mean-severity) --------------------\n")
cat("   The extra effect of burning MORE SEVERELY, holding extent fixed.\n\n")
print(as.data.frame(sev[order(guild, -estimate),
  .(species, guild, est = round(estimate, 4), se = round(se_fire, 4),
    t = round(t_fire, 2), p = round(p_wcb, 4),
    sig = fifelse(p_wcb < 0.05, "*", ""), fires = n_fires)]), row.names = FALSE)

cat("\n-- BY GUILD ---------------------------------------------------------\n\n")
g <- sev[, .(n = .N, n_sig = sum(p_wcb < 0.05), mean_est = round(mean(estimate), 4)), by = guild]
print(as.data.frame(g[order(-mean_est)]), row.names = FALSE)

cat("\n-- THE CONTROLS -----------------------------------------------------\n")
ctl <- sev[guild == "control"]
if (nrow(ctl)) {
  print(as.data.frame(ctl[, .(species, est = round(estimate,4), p = round(p_wcb,4),
                              sig = fifelse(p_wcb < 0.05, "*** PROBLEM ***", "silent"))]),
        row.names = FALSE)
}
cat("\n-- VERDICT ----------------------------------------------------------\n")
snag <- sev[guild == "snag"]
cat(sprintf("  snag guild significant under FIRE clustering : %d of %d\n",
            sum(snag$p_wcb < 0.05), nrow(snag)))
cat(sprintf("  controls significant                        : %d of %d\n",
            sum(ctl$p_wcb < 0.05), nrow(ctl)))
cat("\n  Under CELL clustering it was 4 of 6 snag, 0 of 3 controls.\n")
cat("  If the snag count collapses here, the severity result is a\n")
cat("  clustering artifact like every earlier percentile verdict was.\n")
cat("  If it holds, it is the first result in this project to survive\n")
cat("  both a negative control AND fire-level inference.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
