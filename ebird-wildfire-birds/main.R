# =============================================================================
# main.R -- one entry point for the whole replication.
#
# Open this repository as an RStudio project (or set the working directory to
# the repository root) and source this file. Do NOT call setwd() inside it.
#
#   Rscript main.R              # everything, in order  (~12-15 h)
#   Rscript main.R build        # data build only        (~10-13 h)
#   Rscript main.R analysis     # analyses only, assumes the build exists
#   Rscript main.R figures      # figures only
#   Rscript main.R verify       # regression checks only
#
# Stages are independent once the build exists, so a failed analysis does not
# force a rebuild. See README section 5.
# =============================================================================

STAGE <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(STAGE)) STAGE <- "all"
stopifnot(STAGE %in% c("all", "build", "analysis", "figures", "verify"))

# --- where the data lives ----------------------------------------------------
# Scripts resolve their inputs and outputs as $EBIRD_PROJ_ROOT/eBird. The
# repository root works as a default only if the downloads were unpacked here.
PROJ <- Sys.getenv("EBIRD_PROJ_ROOT", unset = normalizePath(".", winslash = "/"))
CODE <- file.path(PROJ, "eBird")

if (!dir.exists(CODE)) {
  stop("No eBird/ directory under '", PROJ, "'.\n",
       "  Set EBIRD_PROJ_ROOT to the folder that contains eBird/, e.g.\n",
       "  export EBIRD_PROJ_ROOT=/path/to/project", call. = FALSE)
}

setwd(CODE)   # the one permitted setwd: scripts are written to run from here

say <- function(...) cat(sprintf(...), "\n")
run <- function(script, note, halt_on_fail = TRUE) {
  say("\n--- %s  (%s)", script, note)
  t0 <- Sys.time()
  ok <- tryCatch({ system2("Rscript", script) == 0L },
                 error = function(e) { message(conditionMessage(e)); FALSE })
  say("    %s after %.1f min", if (ok) "done" else "returned non-zero",
      as.numeric(difftime(Sys.time(), t0, units = "mins")))
  if (!ok && halt_on_fail) stop("stopped at ", script, call. = FALSE)
  invisible(ok)
}

# --- preflight ---------------------------------------------------------------
# Fail here with a readable message rather than three hours in with a cryptic
# one. The model packages are checked separately because they are only needed
# by the prediction scripts.
need <- c("data.table", "terra", "sf", "lubridate", "auk",
          "ggplot2", "patchwork", "ragg")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop("Missing packages: ", paste(miss, collapse = ", "),
       "\n  install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))",
       call. = FALSE)
}

if (!any(grepl("HDF4", terra::gdal(drivers = TRUE)$name))) {
  warning("GDAL has no HDF4 driver -- the MODIS granules will not open.\n",
          "  The eBird and MTBS steps will still run.", call. = FALSE)
}

say("R           : %s", R.version.string)
say("project root: %s", PROJ)
say("stage       : %s", STAGE)

# --- stages ------------------------------------------------------------------

if (STAGE %in% c("all", "build")) {
  say("\n=== BUILD =====================================================")
  run("zerofill_from_registry.R", "raw eBird -> presence/absence, hours")
  run("join_species_all_years.R", "join to fire and vegetation, ~9 min")
  run("regrid_cells.R",           "aggregate to 5 km x 16 day, ~6 min")
}

if (STAGE %in% c("all", "analysis")) {
  say("\n=== ANALYSIS ==================================================")
  # Run one at a time. Two concurrent jobs will exhaust memory on a 24 GB box.
  run("dose_burn_effect.R", "fire as a dose, ~6 min")
  run("event_study.R",      "year by year, ~6 min")
  run("severity_effect.R",  "severity holding extent fixed, ~6 min")
  run("effort_placebo.R",   "did fire change how people birded, ~1 min")
  run("wild_cluster_boot.R","significance, clustered by fire, ~35 min")
  run("guild_contrast.R",   "group comparisons, ~35 min")
  run("sign_test.R",        "predicted vs observed direction, ~1 min")
}

if (STAGE %in% c("all", "figures")) {
  say("\n=== FIGURES ===================================================")
  run("make_figures.R", "figures/fig[0-6]_*.png and .pdf")
}

if (STAGE %in% c("all", "verify")) {
  say("\n=== VERIFY ====================================================")
  # Expect 40 of 41, so a non-zero exit is NOT a reason to halt here. Check 6c
  # is a known false positive: it looks for scripts calling cfg_*() without
  # sourcing a registry, and does not recognise sys.source(). Read the failure
  # list -- if 6c is the only entry, the run is clean. See README section 5.
  run("verify_published.R", "41 published numbers; expect 40 to pass",
      halt_on_fail = FALSE)
  say("\n    Expected: 40 of 41, with check 6c the only failure.")
  say("    Any OTHER failure means a published number moved -- stop and find out why.")
}

say("\nFinished stage '%s'.", STAGE)
say("Session info written to results/sessionInfo.txt")
dir.create(file.path(PROJ, "results"), showWarnings = FALSE)
capture.output(sessionInfo(), file = file.path(PROJ, "results", "sessionInfo.txt"))
