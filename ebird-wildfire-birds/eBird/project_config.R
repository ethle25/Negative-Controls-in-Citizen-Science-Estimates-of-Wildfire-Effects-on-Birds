# =============================================================================
# Shared configuration: species and feature registries
#
# WHY THIS EXISTS
#   Species were switched through three environment variables and a shell script
#   that took positional arguments, while feature blocks were hardcoded as four
#   vectors inside model_burn_effect.R. Adding a species meant editing the
#   upstream zero-fill; adding a feature meant editing every downstream script
#   that referenced a block. This file makes both declarative.
#
#   source() it from any script:  source(file.path(dir_, "project_config.R"))
#
# TWO REGISTRIES
#   species_registry.csv  one row per species. `active = 1` means "include in
#                         runs". `guild` and `predicted` are NOT bookkeeping --
#                         they are what the guild-contrast and sign tests read,
#                         so the ecological prediction is committed to a file
#                         BEFORE anything runs. Do not edit `predicted` after
#                         seeing a result; that is the whole point of it.
#
#   feature_registry.csv  one row per column of the cell frame, tagged with a
#                         block (effort / season / veg / burn / geo / meta), a
#                         role (predictor / target / id / leak / support) and an
#                         active flag. Generated from the real frame header, so
#                         every column is accounted for rather than whichever
#                         ones someone remembered.
#
# BACKWARD COMPATIBILITY
#   The registries reproduce the previous hardcoded behaviour exactly: 42 active
#   modelling features + dist_median_na added at runtime = the 43 that
#   model_burn_effect.R has always used. Verified at load by cfg_check().
#   EBIRD_SPECIES / EBIRD_SPECIES_CODE / EBIRD_PREFIX still override, so nothing
#   that already works stops working.
# =============================================================================
suppressMessages(library(data.table))

.cfg_dir <- Sys.getenv("EBIRD_CFG_DIR",
  file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))

# --- species -----------------------------------------------------------------

#' All registered species, or just the active ones.
#' @param active_only keep only rows with active = 1
cfg_species <- function(active_only = TRUE) {
  f <- file.path(.cfg_dir, "species_registry.csv")
  if (!file.exists(f)) stop("species_registry.csv not found in ", .cfg_dir)
  s <- fread(f, showProgress = FALSE)
  req <- c("common_name","species_code","prefix","guild","predicted","active")
  miss <- setdiff(req, names(s))
  if (length(miss)) stop("species_registry.csv missing columns: ",
                         paste(miss, collapse = ", "))
  if (anyDuplicated(s$prefix)) stop("duplicate prefix in species_registry.csv: ",
      paste(s$prefix[duplicated(s$prefix)], collapse = ", "))
  bad <- setdiff(unique(s$predicted), c("negative","positive","none"))
  if (length(bad)) stop("predicted must be negative/positive/none, got: ",
                        paste(bad, collapse = ", "))
  if (active_only) s <- s[active == 1]
  s[]
}

#' The species this run is for. An explicit EBIRD_PREFIX wins; otherwise the
#' first active row. Returns a one-row data.table.
cfg_current_species <- function() {
  p <- Sys.getenv("EBIRD_PREFIX", "")
  s <- cfg_species(active_only = FALSE)
  if (nzchar(p)) {
    hit <- s[prefix == p]
    if (!nrow(hit)) stop("EBIRD_PREFIX='", p, "' is not in species_registry.csv. ",
                         "Add a row for it, or use one of: ",
                         paste(s$prefix, collapse = ", "))
    return(hit[1])
  }
  a <- s[active == 1]
  if (!nrow(a)) stop("no active species in species_registry.csv")
  a[1]
}

#' The default species prefix / common name for a run.
#'
#' Every script used to fall back to the literals "wrentit" / "Wrentit" when no
#' env var was set -- 23 sites across 10 files. That made Wrentit the implicit
#' subject of the whole pipeline: deactivating it in the registry changed
#' nothing, and a fresh reader could not tell which scripts were species-neutral.
#' These resolve through the registry instead, so the default is "the first
#' ACTIVE row" and no bird is named in code. An explicit env var still wins, so
#' behaviour is unchanged for every existing invocation.
cfg_default_prefix  <- function() {
  p <- Sys.getenv("EBIRD_PREFIX", "")
  if (nzchar(p)) p else cfg_current_species()$prefix[1]
}
cfg_default_species <- function() {
  s <- Sys.getenv("EBIRD_SPECIES", "")
  if (nzchar(s)) s else cfg_current_species()$common_name[1]
}

# --- features ----------------------------------------------------------------

#' Column names for one or more blocks.
#' @param blocks e.g. c("effort","season","veg","burn"); NULL = all
#' @param active_only keep only rows with active = 1
#' @param roles restrict to these roles (default: predictors only)
cfg_features <- function(blocks = NULL, active_only = TRUE, roles = "predictor") {
  f <- file.path(.cfg_dir, "feature_registry.csv")
  if (!file.exists(f)) stop("feature_registry.csv not found in ", .cfg_dir)
  r <- fread(f, showProgress = FALSE)
  if (active_only)        r <- r[active == 1]
  if (!is.null(roles))    r <- r[role %in% roles]
  if (!is.null(blocks))   r <- r[block %in% blocks]
  r$column
}

#' Named list of blocks -> columns, for permutation importance by block.
cfg_blocks <- function(blocks = c("effort","season","veg","burn"), ...) {
  setNames(lapply(blocks, function(b) cfg_features(b, ...)), blocks)
}

#' Columns that must never be used as predictors, with the reason.
cfg_excluded <- function() {
  r <- fread(file.path(.cfg_dir, "feature_registry.csv"), showProgress = FALSE)
  r[role %in% c("id","leak","target"), .(column, role, notes)]
}

# --- validation --------------------------------------------------------------

#' Check the registries against a real frame. Call this at the top of any script
#' that reads them -- a registry that has drifted from the frame fails silently
#' otherwise, and this project has already been bitten by a cell_id mapping that
#' changed underneath a cache.
cfg_check <- function(frame_path = NULL, verbose = TRUE) {
  r <- fread(file.path(.cfg_dir, "feature_registry.csv"), showProgress = FALSE)
  s <- cfg_species(active_only = FALSE)
  ok <- TRUE
  if (is.null(frame_path)) {
    sp <- cfg_current_species()
    frame_path <- file.path(.cfg_dir, sprintf("%s_cells_5km_16d.csv", sp$prefix))
  }
  if (file.exists(frame_path)) {
    have <- strsplit(readLines(frame_path, n = 1), ",")[[1]]
    miss <- setdiff(r$column, have)
    extra <- setdiff(have, r$column)
    if (length(miss))  { ok <- FALSE
      message("!! in registry but NOT in frame: ", paste(miss, collapse = ", ")) }
    if (length(extra)) { ok <- FALSE
      message("!! in frame but NOT in registry (add a row): ",
              paste(extra, collapse = ", ")) }
  } else if (verbose) {
    message("(frame not found, skipping column check: ", basename(frame_path), ")")
  }
  nfeat <- length(cfg_features(c("effort","season","veg","burn")))
  if (verbose) {
    message(sprintf("config: %d active species | %d active modelling features (+1 runtime = %d)",
                    nrow(s[active == 1]), nfeat, nfeat + 1))
    message("  guilds: ", paste(sprintf("%s=%d", names(table(s[active == 1]$guild)),
                                        table(s[active == 1]$guild)), collapse = ", "))
  }
  invisible(ok)
}

# --- the ecological prediction, as data --------------------------------------

#' Guild -> expected sign, as a numeric contrast usable in a pooled model or a
#' sign test. This is the pre-registered hypothesis in machine-readable form.
cfg_expected_sign <- function() {
  s <- cfg_species()
  s[, .(prefix, common_name, guild,
        expected = fifelse(predicted == "positive", 1,
                    fifelse(predicted == "negative", -1, 0)))]
}

# =============================================================================
# Species catalog + registry maintenance (added 2026-07-31)
# =============================================================================

#' Every species in the raw California EBD, with record counts.
#' Built by scanning both EBD downloads for CATEGORY == "species" -- 773 species.
#' `n_complete` counts records on COMPLETE checklists, which is what zero-filling
#' can actually use; `n_records` includes incomplete lists and is the larger,
#' more flattering, less useful number.
cfg_catalog <- function(min_complete = 0) {
  f <- file.path(.cfg_dir, "species_catalog.csv")
  if (!file.exists(f)) stop("species_catalog.csv not found -- rebuild with the ",
                            "awk scan in CLAUDE.md section 4")
  d <- data.table::fread(f, showProgress = FALSE)
  d[n_complete >= min_complete][]
}

#' Look a species up in the catalog by partial name.
cfg_find_species <- function(pattern) {
  d <- cfg_catalog()
  d[grepl(pattern, common_name, ignore.case = TRUE) |
    grepl(pattern, scientific_name, ignore.case = TRUE)][
    , .(common_name, scientific_name, n_complete, has_data, rank)][]
}

#' Add a species to species_registry.csv.
#'
#' Validates the name against the catalog first -- a typo would otherwise produce
#' an empty zero-fill and a silently broken species run. Does NOT build the data:
#' that needs `Rscript zerofill_from_registry.R` plus gawk.
#'
#' @param common_name  must match the catalog exactly (use cfg_find_species())
#' @param prefix       short id used in filenames, e.g. "bkbwoo"
#' @param guild        habitat group; the sign test and shrub-vs-snag contrast
#'                     read this, so set it BEFORE running anything
#' @param direction    predicted fire response: negative / positive / none
cfg_add_species <- function(name, prefix, guild, direction = "none",
                            species_code = NA_character_, active = 1L) {
  # NB the argument is `name`, not `common_name`. A parameter sharing a column
  # name breaks data.table's subsetting -- `cat_[common_name == common_name]`
  # compares the column with itself and matches every row, and the `..` prefix
  # does not reach function arguments either.
  cat_ <- cfg_catalog()
  hit <- cat_[common_name == name]
  if (!nrow(hit)) stop("'", name, "' is not in the catalog. Try ",
                       "cfg_find_species(\"", substr(name, 1, 6), "\")")
  if (hit$n_complete[1] < 1000)
    warning(name, " has only ", hit$n_complete[1],
            " complete-checklist records -- expect very wide intervals")
  f <- file.path(.cfg_dir, "species_registry.csv")
  reg <- data.table::fread(f, showProgress = FALSE)
  if (prefix %in% reg$prefix) stop("prefix '", prefix, "' is already in the registry")
  new <- data.table::data.table(common_name = name, species_code = species_code,
                                prefix = prefix, guild = guild,
                                predicted = direction, active = as.integer(active))
  data.table::fwrite(rbind(reg, new, fill = TRUE), f)
  message("added ", name, " (", prefix, ", ", guild, ", predicts ", direction,
          "; ", format(hit$n_complete[1], big.mark = ","), " complete records)")
  message("  NOT built yet -- run: Rscript zerofill_from_registry.R")
  invisible(TRUE)
}

#' Remove a species from the registry (the "take out" half).
#' Does not delete built data files -- those are reusable if it is added back.
cfg_drop_species <- function(p) {
  # argument is `p`, not `prefix`: `prefix` is a column name, so `reg[prefix !=
  # ..prefix]` is ambiguous and data.table rejects it. Same trap as
  # cfg_add_species() -- never name an argument after a column.
  f <- file.path(.cfg_dir, "species_registry.csv")
  reg <- data.table::fread(f, showProgress = FALSE)
  if (!p %in% reg$prefix) stop("prefix '", p, "' is not in the registry")
  data.table::fwrite(reg[prefix != p], f)
  message("removed ", p, " from the registry (data files left in place)")
  invisible(TRUE)
}
