# =============================================================================
# Zero-fill every ACTIVE species in species_registry.csv, in ONE pass
#
# WHY THIS EXISTS
#   zerofill_focal_species.R hardcodes three species and points at Windows paths
#   (C:/Users/ethan/...), so it has never run on this machine. That is why the
#   pipeline is species-switchable downstream but stuck at three species
#   upstream: the join can filter to any species you like, but the zero-filled
#   input only ever contained Wrentit, Brown Creeper and Olive-sided Flycatcher.
#
#   This replaces it. Species come from the registry, paths are macOS, and the
#   expensive part -- streaming 52 GB of raw EBD -- happens ONCE for all species
#   rather than once per species.
#
# COST
#   The filter pass is the bottleneck: ~52 GB across two downloads. Measured awk
#   throughput on this machine is ~110 MB/s for a trivial scan, so budget 30-60
#   min per span with species filtering on top. Zero-filling and writing is
#   another ~1 h for a dozen species. Run it once, overnight, then every species
#   in the registry is available to the rest of the pipeline for ~20 min each.
#
# MEMORY
#   auk_zerofill(collapse = TRUE) builds one flat frame of
#   n_checklists x n_species rows. At 19 species that is tens of millions of
#   rows and will not fit comfortably in 24 GB. This uses collapse = FALSE,
#   which returns observations and sampling events separately, then joins and
#   writes ONE FILE PER SPECIES. Downstream that is also faster: the join reads
#   only its own species instead of filtering a multi-GB combined file.
#
# OUTPUT (new filenames -- the existing 3-species files are NOT touched)
#   zerofill_<prefix>_2015_2023.csv
#   zerofill_<prefix>_2024_2026.csv
#   zerofill_registry_manifest.csv    what was written, with row counts
#
# INTEGRATION
#   join_species_all_years.R currently hardcodes the two old zerofill paths at
#   line ~69. Point it at the per-species files (or set EBIRD_ZF_PREFIX) once
#   this has run. Nothing else downstream changes.
# =============================================================================
suppressMessages({ library(data.table); library(auk) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))

DRY <- Sys.getenv("ZF_DRY_RUN", "0") != "0"   # ZF_DRY_RUN=1 -> plan only, no work

# --- AWK check ---------------------------------------------------------------
# auk streams the EBD with awk. macOS ships BSD awk, which auk tolerates for
# simple filters but is markedly slower and has bitten people on large files.
# gawk is the documented recommendation. Fail loudly NOW rather than three
# hours into a 52 GB pass.
awk_path <- auk::auk_get_awk_path()
if (is.na(awk_path) || !nzchar(awk_path)) stop("no awk found; install gawk: brew install gawk")
is_gawk <- grepl("gawk", awk_path) ||
           any(grepl("GNU Awk", tryCatch(system2(awk_path, "--version", stdout = TRUE,
                                                 stderr = FALSE)[1], error = function(e) "")))
if (!is_gawk) {
  say("WARNING: awk at ", awk_path, " is not GNU awk.")
  say("  auk recommends gawk for files this size. If the filter pass is")
  say("  unbearably slow or errors out:  brew install gawk")
  say("  then re-run; auk picks gawk up automatically.")
}

# --- species from the registry ----------------------------------------------
sp <- cfg_species(active_only = TRUE)
say("active species in registry: ", nrow(sp))
print(as.data.frame(sp[, .(common_name, species_code, prefix, guild, predicted)]))

# --- the two raw downloads ---------------------------------------------------
spans <- list(
  list(tag = "2015_2023",
       dir = file.path(dir_, "ebd_US-CA_201501_202312_smp_relJun-2026"),
       ebd = "ebd_US-CA_201501_202312_smp_relJun-2026.txt",
       smp = "ebd_US-CA_201501_202312_smp_relJun-2026_sampling.txt"),
  list(tag = "2024_2026",
       dir = file.path(dir_, "ebd_US-CA_202401_202606_smp_relJun-2026"),
       ebd = "ebd_US-CA_202401_202606_smp_relJun-2026.txt",
       smp = "ebd_US-CA_202401_202606_smp_relJun-2026_sampling.txt"))

for (s in spans) {
  for (f in c(s$ebd, s$smp)) {
    p <- file.path(s$dir, f)
    if (!file.exists(p)) stop("missing raw EBD file: ", p)
  }
}
say("raw EBD present for both spans (",
    sprintf("%.1f GB", sum(sapply(spans, function(s)
      file.size(file.path(s$dir, s$ebd)))) / 1e9), " of observations)")

if (DRY) { say("ZF_DRY_RUN=1 -- stopping before any work"); quit(save = "no") }

# --- per span: one filter pass, then one zero-fill, then split by species ----
manifest <- list()

for (s in spans) {
  say("=== span ", s$tag, " ===")
  ebd_in  <- file.path(s$dir, s$ebd)
  smp_in  <- file.path(s$dir, s$smp)
  ebd_out <- file.path(dir_, sprintf("ebd_registry_filtered_%s.txt", s$tag))
  smp_out <- file.path(dir_, sprintf("ebd_registry_sampling_%s.txt", s$tag))

  # The filtered file is specific to the species set it was built for. Reusing it
  # after a species is ADDED would silently produce a file with no rows for the
  # new bird -- and line ~141 would just say "skipping". So record the set in a
  # sidecar and re-filter whenever it changes. Same class of bug as the cell_id
  # cache fingerprint in regrid_cells.R.
  stamp <- file.path(dir_, sprintf(".zf_species_%s.txt", s$tag))
  want  <- paste(sort(sp$common_name), collapse = "|")
  have  <- if (file.exists(stamp)) readLines(stamp, warn = FALSE)[1] else ""

  if (file.exists(ebd_out) && file.exists(smp_out) && identical(have, want)) {
    say("  filtered files already exist for this exact species set, reusing")
  } else {
    if (file.exists(ebd_out) && !identical(have, want))
      say("  species set CHANGED since the last filter -- re-filtering")
    say("  filtering: ", nrow(sp), " species, US-CA, complete checklists only")
    say("  this is the 52 GB pass -- expect 30-60 min")
    filters <- auk_ebd(ebd_in, file_sampling = smp_in)
    filters <- auk_species(filters, species = sp$common_name)
    filters <- auk_state(filters, state = "US-CA")
    filters <- auk_complete(filters)
    auk_filter(filters, file = ebd_out, file_sampling = smp_out, overwrite = TRUE)
    writeLines(want, stamp)
    say("  filtered -> ", sprintf("%.0f MB obs, %.0f MB sampling",
        file.size(ebd_out)/1e6, file.size(smp_out)/1e6))
  }

  # collapse = FALSE keeps observations and sampling events apart, so peak
  # memory is one copy of each rather than their cross product
  say("  zero-filling")
  zf  <- auk_zerofill(ebd_out, smp_out, collapse = FALSE)
  obs <- as.data.table(zf$observations)
  smp <- as.data.table(zf$sampling_events)
  say("  checklists ", fmt(nrow(smp)), " | observation rows ", fmt(nrow(obs)))

  keep_smp <- intersect(c("checklist_id","latitude","longitude","locality_id",
                          "observation_date","time_observations_started",
                          "protocol_type","duration_minutes","effort_distance_km",
                          "number_observers","all_species_reported"), names(smp))
  smp <- smp[, ..keep_smp]

  # auk_zerofill()'s observations table identifies the species by SCIENTIFIC
  # NAME, not species_code -- its columns are checklist_id, scientific_name,
  # breeding_*, behavior_code, age_sex, observation_count, species_observed.
  # Matching on species_code silently found nothing and every species was
  # "skipped". Resolve each registry row to its scientific name via the catalog
  # (built from the actual download) and fall back to auk's taxonomy.
  tax <- as.data.table(auk::ebird_taxonomy)[, .(species_code, common_name, scientific_name)]
  cat_ <- tryCatch(cfg_catalog()[, .(common_name, scientific_name)],
                   error = function(e) tax[, .(common_name, scientific_name)])

  for (i in seq_len(nrow(sp))) {
    code <- sp$species_code[i]; pfx <- sp$prefix[i]; cn <- sp$common_name[i]
    sci <- cat_[common_name == cn]$scientific_name
    if (!length(sci)) sci <- tax[common_name == cn]$scientific_name
    if (!length(sci)) { say("  !! no scientific name for ", cn, " -- skipping"); next }
    o <- obs[scientific_name == sci[1]]
    if (!nrow(o)) { say("  !! no rows for ", cn, " (", sci[1], ") -- skipping"); next }
    o <- o[, .(checklist_id, scientific_name, observation_count, species_observed)]
    # eBird writes "present but not counted" as X -> NA count, presence flag stays 1
    o[, `:=`(observation_count = suppressWarnings(as.integer(
                fifelse(observation_count == "X", NA_character_, observation_count))),
             present = as.integer(species_observed))]
    o[, species_observed := NULL]
    out <- merge(smp, o, by = "checklist_id", all.x = FALSE)
    out[, common_name := cn]
    setcolorder(out, c("checklist_id","scientific_name","common_name",
                       "present","observation_count"))
    f <- file.path(dir_, sprintf("zerofill_%s_%s.csv", pfx, s$tag))
    fwrite(out, f, na = "")
    manifest[[length(manifest)+1]] <- data.table(
      span = s$tag, common_name = cn, species_code = code, prefix = pfx,
      rows = nrow(out), detections = sum(out$present, na.rm = TRUE),
      rate = round(mean(out$present, na.rm = TRUE), 5),
      file = basename(f), mb = round(file.size(f)/1e6, 1))
    say(sprintf("  %-26s %10s rows  %8s det  rate %.4f",
                cn, fmt(nrow(out)), fmt(sum(out$present, na.rm=TRUE)),
                mean(out$present, na.rm = TRUE)))
  }
  rm(obs, smp, zf); invisible(gc())
}

mf <- rbindlist(manifest)
fwrite(mf, file.path(dir_, "zerofill_registry_manifest.csv"))
say("wrote zerofill_registry_manifest.csv")
print(as.data.frame(mf))

say("DONE")
say("NEXT: point join_species_all_years.R (line ~69) at zerofill_<prefix>_*.csv,")
say("      then ./run_species.sh runs every active species in the registry.")
