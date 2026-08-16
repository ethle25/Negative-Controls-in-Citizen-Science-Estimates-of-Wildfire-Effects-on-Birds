# Zero-fill three focal species from the California EBD for a wildfire study.
#
# Zero-filling turns "which species were reported" into presence/absence: for
# every COMPLETE checklist in California, each focal species gets a row with
# species_observed = TRUE/FALSE and an observation_count (NA where the count was
# reported as "X"). This requires BOTH the observation file (EBD) and the
# sampling event file (the list of all checklists), filtered with the same
# criteria, and it only makes sense for complete checklists.
if (!requireNamespace("auk", quietly = TRUE)) {
  install.packages("auk", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr", repos = "https://cloud.r-project.org")
}
library(auk)
library(dplyr)

data_dir <- file.path(
  "C:/Users/ethan/OneDrive/Documents/eBird",
  "ebd_US-CA_202401_202606_smp_relJun-2026"
)
proj_dir <- "C:/Users/ethan/OneDrive/Documents/eBird"

ebd_file <- file.path(data_dir, "ebd_US-CA_202401_202606_smp_relJun-2026.txt")
smp_file <- file.path(data_dir,
                      "ebd_US-CA_202401_202606_smp_relJun-2026_sampling.txt")

# Filtered intermediates and final output.
ebd_out <- file.path(proj_dir, "ebd_focal_filtered.txt")
smp_out <- file.path(proj_dir, "ebd_focal_sampling_filtered.txt")
csv_out <- file.path(proj_dir, "zerofill_focal_species.csv")

focal_species <- c("Wrentit", "Brown Creeper", "Olive-sided Flycatcher")

# --- AWK setup (space-free gawk copy; see parse_ebd_california.R for why) ------
git_gawk <- "C:/Program Files/Git/usr/bin/gawk.exe"
awk_dir  <- "C:/Users/ethan/ebird_awk"
awk_path <- file.path(awk_dir, "gawk.exe")
awk_dlls <- c(
  "msys-2.0.dll", "msys-gmp-10.dll", "msys-mpfr-6.dll", "msys-intl-8.dll",
  "msys-readline8.dll", "msys-gcc_s-seh-1.dll", "msys-iconv-2.dll",
  "msys-ncursesw6.dll"
)
if (!file.exists(awk_path) && file.exists(git_gawk)) {
  dir.create(awk_dir, showWarnings = FALSE, recursive = TRUE)
  git_bin <- dirname(git_gawk)
  file.copy(git_gawk, awk_path, overwrite = TRUE)
  file.copy(file.path(git_bin, awk_dlls), file.path(awk_dir, awk_dlls),
            overwrite = TRUE)
}
if (file.exists(awk_path)) {
  auk_set_awk_path(awk_path, overwrite = TRUE)
}
message("Using AWK at: ", auk_get_awk_path())

# --- Filter both files: focal species, California, complete checklists only ---
filters <- auk_ebd(ebd_file, file_sampling = smp_file) %>%
  auk_species(species = focal_species) %>%
  auk_state(state = "US-CA") %>%
  auk_complete()
print(filters)

message("Filtering EBD + sampling data (streams both files, several minutes)...")
auk_filter(filters, file = ebd_out, file_sampling = smp_out, overwrite = TRUE)

# --- Zero-fill -----------------------------------------------------------------
# collapse = TRUE returns one flat data frame: one row per checklist x species.
message("Zero-filling...")
zf <- auk_zerofill(ebd_out, smp_out, collapse = TRUE)

# eBird records "present but not counted" as observation_count = "X". Convert to
# a numeric count with NA for those, and add a clear presence flag.
zf <- zf %>%
  mutate(
    observation_count = suppressWarnings(as.integer(
      ifelse(observation_count == "X", NA, observation_count))),
    present = as.integer(species_observed)
  )

# collapse_zerofill() returns scientific_name only; add the common name from the
# packaged eBird taxonomy so the table is readable.
tax <- auk::ebird_taxonomy %>%
  select(scientific_name, common_name) %>%
  distinct()
zf <- left_join(zf, tax, by = "scientific_name")

# Keep the columns useful for a fire x occurrence/migration analysis: species,
# presence, count, location, date/time, and effort covariates (needed to control
# for eBird observer effort when modelling occurrence).
keep <- c(
  "checklist_id", "scientific_name", "common_name",
  "present", "observation_count",  # present = 1/0 presence flag; count NA if "X"
  "latitude", "longitude", "locality_id",
  "observation_date", "time_observations_started",
  "protocol_type", "duration_minutes", "effort_distance_km",
  "number_observers", "all_species_reported"
)
keep <- keep[keep %in% names(zf)]
tbl <- zf[, keep]

readr::write_csv(tbl, csv_out, na = "")

# --- In-memory data variables --------------------------------------------------
# The full zero-filled presence/count table, plus one data frame per focal
# species, available in the session after running/sourcing this script.
focal_data <- tbl

wrentit                <- dplyr::filter(focal_data, common_name == "Wrentit")
brown_creeper          <- dplyr::filter(focal_data, common_name == "Brown Creeper")
olive_sided_flycatcher <- dplyr::filter(focal_data, common_name == "Olive-sided Flycatcher")

# --- Summary -------------------------------------------------------------------
n_checklists <- length(unique(tbl$checklist_id))
message("Zero-filled table written to: ", csv_out)
message("Complete checklists: ", format(n_checklists, big.mark = ","))
message("Rows (checklists x species): ", format(nrow(tbl), big.mark = ","))
message("Detections per focal species:")
print(
  tbl %>%
    group_by(common_name) %>%
    summarise(
      checklists = n(),
      detections = sum(present),
      detection_rate = round(mean(present), 4),
      .groups = "drop"
    )
)
