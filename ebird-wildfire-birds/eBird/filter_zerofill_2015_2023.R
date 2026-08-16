# =============================================================================
# Filter + zero-fill the 2015-2023 California EBD for the three focal species.
#
# Mirrors zerofill_focal_species.R (which handled the 2024-01..2026-06 download)
# but points at the newly supplied 2015-01..2023-12 download. Outputs are written
# with a _2015_2023 suffix so the existing 2024+ intermediates are untouched;
# the two zero-fill tables are concatenated later by join_all_years.R.
#
# Input : C:/Users/ethan/Downloads/ebd_US-CA_201501_202312_smp_relJun-2026/
#           ebd_US-CA_201501_202312_smp_relJun-2026.txt          (33.2 GB)
#           ebd_US-CA_201501_202312_smp_relJun-2026_sampling.txt ( 1.4 GB)
# Output: ebd_focal_filtered_2015_2023.txt
#         ebd_focal_sampling_filtered_2015_2023.txt
#         zerofill_focal_species_2015_2023.csv
# =============================================================================
suppressMessages({ library(auk); library(dplyr) })

log <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                     paste0(...)))

data_dir <- "C:/Users/ethan/Downloads/ebd_US-CA_201501_202312_smp_relJun-2026"
proj_dir <- "C:/eBird"

ebd_file <- file.path(data_dir, "ebd_US-CA_201501_202312_smp_relJun-2026.txt")
smp_file <- file.path(data_dir, "ebd_US-CA_201501_202312_smp_relJun-2026_sampling.txt")

ebd_out <- file.path(proj_dir, "ebd_focal_filtered_2015_2023.txt")
smp_out <- file.path(proj_dir, "ebd_focal_sampling_filtered_2015_2023.txt")
csv_out <- file.path(proj_dir, "zerofill_focal_species_2015_2023.csv")

focal_species <- c("Wrentit", "Brown Creeper", "Olive-sided Flycatcher")

# --- AWK setup: auk shells out to gawk, and a path containing spaces breaks the
# command it builds, so use the space-free copy under C:/Users/ethan/ebird_awk.
git_gawk <- "C:/Program Files/Git/usr/bin/gawk.exe"
awk_dir  <- "C:/Users/ethan/ebird_awk"
awk_path <- file.path(awk_dir, "gawk.exe")
awk_dlls <- c("msys-2.0.dll", "msys-gmp-10.dll", "msys-mpfr-6.dll",
              "msys-intl-8.dll", "msys-readline8.dll", "msys-gcc_s-seh-1.dll",
              "msys-iconv-2.dll", "msys-ncursesw6.dll")
if (!file.exists(awk_path) && file.exists(git_gawk)) {
  dir.create(awk_dir, showWarnings = FALSE, recursive = TRUE)
  git_bin <- dirname(git_gawk)
  file.copy(git_gawk, awk_path, overwrite = TRUE)
  file.copy(file.path(git_bin, awk_dlls), file.path(awk_dir, awk_dlls),
            overwrite = TRUE)
}
if (file.exists(awk_path)) auk_set_awk_path(awk_path, overwrite = TRUE)
log("Using AWK at: ", auk_get_awk_path())

stopifnot(file.exists(ebd_file), file.exists(smp_file))

# --- Filter: focal species, California, complete checklists only --------------
filters <- auk_ebd(ebd_file, file_sampling = smp_file) %>%
  auk_species(species = focal_species) %>%
  auk_state(state = "US-CA") %>%
  auk_complete()
print(filters)

log("Filtering 33 GB EBD + 1.4 GB sampling (streams both through gawk)...")
auk_filter(filters, file = ebd_out, file_sampling = smp_out, overwrite = TRUE)
log("Filter done.")

# --- Zero-fill ----------------------------------------------------------------
log("Zero-filling...")
zf <- auk_zerofill(ebd_out, smp_out, collapse = TRUE)

zf <- zf %>%
  mutate(
    observation_count = suppressWarnings(as.integer(
      ifelse(observation_count == "X", NA, observation_count))),
    present = as.integer(species_observed)
  )

tax <- auk::ebird_taxonomy %>% select(scientific_name, common_name) %>% distinct()
zf <- left_join(zf, tax, by = "scientific_name")

keep <- c("checklist_id", "scientific_name", "common_name",
          "present", "observation_count",
          "latitude", "longitude", "locality_id",
          "observation_date", "time_observations_started",
          "protocol_type", "duration_minutes", "effort_distance_km",
          "number_observers", "all_species_reported")
keep <- keep[keep %in% names(zf)]
tbl  <- zf[, keep]

readr::write_csv(tbl, csv_out, na = "")
log("Wrote ", csv_out)

# --- Summary ------------------------------------------------------------------
log("Complete checklists: ", format(length(unique(tbl$checklist_id)), big.mark = ","))
log("Rows (checklists x species): ", format(nrow(tbl), big.mark = ","))
print(tbl %>% group_by(common_name) %>%
        summarise(checklists = n(), detections = sum(present),
                  detection_rate = round(mean(present), 4), .groups = "drop"))
cat("\nDATE RANGE: ", format(min(tbl$observation_date)), " .. ",
    format(max(tbl$observation_date)), "\n")
log("DONE")
