# Parse the eBird Basic Dataset (EBD) with auk and filter to California.
# The EBD is a huge tab-delimited text file; auk streams it through AWK so it
# never has to be loaded fully into memory.
if (!requireNamespace("auk", quietly = TRUE)) {
  install.packages("auk", repos = "https://cloud.r-project.org")
}
library(auk)

data_dir <- file.path(
  "C:/Users/ethan/OneDrive/Documents/eBird",
  "ebd_US-CA_202401_202606_smp_relJun-2026"
)

ebd_file <- file.path(data_dir, "ebd_US-CA_202401_202606_smp_relJun-2026.txt")
out_file <- file.path(
  "C:/Users/ethan/OneDrive/Documents/eBird",
  "ebd_california_filtered.txt"
)

# auk needs an AWK executable. On Windows we use the gawk that ships with Git
# for Windows. Two Windows-specific gotchas:
#   1. auk builds its command as system(paste(path, ...)) WITHOUT quoting, so a
#      path containing spaces ("C:/Program Files/...") is split and fails.
#   2. auk also runs normalizePath() on the path, which expands an 8.3 short
#      path back to the spaced long form, so the short-path trick doesn't help.
# The reliable fix is a genuinely space-free copy of gawk.exe alongside the msys
# DLLs it depends on. Provision that copy once, then point auk at it.
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

# Build the filter: state = California (US-CA).
ebd <- auk_ebd(ebd_file)
filters <- auk_state(ebd, state = "US-CA")
print(filters)

# Run the filter. This streams the ~16 GB file through AWK and writes only the
# matching rows to out_file. overwrite = TRUE lets the script be re-run.
message("Filtering... this reads the full file and may take several minutes.")
auk_filter(filters, file = out_file, overwrite = TRUE)

message("Done. Filtered file written to: ", out_file)
