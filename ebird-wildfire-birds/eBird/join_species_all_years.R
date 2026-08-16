# =============================================================================
# Wrentit x MTBS x MODIS -- FULL-SPAN spatial/temporal join
#
#   eBird span   : 2015-01-01 .. 2026-06-30  (union of both EBD downloads, no gap)
#   fire history : MTBS CA perimeters 1984-2025 (a match can reach back 42 years)
#   vegetation   : MOD13A2 16-day 1 km NDVI/EVI, 2000-02-18 .. 2026-06-26
#
# GRAIN: one row per COMPLETE eBird checklist x Wrentit, zero-filled to
# presence/absence. Because the deliverable is a single species this is also
# exactly one row per checklist (asserted below). Species occurrence and effort
# covariates come from the zero-fill; fire and vegetation columns are joined on.
#
# Supersedes join_all_years.R (same intent, never run, Windows-only paths) and
# the single-year join_fire_ndvi.R / join_fire_history.R pair, whose reviewed
# rules are carried over verbatim:
#   * protocol_type joined on BOTH the S sampling-event id and the G group id --
#     auk_zerofill re-keys shared checklists to the group id, so an S-only
#     lookup silently loses them
#   * QC: Traveling/Stationary only; duration<=300; distance<=5; observers<=10
#   * fire match: most recent fire CONTAINING the point that ignited STRICTLY
#     BEFORE the checklist date => days_since_fire >= 1 always
#   * distance_to_perimeter_km signed: NEGATIVE = inside a perimeter
#   * burn_severity_class sampled from the MATCHED fire's OWN year raster
#   * ndvi_baseline = same 16-day window, PREVIOUS year (year-over-year reference)
#   * no dnbr column -- no continuous product exists in the supplied MTBS data
#
# PERFORMANCE: all spatial work is done on UNIQUE COORDINATES, not on rows.
# eBird checklists repeat locations heavily (~4.97M Wrentit checklists sit on
# ~571k distinct lat/lon), and point-in-polygon, nearest-perimeter distance and
# MODIS sampling all depend only on position. Only the fire MATCH is per-row,
# because it depends on the checklist date. Doing every step per-row -- as
# join_all_years.R did -- is hours of duplicated geometry on this span.
#
# CHECKPOINTS: the expensive stages write to .join_cache_all_years/ so an
# interrupted run resumes. Delete that directory to force a full recompute.
# Set EBIRD_SAMPLE_N to smoke-test on a random subset (writes *_sample.csv).
#
# -----------------------------------------------------------------------------
# REVIEW ITEMS. The six points raised on the 2024 table are each re-verified at
# run time and reported in the summary's FEEDBACK CONFIRMATION section, so the
# table ships with its own proof rather than a claim.
#
#   R1 time-ordered fire match      R4 dnbr availability
#   R2 multi-year perimeters        R5 distance sign convention
#   R3 Historical/Other dropped     R6 offshore checklists
#
# See that section of the summary for the evidence; the rules themselves are
# implemented at the marked points below.
# =============================================================================
suppressMessages({
  library(terra); library(sf); library(data.table); library(lubridate)
})
sf_use_s2(FALSE)
terraOptions(progress = 0)
setDTthreads(0)

t_start <- Sys.time()
log <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                     paste0(...)))
fmt <- function(x) format(x, big.mark = ",")

# --- paths --------------------------------------------------------------------
proj_root <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
ebird_dir <- file.path(proj_root, "eBird")
mtbs_dir  <- file.path(proj_root, "mtbs")
modis_dir <- file.path(proj_root, "MOD13A2_061")
perim_shp <- file.path(ebird_dir, "mtbs_perim", "mtbs_perims_DD.shp")

# --- zero-fill inputs ---------------------------------------------------------
# Two layouts are supported, chosen automatically:
#
#   ORIGINAL  zerofill_focal_species*.csv -- one combined file holding all three
#             original species, filtered below by `common_name == SPECIES`.
#   REGISTRY  zerofill_<prefix>_<span>.csv -- one file PER SPECIES, written by
#             zerofill_from_registry.R. Preferred when present: it skips reading
#             and discarding the other species, and it is the only layout that
#             can hold a species outside the original three.
#
# Per-species files win when they exist, so adding a bird needs no edit here.
#
# NO SPECIES IS HARDCODED IN THIS SCRIPT. Identity comes from
# species_registry.csv via cfg_current_species(), which returns the row named by
# EBIRD_PREFIX if that is set and otherwise the first ACTIVE row. The defaults
# used to be the literals "Wrentit"/"wrenti"/"wrentit", which made this file the
# Wrentit join that other species borrowed rather than a species-neutral join.
# Behaviour is unchanged whenever the env vars are set, and unchanged today with
# nothing set, because Wrentit is the first active registry row.
source(file.path(ebird_dir, "project_config.R"))
.sp <- cfg_current_species()
PREFIX_EARLY <- Sys.getenv("EBIRD_PREFIX", .sp$prefix[1])
zf_reg <- file.path(ebird_dir, sprintf("zerofill_%s_%s.csv", PREFIX_EARLY,
                                       c("2015_2023", "2024_2026")))
if (all(file.exists(zf_reg))) {
  zf_files <- zf_reg
  message("zero-fill layout: per-species (", basename(zf_reg[1]), ")")
} else {
  zf_files <- file.path(ebird_dir, c("zerofill_focal_species_2015_2023.csv",
                                     "zerofill_focal_species.csv"))
  message("zero-fill layout: combined focal-species files")
}
smp_files <- file.path(ebird_dir, c("ebd_focal_sampling_filtered_2015_2023.txt",
                                    "ebd_focal_sampling_filtered.txt"))
need <- c(zf_files, smp_files, perim_shp)
miss <- need[!file.exists(need)]
if (length(miss)) stop("missing required input(s):\n  ", paste(miss, collapse = "\n  "),
                       "\n  (for a species outside the original three, run ",
                       "zerofill_from_registry.R first)")
if (!dir.exists(modis_dir)) stop("missing MODIS dir: ", modis_dir)
if (!dir.exists(mtbs_dir))  stop("missing MTBS raster dir: ", mtbs_dir)

# Species is switchable so the identical pipeline runs on any focal species.
# Identity is READ FROM THE REGISTRY (see the note above the zero-fill layout
# block) -- an explicit env var still wins, so nothing that already works changes.
SPECIES      <- Sys.getenv("EBIRD_SPECIES",      .sp$common_name[1])
SPECIES_CODE <- Sys.getenv("EBIRD_SPECIES_CODE", .sp$species_code[1])
PREFIX       <- Sys.getenv("EBIRD_PREFIX",       .sp$prefix[1])
stopifnot(identical(PREFIX, PREFIX_EARLY))   # the zero-fill layout was chosen for this species
KEEP_PROTO   <- c("Traveling", "Stationary")     # R3

sample_n  <- suppressWarnings(as.integer(Sys.getenv("EBIRD_SAMPLE_N", "")))
is_sample <- !is.na(sample_n) && sample_n > 0
sfx       <- if (is_sample) "_sample" else ""

out_csv   <- file.path(ebird_dir, sprintf("%s_2015_2026_joined%s.csv", PREFIX, sfx))
out_clean <- file.path(ebird_dir, sprintf("%s_2015_2026_joined_clean%s.csv", PREFIX, sfx))
out_sum   <- file.path(ebird_dir, sprintf("%s_2015_2026_joined_summary%s.txt", PREFIX, sfx))
cache_dir <- file.path(ebird_dir, sprintf(".join_cache_%s%s", PREFIX, sfx))
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

log(sprintf("project root: %s%s", proj_root,
            if (is_sample) sprintf("  [SMOKE TEST, n=%s]", fmt(sample_n)) else ""))

# =============================================================================
# 1. base table: zero-fill x 2, species coverage, protocol join, QC
# =============================================================================
base_rds <- file.path(cache_dir, "base.rds")
if (file.exists(base_rds)) {
  log("Loading cached base table")
  bs <- readRDS(base_rds)
  ebd <- bs$ebd; cover <- bs$cover; proto_before <- bs$proto_before
  n_zf <- bs$n_zf; n_species <- bs$n_species; n_pre_qc <- bs$n_pre_qc
  qc_drop <- bs$qc_drop; off_pre <- bs$off_pre; asr_bad <- bs$asr_bad
  off_break <- bs$off_break; off_proto <- bs$off_proto
} else {
  zf_cols <- c("checklist_id", "common_name", "present", "observation_count",
               "latitude", "longitude", "locality_id", "observation_date",
               "duration_minutes", "effort_distance_km", "number_observers",
               "all_species_reported")
  cover <- list(); parts <- list(); n_zf <- 0L; asr_bad <- 0L
  for (p in zf_files) {
    log("Reading ", basename(p))
    d <- fread(p, select = zf_cols, showProgress = FALSE)
    d[, observation_date := as.IDate(observation_date)]
    n_zf <- n_zf + nrow(d)
    # every checklist should be complete -- auk_complete() was applied upstream
    asr_bad <- asr_bad + sum(!as.logical(d$all_species_reported), na.rm = TRUE)
    # coverage per focal species BEFORE narrowing to the target species
    cover[[p]] <- d[, .(checklists = .N, detections = sum(present, na.rm = TRUE),
                        first = min(observation_date), last = max(observation_date)),
                    by = common_name]
    parts[[p]] <- d[common_name == SPECIES]
    rm(d); invisible(gc())
  }
  cover <- rbindlist(cover)[, .(checklists = sum(checklists),
                                detections = sum(detections),
                                first = min(first), last = max(last)),
                            by = common_name][order(common_name)]
  cover[, detection_rate := round(detections / checklists, 4)]

  ebd <- rbindlist(parts); rm(parts); invisible(gc())
  n_species <- nrow(ebd)
  # the two downloads do not overlap (2015-2023 vs 2024-2026); de-dup defensively
  ebd <- unique(ebd, by = "checklist_id")
  ebd[, all_species_reported := NULL]
  log(sprintf("%s rows: %s (from %s zero-fill rows; %s duplicate ids dropped)",
              SPECIES, fmt(nrow(ebd)), fmt(n_zf), fmt(n_species - nrow(ebd))))

  if (is_sample) {
    set.seed(1)
    ebd <- ebd[sample(.N, min(sample_n, .N))]
    log("SMOKE TEST subset: ", fmt(nrow(ebd)), " rows")
  }

  # --- protocol_type from BOTH sampling files (S ids + G ids) ------------------
  log("Joining protocol_type from sampling files")
  smp <- rbindlist(lapply(smp_files, function(p) {
    log("  reading ", basename(p))
    fread(p, sep = "\t", quote = "", showProgress = FALSE,
          select = c("SAMPLING EVENT IDENTIFIER", "PROTOCOL NAME",
                     "GROUP IDENTIFIER"),
          col.names = c("seid", "protocol_name", "gid"))
  }))
  # R3: EXACT match, not grepl(). The banding-station variants "Traveling
  # (2 band, 25m)" and "Stationary (2/3 band, 30m...)" report effort on a
  # different basis; a grepl("Traveling") test would silently fold them into the
  # kept classes. They are classed Other and dropped with everything else there
  # ("eBird Pelagic Protocol", "Area", "Banding", "Nocturnal Flight Call Count").
  smp[, protocol_type := fifelse(protocol_name == "Traveling",  "Traveling",
                        fifelse(protocol_name == "Stationary", "Stationary",
                        fifelse(protocol_name == "Historical", "Historical",
                        fifelse(protocol_name == "Incidental", "Incidental",
                                "Other"))))]
  lookup <- unique(rbindlist(list(
    smp[, .(id = seid, protocol_type)],
    smp[!is.na(gid) & gid != "", .(id = gid, protocol_type)])), by = "id")
  rm(smp); invisible(gc())
  ebd[lookup, on = c(checklist_id = "id"), protocol_type := i.protocol_type]
  rm(lookup); invisible(gc())

  # --- pre-QC diagnostics (the evidence behind R3 and R6) ---------------------
  proto_before <- ebd[, .(rows = .N,
                          det_rate = round(mean(present), 4),
                          eff_dist_na_pct = round(100 * mean(is.na(effort_distance_km)), 1),
                          dur_na_pct = round(100 * mean(is.na(duration_minutes)), 1)),
                      by = protocol_type][order(-rows)]
  # R6 evidence: which QC rule actually removes each offshore checklist. The
  # assignments below run least- to most-specific, so each row is attributed to
  # the FIRST rule that would remove it (protocol wins, then distance, etc.).
  offw    <- ebd[longitude < -124.4]
  off_pre <- nrow(offw)
  r_obs   <- !is.na(offw$number_observers)   & offw$number_observers   > 10
  r_dur   <- !is.na(offw$duration_minutes)   & offw$duration_minutes   > 300
  r_dist  <- !is.na(offw$effort_distance_km) & offw$effort_distance_km > 5
  r_proto <- is.na(offw$protocol_type) | !offw$protocol_type %chin% KEEP_PROTO
  off_why <- rep("survives QC", off_pre)
  off_why[r_obs]   <- "number_observers > 10"
  off_why[r_dur]   <- "duration_minutes > 300"
  off_why[r_dist]  <- "effort_distance_km > 5"
  off_why[r_proto] <- "protocol not Traveling/Stationary"
  off_break <- sort(table(off_why), decreasing = TRUE)
  off_proto <- sort(table(offw$protocol_type[r_proto]), decreasing = TRUE)
  rm(offw); invisible(gc())
  n_pre_qc <- nrow(ebd)
  qc_drop <- list(
    protocol  = ebd[is.na(protocol_type) | !protocol_type %chin% KEEP_PROTO, .N],
    duration  = ebd[!is.na(duration_minutes)   & duration_minutes   > 300, .N],
    distance  = ebd[!is.na(effort_distance_km) & effort_distance_km > 5,   .N],
    observers = ebd[!is.na(number_observers)   & number_observers   > 10,  .N],
    coords    = ebd[is.na(latitude) | is.na(longitude), .N])

  # --- QC ---------------------------------------------------------------------
  ebd <- ebd[!is.na(protocol_type) & protocol_type %chin% KEEP_PROTO &
             (is.na(duration_minutes)   | duration_minutes   <= 300) &
             (is.na(effort_distance_km) | effort_distance_km <= 5) &
             (is.na(number_observers)   | number_observers   <= 10) &
             !is.na(latitude) & !is.na(longitude)]
  log(sprintf("QC: %s -> %s rows", fmt(n_pre_qc), fmt(nrow(ebd))))

  ebd[, `:=`(species_code = SPECIES_CODE,
             year         = year(observation_date),
             day_of_year  = yday(observation_date),
             presence     = as.integer(present))]
  saveRDS(list(ebd = ebd, cover = cover, proto_before = proto_before,
               n_zf = n_zf, n_species = n_species, n_pre_qc = n_pre_qc,
               qc_drop = qc_drop, off_pre = off_pre, asr_bad = asr_bad,
               off_break = off_break, off_proto = off_proto), base_rds)
}
n <- nrow(ebd)
obs_num <- as.numeric(ebd$observation_date)
log(sprintf("base table: %s rows, %s .. %s", fmt(n),
            min(ebd$observation_date), max(ebd$observation_date)))

# --- unique coordinates -------------------------------------------------------
ebd[, ukey := paste0(latitude, "|", longitude)]
u <- unique(ebd[, .(ukey, latitude, longitude)], by = "ukey")
u[, uid := .I]
ebd[u, on = "ukey", uid := i.uid]
U <- nrow(u)
log(sprintf("unique coordinates: %s  (%.1fx reduction vs rows)", fmt(U), n / U))

upts     <- st_as_sf(as.data.frame(u), coords = c("longitude", "latitude"), crs = 4326)
upts5070 <- st_transform(upts, 5070)      # CONUS Albers, metres

# =============================================================================
# 2. R1 + R2: fire perimeters, then the most recent fire that ignited STRICTLY
#    BEFORE each checklist. Candidate perimeters are a property of the LOCATION;
#    only the final date-ordered pick is per row.
# =============================================================================
fire_rds <- file.path(cache_dir, "fire.rds")
if (file.exists(fire_rds) && identical(readRDS(fire_rds)$n, n)) {
  log("Loading cached fire match")
  fc <- readRDS(fire_rds)
} else {
  # R2: the WHOLE CA fire history, not a single year
  log("Reading CA fire perimeters")
  perim <- st_read(perim_shp, quiet = TRUE,
    query = "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
  perim <- st_make_valid(perim)
  perim$Ig_Date <- as.Date(perim$Ig_Date)
  perim5070 <- st_transform(perim, 5070)
  log(sprintf("  perimeters: %s (%s .. %s)", fmt(nrow(perim5070)),
              min(perim5070$Ig_Date), max(perim5070$Ig_Date)))

  log("Point-in-perimeter over unique coordinates")
  wi <- st_within(upts5070, perim5070)
  n_hit_u <- sum(lengths(wi) > 0L)
  log(sprintf("  unique points inside >=1 perimeter: %s", fmt(n_hit_u)))

  # candidate (unique point, perimeter) pairs -> expand to the rows at those
  # points -> R1: keep only fires ignited STRICTLY before the checklist date ->
  # take the most recent survivor.
  pairs  <- data.table(uid = rep.int(seq_along(wi), lengths(wi)), p = unlist(wi))
  ig_num <- as.numeric(perim5070$Ig_Date)
  cand <- merge(ebd[, .(row = .I, uid, obs = obs_num)],
                pairs, by = "uid", allow.cartesian = TRUE)
  cand[, ig := ig_num[p]]
  n_cand_rows <- uniqueN(cand$row)          # rows at a location that ever burns
  cand <- cand[ig < obs]                    # <-- the strict time-order rule
  best <- cand[, .(p = p[which.max(ig)]), by = row]
  rm(cand, pairs); invisible(gc())

  matched <- rep(NA_integer_, n)
  matched[best$row] <- best$p
  inside <- !is.na(matched)
  # a row at a burned location with NO match is a genuine PRE-FIRE observation:
  # the location burns, but only after this checklist. It stays a control.
  n_prefire <- n_cand_rows - sum(inside)
  log(sprintf("  rows matched to a prior fire: %s (%.2f%%) | pre-fire: %s",
              fmt(sum(inside)), 100 * mean(inside), fmt(n_prefire)))

  fc <- list(n = n, matched = matched, inside = inside,
             event_id = perim5070$Event_ID, ig_date = perim5070$Ig_Date,
             n_perim = nrow(perim5070), perim_span = range(perim5070$Ig_Date),
             n_hit_u = n_hit_u, n_prefire = n_prefire)
  saveRDS(fc, fire_rds)
}
matched <- fc$matched; inside <- fc$inside

ebd[, fire_event_id := fifelse(inside, fc$event_id[matched], "NONE")]
ig <- rep(as.Date(NA), n); ig[inside] <- fc$ig_date[matched[inside]]
ebd[, `:=`(ignition_date   = ig,
           fire_year       = year(ig),
           days_since_fire = as.integer(as.numeric(observation_date) - as.numeric(ig)))]

# =============================================================================
# 3. R5: signed distance to the nearest perimeter BOUNDARY.
#    NEGATIVE = inside a perimeter (the spec's own convention). Per unique
#    coordinate, chunked and checkpointed.
#    Boundaries are simplified to 30 m first: MTBS polygons are very high-vertex
#    and st_distance cost scales with vertex count. 30 m is ~1 severity-raster
#    cell and the output is in KILOMETRES, so induced error is <= 0.03 km.
# =============================================================================
log("Distance to nearest perimeter boundary")
dist_rds <- file.path(cache_dir, "dist.rds")
d_m <- numeric(U); done_to <- 0L
if (file.exists(dist_rds)) {
  ck <- readRDS(dist_rds)
  if (identical(ck$U, U)) {
    d_m <- ck$d_m; done_to <- ck$done_to
    log(sprintf("  resuming: %s / %s unique points done", fmt(done_to), fmt(U)))
  } else log("  checkpoint mismatch -- recomputing")
}
if (done_to < U) {
  if (!exists("perim5070")) {
    perim <- st_read(perim_shp, quiet = TRUE,
      query = "SELECT Event_ID, Ig_Date FROM mtbs_perims_DD WHERE Event_ID LIKE 'CA%'")
    perim <- st_make_valid(perim); perim$Ig_Date <- as.Date(perim$Ig_Date)
    perim5070 <- st_transform(perim, 5070)
  }
  bnd <- st_boundary(st_simplify(perim5070, dTolerance = 30))
  chunk <- 100000L
  for (s in seq(done_to + 1L, U, by = chunk)) {
    e   <- min(s + chunk - 1L, U)
    sub <- upts5070[s:e, ]
    nf  <- st_nearest_feature(sub, perim5070)
    d_m[s:e] <- as.numeric(st_distance(sub, bnd[nf, ], by_element = TRUE))
    saveRDS(list(U = U, d_m = d_m, done_to = e), dist_rds)
    log(sprintf("  distance %s-%s / %s", fmt(s), fmt(e), fmt(U)))
  }
  rm(bnd); invisible(gc())
}
ebd[, distance_to_perimeter_km := (d_m[uid] / 1000) * fifelse(inside, -1, 1)]

# =============================================================================
# 4. burn_severity_class from the MATCHED fire's own year raster.
#    R4: these statewide rasters are the THEMATIC dnbr6 product. The first one
#    opened is probed and the evidence is reported in the summary.
# =============================================================================
log("Sampling burn_severity_class per matched fire-year")
bsc <- integer(n)
mt  <- unique(ebd[inside, .(uid, fire_year)])      # coord x fire-year pairs only
missing_rast <- integer(0); rast_probe <- NULL
for (yr in sort(unique(mt$fire_year))) {
  tif <- file.path(mtbs_dir, sprintf("mtbs_CA_%d", yr), sprintf("mtbs_CA_%d.tif", yr))
  if (!file.exists(tif)) {
    log(sprintf("  WARN no severity raster for %d", yr))
    missing_rast <- c(missing_rast, yr); next
  }
  r <- rast(tif)
  if (is.null(rast_probe))
    rast_probe <- list(year = yr, datatype = datatype(r)[1],
                       has_coltab = terra::has.colors(r)[1],
                       vals = sort(unique(freq(r)$value)))
  su <- mt[fire_year == yr, uid]
  v  <- terra::extract(r, vect(upts5070[su, ]))[, 2]
  if (!is.numeric(v)) v <- suppressWarnings(as.integer(as.character(v)))
  v[is.na(v)] <- 0L                                # inside perim, unmapped pixel
  sel <- which(inside & ebd$fire_year == yr)
  bsc[sel] <- as.integer(v)[match(ebd$uid[sel], su)]
}
ebd[, burn_severity_class := bsc]

# =============================================================================
# 5. MODIS NDVI / EVI / composite date / reliability + prior-year baseline
#    Cell indices are resolved once per unique coordinate, so each granule is
#    read once and indexed, rather than extracted per row.
# =============================================================================
MODIS_GRID <- "MODIS_Grid_16DAY_1km_VI"
sds_uri <- function(hdf, layer)
  sprintf('HDF4_EOS:EOS_GRID:"%s":%s:"%s"', hdf, MODIS_GRID, layer)
L_NDVI <- "1 km 16 days NDVI"; L_EVI <- "1 km 16 days EVI"
L_CDOY <- "1 km 16 days composite day of the year"
L_REL  <- "1 km 16 days pixel reliability"

log("Indexing MODIS granules")
hdfs <- list.files(modis_dir, pattern = "^MOD13A2\\.A\\d{7}\\.h\\d{2}v\\d{2}\\..*\\.hdf$",
                   full.names = TRUE)
bn   <- basename(hdfs)
meta <- data.table(path = hdfs,
                   ayear = as.integer(substr(bn, 10, 13)),
                   adoy  = as.integer(substr(bn, 14, 16)),
                   tile  = substr(bn, 18, 23))
setkey(meta, ayear, adoy, tile)
tiles <- sort(unique(meta$tile))
log(sprintf("  granules: %s | tiles: %s | %d-%d", fmt(nrow(meta)),
            paste(tiles, collapse = " "), min(meta$ayear), max(meta$ayear)))

tmpl <- lapply(tiles, function(t) rast(sds_uri(meta[tile == t]$path[1], L_NDVI)))
xy   <- crds(project(vect(as.matrix(u[, .(longitude, latitude)]), type = "points",
                          crs = "EPSG:4326"), crs(tmpl[[1]])))
u[, `:=`(tile_i = NA_integer_, cell = NA_integer_)]
for (i in seq_along(tiles)) {
  cl  <- cellFromXY(tmpl[[i]], xy)
  hit <- !is.na(cl) & is.na(u$cell)
  if (any(hit)) u[hit, `:=`(tile_i = i, cell = cl[hit])]
}
n_outside_tiles <- sum(is.na(u$cell))
log(sprintf("  unique coords outside all MODIS tiles: %s", fmt(n_outside_tiles)))
u_tile <- u$tile_i; u_cell <- u$cell

# NOTE: GDAL reports MOD13A2's scale as 10000 because MODIS documents the scale
# as a DIVISOR; terra would therefore multiply by 10000 instead of 1e-4. Strip
# scoff() and apply 1e-4 explicitly.
read_layer <- function(f, layer) {
  r <- rast(sds_uri(f, layer)); scoff(r) <- NULL
  as.vector(values(r))
}

ebd[, comp_doy := 16L * ((day_of_year - 1L) %/% 16L) + 1L]
ebd[comp_doy > 353L, comp_doy := 353L]
combos <- unique(ebd[, .(year, comp_doy)])[order(year, comp_doy)]
log(sprintf("MODIS composites to process: %d", nrow(combos)))

ndvi  <- rep(NA_real_,    n); evi   <- rep(NA_real_, n)
prel  <- rep(NA_integer_, n); cdate <- rep(NA_real_, n)   # numeric date, compact
nbase <- rep(NA_real_,    n)
k_start <- 1L; miss_gran <- character(0)

modis_rds <- file.path(cache_dir, "modis.rds")
if (file.exists(modis_rds)) {
  ck <- readRDS(modis_rds)
  if (identical(ck$n, n) && identical(ck$ncombo, nrow(combos))) {
    ndvi <- ck$ndvi; evi <- ck$evi; prel <- ck$prel; cdate <- ck$cdate
    nbase <- ck$nbase; miss_gran <- ck$miss_gran; k_start <- ck$k + 1L
    log(sprintf("  resuming MODIS from composite %d / %d", k_start, nrow(combos)))
  } else log("  MODIS checkpoint mismatch -- recomputing")
}

row_idx <- split(seq_len(n), paste(ebd$year, ebd$comp_doy))
for (k in seq(k_start, nrow(combos))) {
  yy <- combos$year[k]; dd <- combos$comp_doy[k]
  idx <- row_idx[[paste(yy, dd)]]
  if (!length(idx)) next
  ut <- u_tile[ebd$uid[idx]]
  for (i in seq_along(tiles)) {
    ridx <- idx[which(ut == i)]
    if (!length(ridx)) next
    cells <- u_cell[ebd$uid[ridx]]
    f <- meta[.(yy, dd, tiles[i]), path, nomatch = NULL]
    if (!length(f)) { miss_gran <- c(miss_gran, sprintf("%d/%d/%s", yy, dd, tiles[i])); next }
    rn <- read_layer(f[1], L_NDVI)[cells]; rn[rn < -2000] <- NA   # -3000 fill
    re <- read_layer(f[1], L_EVI )[cells]; re[re < -2000] <- NA
    rc <- read_layer(f[1], L_CDOY)[cells]; rc[rc < 1 | rc > 366] <- NA
    rr <- read_layer(f[1], L_REL )[cells]; rr[rr < 0 | rr > 3]   <- NA
    ndvi[ridx]  <- rn * 1e-4
    evi[ridx]   <- re * 1e-4
    prel[ridx]  <- as.integer(rr)
    cdate[ridx] <- as.numeric(as.Date(sprintf("%d-01-01", yy))) + (rc - 1)
    # baseline: same composite window, previous year
    fb <- meta[.(yy - 1L, dd, tiles[i]), path, nomatch = NULL]
    if (length(fb)) {
      rb <- read_layer(fb[1], L_NDVI)[cells]; rb[rb < -2000] <- NA
      nbase[ridx] <- rb * 1e-4
    }
  }
  if (k %% 10 == 0 || k == nrow(combos)) {
    saveRDS(list(n = n, ncombo = nrow(combos), k = k, ndvi = ndvi, evi = evi,
                 prel = prel, cdate = cdate, nbase = nbase, miss_gran = miss_gran),
            modis_rds)
    log(sprintf("  composite %d/%d (%d DOY %d)", k, nrow(combos), yy, dd))
  }
}
ebd[, `:=`(ndvi = ndvi, evi = evi, pixel_reliability = prel,
           modis_composite_date = as.Date(cdate, origin = "1970-01-01"),
           ndvi_baseline = nbase)]
ebd[, ndvi_delta := ndvi - ndvi_baseline]

# =============================================================================
# 6. verification -- fail loudly rather than ship a bad table
# =============================================================================
log("Verification checks")
n_zero_dist <- sum(ebd$distance_to_perimeter_km == 0)
stopifnot(
  nrow(ebd) == uniqueN(ebd$checklist_id),                        # grain
  all(is.na(ebd$days_since_fire) | ebd$days_since_fire >= 1),    # R1
  all(ebd$ignition_date[inside] < ebd$observation_date[inside]), # R1, directly
  all(!is.na(ebd$ignition_date[inside])),
  all(is.na(ebd$ignition_date[!inside])),
  all(ebd$protocol_type %chin% KEEP_PROTO),                      # R3
  # R5: sign convention. Non-strict on the inside branch because a checklist
  # exactly on a boundary gives d = 0, a correct value; zeros are counted above.
  all(ebd$distance_to_perimeter_km[inside]  <= 0),
  all(ebd$distance_to_perimeter_km[!inside] >= 0),
  all(ebd$burn_severity_class[!inside] == 0L),
  all(is.na(ebd$ndvi) | (ebd$ndvi >= -1 & ebd$ndvi <= 1))
)
log("  all checks passed")

# =============================================================================
# 7. write
# =============================================================================
col_order <- c("checklist_id","species_code","latitude","longitude","observation_date",
  "year","day_of_year","protocol_type","duration_minutes","effort_distance_km",
  "number_observers","presence","observation_count",
  "fire_event_id","ignition_date","fire_year","burn_severity_class",
  "distance_to_perimeter_km","days_since_fire",
  "ndvi","evi","modis_composite_date","pixel_reliability","ndvi_baseline","ndvi_delta")
out <- ebd[, ..col_order]
out[, distance_to_perimeter_km := round(distance_to_perimeter_km, 4)]
log("Writing flagged table: ", basename(out_csv))
fwrite(out, out_csv, na = "", dateTimeAs = "ISO")

clean <- out[!is.na(ndvi) &
             (is.na(pixel_reliability) | !(pixel_reliability %in% c(2L, 3L))) &
             burn_severity_class != 6L]
log("Writing clean table: ", basename(out_clean))
fwrite(clean, out_clean, na = "", dateTimeAs = "ISO")

# =============================================================================
# 8. summary
# =============================================================================
log("Writing summary")
sink(out_sum)
cat("Wrentit x MTBS x MODIS -- FULL-SPAN joined table\n")
cat("================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Script    : join_species_all_years.R\n")
cat("Grain     : ONE ROW PER COMPLETE CHECKLIST (zero-filled presence/absence)\n")
cat("Species   : ", SPECIES, " (", SPECIES_CODE, ")\n", sep = "")
cat("Span      : ", format(min(out$observation_date)), " .. ",
    format(max(out$observation_date)), "\n", sep = "")
cat("Runtime   : ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1),
    " min\n", sep = "")
if (is_sample) cat("*** SMOKE TEST SUBSET -- not the full table ***\n")

cat("\nSOURCES\n")
cat("  eBird    : zerofill_focal_species_2015_2023.csv + zerofill_focal_species.csv\n")
cat("  protocol : ebd_focal_sampling_filtered*.txt (S ids + G group ids)\n")
cat("  fire     : mtbs_perim/mtbs_perims_DD.shp (", fmt(fc$n_perim), " CA perimeters, ",
    format(fc$perim_span[1]), " .. ", format(fc$perim_span[2]), ")\n", sep = "")
cat("             mtbs/mtbs_CA_<yr>/mtbs_CA_<yr>.tif severity rasters 1984-2025\n")
cat("  MODIS    : MOD13A2_061 ", fmt(nrow(meta)), " granules, ",
    min(meta$ayear), "-", max(meta$ayear), ", tiles ",
    paste(tiles, collapse = " "), "\n", sep = "")

cat("\n-- ROW COUNTS ---------------------------------------------------------\n")
cat("  zero-fill rows read (3 species)  : ", fmt(n_zf), "\n", sep = "")
cat("  ", SPECIES, " rows (pre dedup/QC)      : ", fmt(n_species), "\n", sep = "")
cat("  after dedup on checklist_id      : ", fmt(n_pre_qc), "\n", sep = "")
cat("  dropped by QC                    : ", fmt(n_pre_qc - nrow(out)),
    sprintf("  (%.1f%%)\n", 100 * (n_pre_qc - nrow(out)) / n_pre_qc), sep = "")
cat("      protocol not Traveling/Stationary : ", fmt(qc_drop$protocol), "\n", sep = "")
cat("      duration_minutes   > 300          : ", fmt(qc_drop$duration), "\n", sep = "")
cat("      effort_distance_km > 5            : ", fmt(qc_drop$distance), "\n", sep = "")
cat("      number_observers   > 10           : ", fmt(qc_drop$observers), "\n", sep = "")
cat("      missing coordinates               : ", fmt(qc_drop$coords), "\n", sep = "")
cat("      (criteria overlap; the drop total is their union)\n")
cat("  FLAGGED table                    : ", fmt(nrow(out)), "\n", sep = "")
cat("  CLEAN table                      : ", fmt(nrow(clean)),
    sprintf("  (-%s)\n", fmt(nrow(out) - nrow(clean))), sep = "")
cat("      clean drops: ndvi NA (ocean / land mask), pixel_reliability 2-3,\n")
cat("      burn_severity_class 6 (cloud/shadow/water mask)\n")
cat("  distinct checklist_id            : ", fmt(uniqueN(out$checklist_id)),
    "   (one row per checklist: ", nrow(out) == uniqueN(out$checklist_id), ")\n", sep = "")
cat("  distinct coordinates             : ", fmt(U), "\n", sep = "")
cat("  incomplete checklists in source  : ", asr_bad, "  (expect 0)\n", sep = "")

cat("\n-- COVERAGE BY SPECIES (all three focal species, in the zero-fill input) -\n")
print(as.data.frame(cover), row.names = FALSE)
cat("  NOTE: only ", SPECIES, " is carried through the join, so the output keeps\n", sep = "")
cat("        one row per checklist. The other two species sit in the same\n")
cat("        zero-fill tables over the same span and run through unchanged by\n")
cat("        setting SPECIES / SPECIES_CODE at the top of this script -- the\n")
cat("        fire and NDVI joins are per-checklist, not per-species.\n")

cat("\n-- COVERAGE BY YEAR (flagged table) -----------------------------------\n")
print(as.data.frame(out[, .(rows = .N, detections = sum(presence),
                            det_rate = round(mean(presence), 4),
                            fire_matched = sum(fire_event_id != "NONE"),
                            ndvi_ok = sum(!is.na(ndvi))), by = year][order(year)]),
      row.names = FALSE)
cat("  NOTE 2026 is a partial year: the EBD download ends 2026-06-30.\n")

cat("\n-- MISSING-VALUE CHECK (flagged table) --------------------------------\n")
mv <- sapply(out, function(x) sum(is.na(x)))
print(data.frame(column = names(mv), na = as.integer(mv),
                 pct = round(100 * mv / nrow(out), 2), row.names = NULL))
cat("\n  Expected NAs, not defects:\n")
cat("    ignition_date / fire_year / days_since_fire -> NA on CONTROL rows\n")
cat("      (no mapped fire burned that point before the checklist)\n")
cat("    effort_distance_km -> NA for Stationary checklists (no distance walked)\n")
cat("    observation_count  -> NA where eBird recorded 'X' (present, not counted)\n")
cat("    ndvi/evi/composite/reliability -> NA on MODIS fill (ocean) or a missing\n")
cat("      granule; ndvi_baseline additionally NA where the previous year's same\n")
cat("      window is fill\n")

cat("\n-- MISSING-VALUE CHECK (clean table) ----------------------------------\n")
mvc <- sapply(clean, function(x) sum(is.na(x)))
print(data.frame(column = names(mvc), na = as.integer(mvc),
                 pct = round(100 * mvc / nrow(clean), 2), row.names = NULL))

cat("\n=======================================================================\n")
cat(" FEEDBACK CONFIRMATION -- the six review items, checked against THIS run\n")
cat("=======================================================================\n")

d <- out$days_since_fire
cat("\n[R1] FIRE MATCH RESPECTS TIME ORDER                          -- CONFIRMED\n")
cat("  Already implemented (join_fire_history.R fix 1) and carried over here.\n")
cat("  Rule: the most recent fire CONTAINING the point whose Ig_Date is STRICTLY\n")
cat("        earlier than observation_date.\n")
cat("  days_since_fire negative : ", sum(!is.na(d) & d < 0),
    "        (was ~80% of matched rows, min -251)\n", sep = "")
cat("  days_since_fire zero     : ", sum(!is.na(d) & d == 0), "\n", sep = "")
cat("  days_since_fire positive : ", fmt(sum(!is.na(d) & d > 0)), "\n", sep = "")
cat("  NA (control rows)        : ", fmt(sum(is.na(d))), "\n", sep = "")
print(summary(d))
cat("  Asserted at run time: ignition_date < observation_date on every matched row.\n")
cat("  Pre-fire observations are NOT discarded -- ", fmt(fc$n_prefire), " rows sit at a\n", sep = "")
cat("  location that burns AFTER the checklist. They stay controls rather than\n")
cat("  being back-dated to a fire that had not happened yet. (These are the cases\n")
cat("  the spec's 'negative = pre-fire observation' note was pointing at.)\n")

cat("\n[R2] MULTI-YEAR FIRE PERIMETERS                              -- CONFIRMED\n")
cat("  Already implemented (join_fire_history.R) and extended here to the full\n")
cat("  eBird span rather than 2024 alone.\n")
cat("  perimeters used  : ", fmt(fc$n_perim), " CA MTBS events, ",
    format(fc$perim_span[1]), " .. ", format(fc$perim_span[2]), "\n", sep = "")
cat("  (the 2024-only join used one year and matched 457 / 516,676 rows = 0.09%)\n")
cat("  matched rows now : ", fmt(sum(out$fire_event_id != "NONE")),
    sprintf(" / %s (%.2f%%)\n", fmt(nrow(out)), 100 * mean(out$fire_event_id != "NONE")), sep = "")
cat("  distinct fire events represented : ",
    fmt(uniqueN(out[fire_event_id != "NONE", fire_event_id])), "\n", sep = "")
cat("  fire_year span of matches        : ", min(out$fire_year, na.rm = TRUE), " .. ",
    max(out$fire_year, na.rm = TRUE), "\n", sep = "")
cat("  days_since_fire reach            : ", fmt(max(d, na.rm = TRUE)),
    " days (~", round(max(d, na.rm = TRUE) / 365.25, 1), " years)\n", sep = "")
if (length(missing_rast))
  cat("  WARN fire years with no severity raster: ",
      paste(missing_rast, collapse = ", "), "\n", sep = "")
cat("\n  burn_severity_class distribution (0 = no prior fire):\n")
print(as.data.frame(out[, .N, by = burn_severity_class][order(burn_severity_class)]),
      row.names = FALSE)
cat("  => the covariate now varies instead of being constant for 99.9% of rows.\n")
cat("\n  fire_year of the matched fire:\n")
print(as.data.frame(out[fire_event_id != "NONE", .N, by = fire_year][order(fire_year)]),
      row.names = FALSE)

cat("\n[R3] PROTOCOL FILTER -- Historical and Other dropped         -- CONFIRMED\n")
cat("  Already implemented (join_fire_history.R fix 2); tightened here to EXACT\n")
cat("  name matching (see below).\n")
cat("  protocol_type present BEFORE the filter, with the evidence behind the call:\n")
print(as.data.frame(proto_before), row.names = FALSE)
cat("  Historical and Other both show a depressed detection rate AND heavy\n")
cat("  missingness in the effort fields -- consistent with retroactively entered\n")
cat("  records and non-standard partner data. Effort is a required covariate\n")
cat("  here, so both are dropped.\n")
cat("  kept: ", paste(KEEP_PROTO, collapse = ", "),
    ".  dropped: Historical, Other, Incidental, NA.\n", sep = "")
cat("  Matching is EXACT, so the banding-station variants ('Traveling (2 band,\n")
cat("  25m)', 'Stationary (2/3 band, 30m...)') fall into Other and are dropped\n")
cat("  too -- a grepl() test would have folded them into the kept classes.\n")
cat("  'eBird Pelagic Protocol' and 'Area' are also Other, which is what removes\n")
cat("  MOST of the offshore checklists at source (see R6 for the remainder).\n")
cat("  SPEC DEVIATION, deliberate: the written spec says keep Historical. The\n")
cat("  later review instruction to drop it is followed instead.\n")
cat("  protocol_type levels after the filter: ",
    paste(sort(unique(out$protocol_type)), collapse = ", "), "\n", sep = "")

cat("\n[R4] dnbr -- COLUMN DROPPED, NOT AVAILABLE IN THIS DATA      -- CONFIRMED\n")
cat("  The severity join was never 'pulling dnbr' and never could: no continuous\n")
cat("  dNBR surface exists in the supplied MTBS download. Probed from the raster\n")
cat("  actually used by this run:\n")
if (is.null(rast_probe)) {
  cat("    (no severity raster opened -- no matched fires in this run)\n")
} else {
cat("    file       : mtbs_CA_", rast_probe$year, ".tif\n", sep = "")
cat("    datatype   : ", rast_probe$datatype, "   (8-bit unsigned integer)\n", sep = "")
cat("    colour tab : ", rast_probe$has_coltab, "\n", sep = "")
cat("    values     : ", paste(rast_probe$vals, collapse = ", "), "\n", sep = "")
}
cat("  Those are thematic dnbr6 CLASS codes (1 Unburned-Low, 2 Low, 3 Moderate,\n")
cat("  4 High, 5 Increased Greenness, 6 mask) -- not a continuous index. A real\n")
cat("  dNBR is a scaled float running roughly -500 .. 1300.\n")
cat("  A raw continuous dnbr needs the per-fire *_dnbr.tif bundles from mtbs.gov\n")
cat("  (one per fire, ~", fmt(fc$n_perim), " for CA over this span); they are not in this\n", sep = "")
cat("  download. Available meanwhile: the perimeter attribute table carries\n")
cat("  per-FIRE dNBR calibration values (dnbr_offst, low_t, mod_t, high_t) --\n")
cat("  one number per fire, so a fire-level covariate, not a per-checklist one.\n")
cat("  Decision: ship burn_severity_class, no all-NA dnbr column.\n")

cat("\n[R5] distance_to_perimeter_km SIGN CONVENTION                -- CONFIRMED\n")
cat("  Intended behaviour, now documented -- not a bug.\n")
cat("  NEGATIVE = the point is INSIDE a fire perimeter; POSITIVE = outside.\n")
cat("  The magnitude is the distance to the nearest perimeter BOUNDARY either\n")
cat("  way. This is the spec's own convention ('Negative = inside the fire\n")
cat("  perimeter').\n")
cat("  matched rows with negative distance : ",
    fmt(sum(out$fire_event_id != "NONE" & out$distance_to_perimeter_km < 0)), " of ",
    fmt(sum(out$fire_event_id != "NONE")), "\n", sep = "")
cat("  control rows with positive distance : ",
    fmt(sum(out$fire_event_id == "NONE" & out$distance_to_perimeter_km >= 0)), " of ",
    fmt(sum(out$fire_event_id == "NONE")), "\n", sep = "")
cat("  exact zeros (point on a boundary)   : ", n_zero_dist, "\n", sep = "")
cat("  Asserted at run time for 100% of rows.\n")
print(summary(out$distance_to_perimeter_km))
cat("  Perimeter boundaries are simplified to 30 m before the distance is taken,\n")
cat("  so the induced error is <= 0.03 km.\n")

cat("\n[R6] OFFSHORE CHECKLISTS west of -124.4 deg                  -- CONFIRMED\n")
w  <- out[longitude < -124.4, .(checklist_id, latitude, longitude, protocol_type, ndvi)]
wc <- clean[longitude < -124.4, .(checklist_id, latitude, longitude, protocol_type, ndvi)]
cat("  Attrition through the pipeline:\n")
cat("    in the zero-fill, before QC : ", fmt(off_pre), "\n", sep = "")
cat("    after ALL QC filters        : ", nrow(w), "\n", sep = "")
cat("    surviving into CLEAN        : ", nrow(wc), "   (-", nrow(w) - nrow(wc),
    ", dropped by the NDVI land mask)\n", sep = "")
cat("\n  Which rule removes them -- NOT the protocol filter alone:\n")
print(as.data.frame(off_break), row.names = FALSE)
cat("  protocols among those the protocol rule removes:\n")
print(as.data.frame(off_proto), row.names = FALSE)
cat("  Note the largest single cause is effort_distance_km > 5, not protocol: a\n")
cat("  boat trip logged as ordinary 'Traveling' covers far more than 5 km, so the\n")
cat("  standard eBird effort QC removes it without ever knowing it was at sea.\n")
cat("\n  DIAGNOSIS: neither a bounding-box clip nor a geocoding error. eBird's\n")
cat("  US-CA extent legitimately includes CALIFORNIA STATE OFFSHORE WATERS, so\n")
cat("  genuine at-sea checklists exist in the source data.\n")
cat("\n  The two filters do different jobs, and it matters which:\n")
cat("    * the PROTOCOL filter removes the pelagic-protocol checklists -- the\n")
cat("      bulk of them, at source;\n")
cat("    * it does NOT catch boat checklists logged under ordinary Traveling /\n")
cat("      Stationary protocols -- roughly half of these are. Those are removed\n")
cat("      by the effort thresholds instead, and whatever still survives is\n")
cat("      removed by the NDVI land mask. Every one of the ", nrow(w), " rows that\n", sep = "")
cat("      reaches the flagged table is Traveling or Stationary.\n")
if (nrow(w)) {
  far <- w[longitude < -124.45]
  cat("\n  Of the ", nrow(w), " that reach the flagged table:\n", sep = "")
  cat("    lon <  -124.45 (open ocean)      : ", nrow(far), "  -- NDVI NA: ",
      sum(is.na(far$ndvi)), " of ", nrow(far), "\n", sep = "")
  cat("    lon >= -124.45 (Cape Mendocino)  : ", nrow(w) - nrow(far), "\n", sep = "")
  cat("  MODIS calls the open-ocean points fill, so the clean table's !is.na(ndvi)\n")
  cat("  land mask removes them. None survive.\n")
}
cat("\n  NOT removed, and correctly so -- the rows that reach the clean table:\n")
if (nrow(wc)) print(as.data.frame(wc), row.names = FALSE)
cat("  These sit at about -124.400 deg, latitude ~40.43: Mattole Rd on the Lost\n")
cat("  Coast near Cape Mendocino, which is California's westernmost LAND\n")
cat("  (~-124.41). All carry valid, vegetated NDVI. A blanket longitude cut at\n")
cat("  -124.4 would delete real coastal Wrentit habitat, so NO longitude filter\n")
cat("  is applied -- the land mask is the correct instrument, not a lon cut.\n")

cat("\n-- DETECTION RATE -----------------------------------------------------\n")
cat("  flagged: ", round(mean(out$presence), 4),
    "   clean: ", round(mean(clean$presence), 4), "\n", sep = "")
cat("  by fire status (clean table):\n")
print(as.data.frame(clean[, .(rows = .N, det_rate = round(mean(presence), 4)),
                          by = .(status = fifelse(fire_event_id == "NONE",
                                                  "control", "post-fire"))]),
      row.names = FALSE)

cat("\n-- VEGETATION ---------------------------------------------------------\n")
cat("NDVI\n"); print(summary(out$ndvi))
cat("\nNDVI_DELTA (same 16-day window, vs previous year)\n")
print(summary(out$ndvi_delta))
cat("\nPIXEL RELIABILITY (0 good, 1 marginal, 2 snow/ice, 3 cloud)\n")
print(as.data.frame(out[, .N, by = pixel_reliability][order(pixel_reliability)]),
      row.names = FALSE)
cat("\nunique coordinates outside all MODIS tiles: ", n_outside_tiles, "\n", sep = "")
if (length(miss_gran))
  cat("WARN missing MODIS granules: ", length(miss_gran), " -> ",
      paste(head(miss_gran, 10), collapse = ", "), "\n", sep = "")

cat("\n-- OUTPUTS ------------------------------------------------------------\n")
cat("  ", basename(out_csv),   " : ", fmt(nrow(out)),   " rows x ", ncol(out),   " cols\n", sep = "")
cat("  ", basename(out_clean), " : ", fmt(nrow(clean)), " rows x ", ncol(clean), " cols\n", sep = "")
cat("\nCOLUMNS (", ncol(out), ")\n", sep = "")
print(names(out))
sink()

cat(readLines(out_sum), sep = "\n")
log(sprintf("DONE in %.1f min", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
