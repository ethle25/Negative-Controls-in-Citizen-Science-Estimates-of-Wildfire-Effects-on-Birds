# =============================================================================
# OPTION B: can severity CLASSES be turned back into continuous dNBR?
#
# The download has no continuous dNBR raster -- only the thematic dnbr6 classes
# (1-6). But the perimeter attributes carry each fire's own dNBR class
# THRESHOLDS, so a class can be inverted into a dNBR RANGE:
#
#     class 5  dNBR <  incgreen_t                (increased greenness)
#     class 1  incgreen_t <= dNBR < low_t        (unburned to low)
#     class 2  low_t      <= dNBR < mod_t        (low)
#     class 3  mod_t      <= dNBR < high_t       (moderate)
#     class 4  dNBR >= high_t                    (high)
#
# Taking the MIDPOINT of each range gives an approximate continuous dNBR that is
# calibrated with that fire's own thresholds -- which also puts every fire on a
# common physical scale instead of 1,407 separate grading schemes.
#
# THIS SCRIPT TESTS WHETHER THAT WORKS, before any of it is used.
#
# The CBI field-plot dataset (bsp_FODpoints_DD) is the ground truth: 53,383 plots
# that carry BOTH the measured dNBR (dnbr_val) AND their fire's thresholds. So
# for every plot we can:
#   1. classify its measured dNBR using its own fire's thresholds
#   2. reconstruct a midpoint dNBR from that class
#   3. compare the reconstruction with the measurement
# and separately check whether the reconstruction still tracks CBI -- the
# field-crew severity score, which is independent of any satellite index.
#
# Class 4 is open-ended (no upper bound), so its midpoint needs an assumption:
# it is set to high_t + half the width of class 3. That assumption is tested
# here like everything else.
#
# WRITES NEW FILES ONLY.
# =============================================================================
suppressMessages({ library(data.table); library(sf) })
setDTthreads(0)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
out_sum <- file.path(dir_, "dnbr_validation_summary.txt")
out_csv <- file.path(dir_, "dnbr_class_midpoints.csv")

say("Reading CBI field plots")
p <- st_read(file.path(dir_, "mtbs_bsp", "bsp_FODpoints_DD.shp"), quiet = TRUE)
d <- as.data.table(st_drop_geometry(p))
say("  plots: ", format(nrow(d), big.mark = ","),
    " | with measured dNBR: ", format(sum(!is.na(d$dnbr_val)), big.mark = ","))

# TWO problems with the points table:
#   1. its threshold columns are filled with zeros -- only the PERIMETER table
#      has the real per-fire thresholds
#   2. its event_id is a plot/project code ("cbi_quaking_00003"), NOT an MTBS
#      fire id ("UT41764111488..."), so the two tables share no key at all
# The link has to be made SPATIALLY: find which fire perimeter each plot falls
# inside, restricted to fires that ignited before the plot was surveyed.
say("Linking plots to fires SPATIALLY (no shared id exists)")
pts <- p[!is.na(p$dnbr_val), ]
say("  plots with measured dNBR: ", format(nrow(pts), big.mark = ","))
per <- st_read(file.path(dir_, "mtbs_perim", "mtbs_perims_DD.shp"), quiet = TRUE)
per <- st_make_valid(per)
sf::sf_use_s2(FALSE)
hit <- st_within(st_transform(pts, st_crs(per)), per)
say("  plots inside >=1 perimeter: ", format(sum(lengths(hit) > 0), big.mark = ","))

pd <- as.data.table(st_drop_geometry(pts))
pd[, .row := .I]
pdate <- suppressWarnings(as.Date(pd$field_date))
perd  <- as.data.table(st_drop_geometry(per))
perig <- suppressWarnings(as.Date(perd$ig_date))

# for each plot keep the most recent fire that ignited before the field survey
pick <- rep(NA_integer_, nrow(pd))
for (i in which(lengths(hit) > 0)) {
  cand <- hit[[i]]
  if (!is.na(pdate[i])) cand <- cand[!is.na(perig[cand]) & perig[cand] <= pdate[i]]
  if (length(cand)) pick[i] <- cand[which.max(perig[cand])]
}
say("  plots matched to a prior fire: ", format(sum(!is.na(pick)), big.mark = ","))
pd[, fire := pick]
pd <- pd[!is.na(fire)]
for (k in c("low_t","mod_t","high_t","incgreen_t","dnbr_offst"))
  set(pd, NULL, k, perd[[k]][pd$fire])
pd[, mtbs_event := perd$event_id[pd$fire]]
d <- pd

num <- function(x) { v <- suppressWarnings(as.numeric(as.character(x)))
                     v[!is.finite(v) | abs(v) >= 9999] <- NA_real_; v }
for (k in c("dnbr_val","dnbr_offst","low_t","mod_t","high_t","incgreen_t","cbi"))
  if (k %in% names(d)) set(d, NULL, k, num(d[[k]]))

d[, ca := grepl("^CA", toupper(mtbs_event))]
ok <- d[!is.na(dnbr_val) & !is.na(low_t) & !is.na(mod_t) & !is.na(high_t) &
        low_t < mod_t & mod_t < high_t]
say("  usable plots (measured dNBR + ordered thresholds): ",
    format(nrow(ok), big.mark = ","), "  (CA: ", format(sum(ok$ca), big.mark=","), ")")

# incgreen_t is missing for some fires; fall back to a conventional -100
ok[, ig_t := fifelse(is.na(incgreen_t), -100, incgreen_t)]

# 1. classify the MEASURED dNBR with the fire's own thresholds
ok[, cls := fifelse(dnbr_val <  ig_t,  5L,
            fifelse(dnbr_val <  low_t, 1L,
            fifelse(dnbr_val <  mod_t, 2L,
            fifelse(dnbr_val <  high_t,3L, 4L))))]

# 2. reconstruct the class midpoint (the quantity a raster-only workflow gets)
ok[, mid := fifelse(cls == 1L, (ig_t + low_t) / 2,
            fifelse(cls == 2L, (low_t + mod_t) / 2,
            fifelse(cls == 3L, (mod_t + high_t) / 2,
            fifelse(cls == 4L, high_t + (high_t - mod_t) / 2, NA_real_))))]
ok[, mid_adj := mid - fifelse(is.na(dnbr_offst), 0, dnbr_offst)]
ok[, val_adj := dnbr_val - fifelse(is.na(dnbr_offst), 0, dnbr_offst)]
ok[, err := mid - dnbr_val]

v <- ok[cls %in% 1:4 & !is.na(mid)]
say("  plots in severity classes 1-4: ", format(nrow(v), big.mark = ","))

fit <- function(x, y) { k <- is.finite(x) & is.finite(y)
  if (sum(k) < 30) return(c(NA, NA, NA))
  c(cor(x[k], y[k]), sqrt(mean((x[k]-y[k])^2)), mean(x[k]-y[k])) }

overall <- fit(v$mid, v$dnbr_val)
byclass <- v[, .(plots = .N,
                 measured_median = round(median(dnbr_val), 1),
                 midpoint = round(median(mid), 1),
                 bias = round(mean(mid - dnbr_val), 1),
                 rmse = round(sqrt(mean((mid - dnbr_val)^2)), 1),
                 sd_within = round(sd(dnbr_val), 1)), by = cls][order(cls)]

# 3. does the reconstruction preserve the relationship with FIELD severity (CBI)?
cb <- v[!is.na(cbi)]
cbi_meas <- fit(cb$dnbr_val, cb$cbi)
cbi_recon <- fit(cb$mid, cb$cbi)
cbi_cls  <- fit(as.numeric(cb$cls), cb$cbi)

# CA-only replication
vca <- v[ca == TRUE]
ca_fit <- fit(vca$mid, vca$dnbr_val)

fwrite(unique(ok[, .(mtbs_event, ig_t, low_t, mod_t, high_t, dnbr_offst)]), out_csv)

sink(out_sum)
cat("Can severity CLASSES be inverted back into continuous dNBR?\n")
cat("===========================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Ground truth: MTBS CBI field plots (bsp_FODpoints_DD), which carry both a\n")
cat("              MEASURED dNBR and their fire's class thresholds.\n")
cat("Method    : classify the measured dNBR with the fire's own thresholds, then\n")
cat("            reconstruct the class MIDPOINT and compare with the measurement.\n\n")
cat("  plots spatially matched to a fire     : ", format(nrow(d), big.mark=","), "\n", sep="")
cat("  usable (dNBR + ordered thresholds)   : ", format(nrow(ok), big.mark=","), "\n", sep="")
cat("  in severity classes 1-4              : ", format(nrow(v), big.mark=","), "\n", sep="")
cat("  of which California                  : ", format(nrow(vca), big.mark=","), "\n", sep="")

cat("\n-- 1. HOW CLOSE IS THE RECONSTRUCTION TO MEASURED dNBR? --------------\n")
cat(sprintf("  correlation midpoint vs measured : %.3f\n", overall[1]))
cat(sprintf("  RMSE                             : %.1f dNBR units\n", overall[2]))
cat(sprintf("  bias (midpoint - measured)       : %+.1f\n", overall[3]))
cat(sprintf("  California only, correlation     : %.3f  (RMSE %.1f)\n", ca_fit[1], ca_fit[2]))
cat("\n  For scale: dNBR runs roughly -200 (greener) to 1300 (total consumption),\n")
cat("  and the class thresholds sit near 75 / 296 / 550.\n")

cat("\n-- 2. BY SEVERITY CLASS ---------------------------------------------\n")
cat("  'sd_within' is the spread of REAL dNBR inside each class -- the\n")
cat("  information a midpoint necessarily throws away.\n\n")
print(as.data.frame(byclass), row.names = FALSE)

cat("\n-- 3. DOES IT STILL TRACK FIELD-MEASURED SEVERITY (CBI)? -------------\n")
cat("  CBI is scored by field crews on the ground, independent of any satellite\n")
cat("  index. If the reconstruction correlates with CBI nearly as well as the\n")
cat("  measured dNBR does, little of value was lost.\n\n")
cat(sprintf("  plots with CBI                       : %s\n", format(nrow(cb), big.mark=",")))
cat(sprintf("  measured dNBR      vs CBI  r = %.3f\n", cbi_meas[1]))
cat(sprintf("  reconstructed dNBR vs CBI  r = %.3f\n", cbi_recon[1]))
cat(sprintf("  raw class (1-4)    vs CBI  r = %.3f\n", cbi_cls[1]))
cat("\n  The third line is the baseline: what you get WITHOUT reconstructing.\n")
cat("  Reconstruction is worth doing only if line 2 beats line 3.\n")

cat("\n-- VERDICT ----------------------------------------------------------\n")
gain <- cbi_recon[1] - cbi_cls[1]
cat(sprintf("  reconstruction vs raw class, gain in CBI correlation: %+.3f\n", gain))
# braces required: an unbraced else-if branch cannot hold two statements
if (is.na(gain)) {
  cat("  INSUFFICIENT DATA to judge.\n")
} else if (gain > 0.02) {
  cat("  USE IT: the midpoint carries more severity information than the class.\n")
} else if (gain > -0.05) {
  cat("  DO NOT BOTHER: the midpoint is no better than the plain class number,\n")
  cat("  and here it is slightly worse. Converting classes to pseudo-continuous\n")
  cat("  dNBR adds false precision -- the spread of real dNBR WITHIN a class is\n")
  cat("  large (see sd_within above), so a single midpoint misrepresents it.\n")
  cat("  Keep using burn_severity_class / the per-class area fractions.\n")
} else {
  cat("  DO NOT USE: reconstruction actively loses information.\n")
}

cat("\n-- CAVEATS ----------------------------------------------------------\n")
cat("  * Field plots are not a random sample of burned area -- they are placed\n")
cat("    for calibration, so they over-represent well-studied fires.\n")
cat("  * Class 4 is open-ended; its midpoint uses high_t + half the class-3\n")
cat("    width. Check the class-4 bias row: a large positive or negative bias\n")
cat("    there means that assumption is wrong.\n")
cat("  * These plots are national. The CA-only correlation is reported above as\n")
cat("    a check that California behaves like the whole.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
