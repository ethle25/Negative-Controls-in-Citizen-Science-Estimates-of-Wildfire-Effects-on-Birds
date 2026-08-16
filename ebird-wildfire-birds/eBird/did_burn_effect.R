# =============================================================================
# What does fire actually do to Wrentit detection? -- within-cell estimates
#
# Separate from model_burn_effect.R on purpose. That script answers a PREDICTIVE
# question (how much does burn history improve forecasts). This one answers a
# CAUSAL-flavoured question (what happens to detection after a cell burns), and
# the two must not be quoted interchangeably.
#
# THE PROBLEM WITH THE RAW CONTRAST
#   Burned cells show far higher detection than never-burned cells. That is
#   mostly habitat: MTBS fires concentrate in chaparral, which is precisely
#   Wrentit habitat, while never-burned cells include desert, farmland and city.
#   Comparing burned vs unburned CELLS compares chaparral vs California.
#
# THE DESIGN
#   Treated : cells whose burn status flips 0 -> 1 inside the study window, so
#             the same cell is observed both before and after its own fire.
#             Cell identity -- and therefore habitat -- cancels.
#   Control : cells that never burn in the window, each assigned a PSEUDO event
#             year drawn from the treated cells' event-year distribution. This
#             is what makes a difference-in-differences possible: without it the
#             controls have no "before"/"after" to speak of. (The previous
#             attempt used min(event year) = 2015 as a single cutoff, which left
#             the controls with no pre-period at all and no DiD.)
#
#   Outcome : det_rate = detections / checklists, i.e. ALREADY EFFORT-NORMALISED.
#             Using any_detection here would confound with how many checklists
#             a cell-time received.
#
#   Estimators
#     1. 2x2 DiD, checklist-weighted
#     2. two-way fixed effects (cell + year), checklist-weighted, via iterative
#        demeaning -- absorbs every time-invariant cell trait and every
#        statewide year shock
#     3. EVENT-TIME profile: detection by years since the cell's own fire,
#        which is the part a single before/after average throws away
#
#   Uncertainty: a cell-level bootstrap, because cells -- not cell-time units --
#   are the independent sampling unit.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)
set.seed(1)

say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }
fmt <- function(x) format(x, big.mark = ",")

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")

source(file.path(dir_, "panel_utils.R"))
PREFIX  <- cfg_default_prefix()
SPECIES <- cfg_default_species()
DOSE    <- Sys.getenv("DID_DOSE", "frac_max")

# Output is per species. It used to be a single burn_did_summary.txt, so running
# a second species silently overwrote the first one's report with numbers
# labelled "Wrentit".
out_sum <- file.path(dir_, sprintf("burn_did_summary_%s.txt", PREFIX))

say("Reading frame")
# MUST be the PREFIX species' frame. This line hardcoded wrentit_cells_5km_16d.csv
# while build_panel() below correctly used PREFIX, so EBIRD_PREFIX=olsfly produced
# sections 1-3 from olsfly and section 4 (severity) from Wrentit -- one report,
# two species, no error.
d <- fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", PREFIX)), showProgress = FALSE,
           select = c("cell_id","year","comp_doy","ever_burned","yrs_since_burn",
                      "n_checklists","n_detections","any_detection","burn_frac_cum",
                      "max_sev_prior"))

# --- who is treated, who is control -------------------------------------------
# Shared definition (panel_utils.R): treated = burn status flips inside the
# window; controls never burn and get a pseudo event year drawn from the treated
# distribution, so the same before/after split is applied to a group that never
# actually burned. Identical to dose_burn_effect.R by construction.
#
# The dose is irrelevant to THIS script -- treatment is binary -- but build_panel
# needs one, and DID_DOSE lets the binary subset be defined by any registry
# measure if that is ever wanted.
cy <- build_panel(PREFIX, dose = DOSE)
cy[, evt := year - ev_year]
# later sections work off the raw 16-day frame, so they still need the cell lists
treated <- unique(cy[treat == 1]$cell_id)
control <- unique(cy[treat == 0]$cell_id)
ev      <- unique(cy[treat == 1, .(cell_id, ev_year)])
say("treated cells: ", fmt(uniqueN(cy[treat == 1]$cell_id)),
    " | control cells: ", fmt(uniqueN(cy[treat == 0]$cell_id)))
say("treated event years: ",
    paste(range(cy[treat == 1]$ev_year), collapse = " .. "))
say("cell-year panel: ", fmt(nrow(cy)), " rows")

# =============================================================================
# 1. 2x2 difference-in-differences, checklist-weighted
# =============================================================================
cell2 <- cy[, .(rate = sum(det) / sum(chk), chk = sum(chk), cells = uniqueN(cell_id)),
            by = .(treat, post)][order(-treat, post)]
g <- function(tr, po) cell2[treat == tr & post == po]$rate
did_22 <- (g(1,1) - g(1,0)) - (g(0,1) - g(0,0))

# =============================================================================
# 2. two-way fixed effects (cell + year), weighted, by iterative demeaning
# =============================================================================
say("Two-way fixed effects")
tw <- copy(cy)[, .(cell_id, year, rate = as.numeric(rate),
                   post = as.numeric(post), w = as.numeric(w))]
# post MUST be double: data.table `:=` into an integer column truncates the
# demeaned values to 0, which silently zeroes the denominator and yields NaN
for (i in 1:40) {
  tw[, rate := rate - weighted.mean(rate, w), by = cell_id]
  tw[, post := post - weighted.mean(post, w), by = cell_id]
  tw[, rate := rate - weighted.mean(rate, w), by = year]
  tw[, post := post - weighted.mean(post, w), by = year]
}
beta_twfe <- sum(tw$w * tw$post * tw$rate) / sum(tw$w * tw$post^2)

# cell-level bootstrap: cells are the independent unit, not cell-years
say("Bootstrapping (cell-level, 200 draws)")
cells_all <- unique(cy$cell_id)
boot <- numeric(200)
for (b in seq_along(boot)) {
  smp <- data.table(cell_id = sample(cells_all, length(cells_all), replace = TRUE))
  smp[, newid := .I]
  bb <- merge(smp, cy, by = "cell_id", allow.cartesian = TRUE)
  t2 <- bb[, .(cell_id = newid, year, rate = as.numeric(rate),
               post = as.numeric(post), w = as.numeric(w))]
  for (i in 1:25) {
    t2[, rate := rate - weighted.mean(rate, w), by = cell_id]
    t2[, post := post - weighted.mean(post, w), by = cell_id]
    t2[, rate := rate - weighted.mean(rate, w), by = year]
    t2[, post := post - weighted.mean(post, w), by = year]
  }
  den <- sum(t2$w * t2$post^2)
  boot[b] <- if (den > 0) sum(t2$w * t2$post * t2$rate) / den else NA_real_
}
ci <- quantile(boot, c(0.025, 0.975), na.rm = TRUE)

# =============================================================================
# 3. event-time profile
# =============================================================================
prof <- cy[evt >= -6 & evt <= 10,
           .(rate = sum(det) / sum(chk), chk = sum(chk), cells = uniqueN(cell_id)),
           by = .(treat, evt)][order(treat, evt)]
pw <- dcast(prof, evt ~ treat, value.var = c("rate", "chk", "cells"))
setnames(pw, c("evt","rate_ctrl","rate_treat","chk_ctrl","chk_treat","cells_ctrl","cells_treat"))
# normalise each arm to its own pre-period (evt in -3..-1), then difference
base_t <- prof[treat == 1 & evt %in% -3:-1, sum(rate * chk) / sum(chk)]
base_c <- prof[treat == 0 & evt %in% -3:-1, sum(rate * chk) / sum(chk)]
pw[, `:=`(treat_vs_base = rate_treat - base_t,
          ctrl_vs_base  = rate_ctrl  - base_c)]
pw[, did_evt := treat_vs_base - ctrl_vs_base]

# =============================================================================
# 4. by severity of the cell's own fire
# =============================================================================
sev <- d[cell_id %in% treated]
sev <- merge(sev, ev, by = "cell_id")
sev[, post := as.integer(year >= ev_year)]
sev_cell <- d[cell_id %in% treated & ever_burned == 1,
              .(sev_max = max(max_sev_prior), frac = max(burn_frac_cum)), by = cell_id]
sev <- merge(sev, sev_cell, by = "cell_id")
sev[, sev_band := cut(sev_max, c(-0.1, 1.5, 2.5, 3.5, 4.5),
                      labels = c("1 unburned-low","2 low","3 moderate","4 high"))]
sev_tab <- sev[, .(cells = uniqueN(cell_id), chk = sum(n_checklists),
                   rate = sum(n_detections) / sum(n_checklists)),
               by = .(sev_band, period = ifelse(post == 1, "after", "before"))]
sev_w <- dcast(sev_tab, sev_band ~ period, value.var = c("rate","chk","cells"))
sev_w[, change := rate_after - rate_before]

# =============================================================================
# report
# =============================================================================
sink(out_sum)
cat("What does fire do to ", SPECIES, " detection? -- within-cell estimates\n", sep = "")
cat("===============================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Outcome   : det_rate = detections / checklists (already effort-normalised)\n")
cat("Unit      : cell-year, checklist-weighted\n")

cat("\n-- DESIGN ------------------------------------------------------------\n")
cat("  treated cells (burn status flips inside the window) : ", fmt(length(treated)), "\n", sep = "")
cat("  control cells (never burn in the window)            : ", fmt(length(control)), "\n", sep = "")
cat("  treated event years                                 : ",
    paste(range(ev$ev_year), collapse = " .. "), "\n", sep = "")
cat("  controls carry a PSEUDO event year drawn from the treated distribution,\n")
cat("  so both arms are split into before/after on the same calendar pattern.\n")

cat("\n-- 1. RAW 2x2 DIFFERENCE-IN-DIFFERENCES ------------------------------\n")
print(as.data.frame(cell2[, .(arm = ifelse(treat == 1, "treated (burns)", "control (never burns)"),
                              period = ifelse(post == 1, "after", "before"),
                              cells, checklists = chk, det_rate = round(rate, 4))]),
      row.names = FALSE)
cat("\n  treated change : ", sprintf("%+.4f", g(1,1) - g(1,0)), "\n", sep = "")
cat("  control change : ", sprintf("%+.4f", g(0,1) - g(0,0)), "\n", sep = "")
cat("  DiD            : ", sprintf("%+.4f", did_22), " detections per checklist\n", sep = "")

cat("\n-- 2. TWO-WAY FIXED EFFECTS (cell + year) ----------------------------\n")
cat("  beta(post)     : ", sprintf("%+.4f", beta_twfe), "\n", sep = "")
cat("  95% CI (cell bootstrap, 200 draws) : [",
    sprintf("%+.4f", ci[1]), ", ", sprintf("%+.4f", ci[2]), "]\n", sep = "")
cat("  baseline det_rate in treated cells BEFORE their fire : ",
    sprintf("%.4f", g(1,0)), "\n", sep = "")
cat("  => relative change : ", sprintf("%+.1f%%", 100 * beta_twfe / g(1,0)), "\n", sep = "")
cat("\n  This absorbs every time-invariant property of a cell (habitat, access,\n")
cat("  observer community) and every statewide year effect. What is left is the\n")
cat("  within-cell change around its own fire, net of the control trend.\n")

cat("\n-- 3. EVENT-TIME PROFILE ---------------------------------------------\n")
cat("  det_rate by years since the cell's own fire (evt 0 = fire year);\n")
cat("  each arm is expressed relative to its own evt -3..-1 baseline.\n\n")
print(as.data.frame(pw[, .(evt, treat_rate = round(rate_treat, 4),
                           ctrl_rate = round(rate_ctrl, 4),
                           treat_vs_base = round(treat_vs_base, 4),
                           ctrl_vs_base = round(ctrl_vs_base, 4),
                           did = round(did_evt, 4),
                           treat_cells = cells_treat)]), row.names = FALSE)
cat("\n  Pre-period rows (evt < 0) are the parallel-trends check: if the treated\n")
cat("  and control 'vs_base' columns track each other BEFORE the fire, the DiD\n")
cat("  is credible. If they diverge pre-fire, it is not.\n")

cat("\n-- 4. BY SEVERITY OF THE CELL'S OWN FIRE -----------------------------\n")
print(as.data.frame(sev_w[, .(sev_band, cells = cells_after,
                              rate_before = round(rate_before, 4),
                              rate_after = round(rate_after, 4),
                              change = round(change, 4))]), row.names = FALSE)

cat("\n-- INTERPRETATION ----------------------------------------------------\n")
cat("  The cross-sectional contrast (burned cells vs never-burned cells) and\n")
cat("  the within-cell estimate answer different questions and can easily carry\n")
cat("  OPPOSITE signs. The cross-sectional gap is dominated by habitat; the\n")
cat("  within-cell estimate is the one to quote for 'what does fire do'.\n")
cat("\n  Remaining limitations, none of them fixed by this design:\n")
cat("   * MTBS under-maps 2023-2025, so recent fires are missing and some\n")
cat("     'control' cells may have burned without being recorded\n")
cat("   * cells burning early in the window have a very short pre-period\n")
cat("   * detection is not occupancy: a change here can be birds leaving OR\n")
cat("     birds becoming harder to detect in changed vegetation\n")
cat("   * observer behaviour may change after a fire (access, interest)\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
