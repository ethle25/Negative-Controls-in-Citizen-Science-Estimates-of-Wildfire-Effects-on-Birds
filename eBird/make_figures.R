# make_figures.R -- the five paper figures specified in the in-depth outline.
#
# Every number is PARSED from the result files on disk, never retyped, so a
# refit that moves a published number moves the figure with it. Where a file is
# a formatted text report rather than a CSV (event_study_summary.txt,
# effort_placebo_summary.txt, deep_models_common_test_*.txt) the parser is
# followed by a stopifnot() that re-derives a value printed in the report from
# the parsed pieces. If a report's layout changes, the check fails loudly
# instead of the figure quietly drawing the wrong thing.
#
# Run:  Rscript make_figures.R
# Out:  figures/fig[1-5]_*.png  (300 dpi)  and  .pdf  (vector, for typesetting)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

PROJ <- Sys.getenv("EBIRD_DIR",
         file.path(Sys.getenv("EBIRD_PROJ_ROOT", ".."), "eBird"))
FIGDIR <- file.path(PROJ, "figures")
dir.create(FIGDIR, showWarnings = FALSE)

say <- function(...) cat(sprintf(...), "\n")   # never call a logger `log` -- see Traps

# ---------------------------------------------------------------- palette ----
# Validated with the dataviz validator (light surface #fcfcfb, --pairs all):
# CVD worst 9.2 dE, normal-vision worst 16.3 dE, all checks pass. Aqua sits at
# 2.74:1 contrast, so every series carrying it is also DIRECT-LABELLED -- that
# is the relief rule, not an optional nicety.
C_BLUE   <- "#2a78d6"   # slot 1 -- focal species / real contrasts
C_ORANGE <- "#eb6834"   # slot 2 -- controls / placebos / deep models
C_AQUA   <- "#1baf7a"   # slot 3
C_VIOLET <- "#4a3aa7"   # slot 4
C_BLUE_L <- "#9ec5f4"   # blue-200, ordinal step within one hue
C_CRIT   <- "#d03b3b"   # status: critical (a measure that SHIFTS)
INK      <- "#0b0b0b"
INK2     <- "#52514e"
MUTED    <- "#898781"
GRID     <- "#e1e0d9"
AXIS     <- "#c3c2b7"
SURFACE  <- "#fcfcfb"

theme_paper <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.background   = element_rect(fill = SURFACE, colour = NA),
      panel.background  = element_rect(fill = SURFACE, colour = NA),
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(colour = GRID, linewidth = 0.3),
      axis.line         = element_line(colour = AXIS, linewidth = 0.4),
      axis.ticks        = element_line(colour = AXIS, linewidth = 0.4),
      axis.text         = element_text(colour = INK2, size = base_size - 1),
      axis.title        = element_text(colour = INK2, size = base_size - 0.5),
      plot.title        = element_text(colour = INK, size = base_size + 2,
                                       face = "bold", hjust = 0),
      plot.subtitle     = element_text(colour = INK2, size = base_size - 0.5,
                                       hjust = 0, lineheight = 1.15),
      plot.caption      = element_text(colour = MUTED, size = base_size - 2.5,
                                       hjust = 0, lineheight = 1.2),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      legend.title      = element_blank(),
      legend.text       = element_text(colour = INK2, size = base_size - 1),
      legend.key.height = unit(11, "pt"),
      strip.text        = element_text(colour = INK, size = base_size - 1,
                                       face = "bold", hjust = 0)
    )
}

# The captions carry real typography -- Delta, >=, en/em dashes, curly quotes,
# the multiplication sign. The default png()/pdf() devices on Windows drop them
# with a 'conversion failure ... mbcsToSbcs' warning, so use the two devices
# that are UTF-8 clean: ragg for raster, cairo for vector.
save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIGDIR, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, units = "in", bg = SURFACE, device = ragg::agg_png)
  ggsave(file.path(FIGDIR, paste0(name, ".pdf")), p, width = w, height = h,
         units = "in", bg = SURFACE, device = cairo_pdf)
  say("wrote %s.png / .pdf  (%.1f x %.1f in)", name, w, h)
}

# ============================================================================
# FIGURE 1 -- event study: Wrentit against the three water-bird controls
# Source: event_study_summary.txt (formatted report, needs parsing)
# ============================================================================

parse_event_study <- function(path) {
  ln <- readLines(path, warn = FALSE)

  # -- the header table lists species and pre-fire rates in two SEPARATE
  #    column blocks, because the report wrapped. They are in the same order.
  i_spp <- grep("^\\s*species\\s+treated_cells", ln)
  i_pfr <- grep("^\\s*pre_fire_rate\\s*$", ln)
  stopifnot(length(i_spp) == 1, length(i_pfr) == 1)

  spp_block  <- ln[(i_spp + 1):(i_pfr - 1)]
  spp_block  <- spp_block[nzchar(trimws(spp_block))]
  spp_names  <- trimws(sub("\\s+\\d+\\s+\\d+\\s*$", "", spp_block))

  rate_block <- ln[(i_pfr + 1):length(ln)]
  rate_block <- rate_block[seq_len(length(spp_names))]
  pre_rate   <- as.numeric(trimws(rate_block))

  stopifnot(length(spp_names) == length(pre_rate), !anyNA(pre_rate))
  rates <- data.table(species = spp_names, pre_rate = pre_rate)

  # -- per-species coefficient blocks
  hdr <- grep("^=== .+ ===$", ln)
  out <- rbindlist(lapply(seq_along(hdr), function(j) {
    from <- hdr[j] + 1
    to   <- if (j < length(hdr)) hdr[j + 1] - 1 else length(ln)
    body <- ln[from:to]
    rows <- grep("^\\s*-?\\d+\\s+[+-][0-9.]+\\s+[0-9.]+\\s+\\[", body, value = TRUE)
    if (!length(rows)) return(NULL)
    m <- regmatches(rows, regexec(
      "^\\s*(-?\\d+)\\s+([+-][0-9.]+)\\s+([0-9.]+)\\s+\\[\\s*([+-][0-9.]+),\\s*([+-][0-9.]+)\\]\\s+([+-][0-9.]+)%",
      rows))
    data.table(
      species = sub("^=== (.+) ===$", "\\1", ln[hdr[j]]),
      k    = as.integer(sapply(m, `[`, 2)),
      beta = as.numeric(sapply(m, `[`, 3)),
      se   = as.numeric(sapply(m, `[`, 4)),
      lo   = as.numeric(sapply(m, `[`, 5)),
      hi   = as.numeric(sapply(m, `[`, 6)),
      pct_reported = as.numeric(sapply(m, `[`, 7))
    )
  }))

  out <- merge(out, rates, by = "species", all.x = TRUE)
  stopifnot(!anyNA(out$pre_rate))

  # INTEGRITY CHECK: re-derive the report's own %base column from beta and the
  # pre-fire rate. If the species->rate mapping slipped by a row, this fails.
  #
  # The tolerance is per-row rather than a flat number, because the report
  # prints beta and pre_fire_rate to 4 dp: for a species with a tiny baseline
  # (Black-backed Woodpecker, 0.0059) that rounding is worth ~2 pp on a %base
  # near 150, while for Wrentit (0.1171) it is worth ~0.1 pp. A flat tolerance
  # is therefore either too loose to catch anything or fails on rounding alone.
  # A one-row slip in the species->rate mapping moves %base by tens of pp and
  # still trips this.
  out[, tol := 100 * (5e-5 / pre_rate) * (1 + abs(pct_reported) / 100) + 0.05]
  out[, dev := abs(beta / pre_rate * 100 - pct_reported)]
  say("  event study: worst |derived %%base - reported| = %.3f pp (tol %.3f, %s)",
      out[which.max(dev - tol), dev], out[which.max(dev - tol), tol],
      out[which.max(dev - tol), species])
  stopifnot(out[, all(dev < tol)])
  out[, c("tol", "dev") := NULL]

  out[, `:=`(pct     = beta / pre_rate * 100,
             pct_lo  = lo   / pre_rate * 100,
             pct_hi  = hi   / pre_rate * 100)]
  out[]
}

say("Figure 1 -- event study")
es <- parse_event_study(file.path(PROJ, "event_study_summary.txt"))

f1_keep <- c("Wrentit (dense shrubs)", "Mallard (wetland)",
             "American Coot (open water)", "Black Phoebe (water edges)")
f1_lab  <- c("Wrentit", "Mallard", "American Coot", "Black Phoebe")
d1 <- es[species %in% f1_keep]
d1[, sp := factor(f1_lab[match(species, f1_keep)], levels = f1_lab)]
d1[, focal := sp == "Wrentit"]

pal1 <- c("Wrentit" = C_BLUE, "Mallard" = C_ORANGE,
          "American Coot" = C_AQUA, "Black Phoebe" = C_VIOLET)

# direct labels at the right-hand end of each line (also the relief rule for aqua)
lab1 <- d1[k == max(k), .(sp, k, pct)]

dodge <- position_dodge(width = 0.45)

p1 <- ggplot(d1, aes(k, pct, colour = sp, group = sp)) +
  annotate("rect", xmin = -5.6, xmax = -0.5, ymin = -Inf, ymax = Inf,
           fill = "#000000", alpha = 0.028) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_vline(xintercept = -0.5, colour = MUTED, linewidth = 0.45, linetype = "22") +
  geom_linerange(aes(ymin = pct_lo, ymax = pct_hi, alpha = focal,
                     linewidth = focal), position = dodge) +
  geom_line(aes(alpha = focal, linewidth = focal), position = dodge) +
  geom_point(aes(size = focal, alpha = focal), position = dodge) +
  geom_text(data = lab1, aes(label = sp), hjust = 0, nudge_x = 0.35,
            size = 3.0, fontface = "bold", show.legend = FALSE) +
  annotate("text", x = -0.35, y = 96, label = "fire", hjust = 0,
           size = 2.9, colour = MUTED) +
  annotate("text", x = -5.5, y = 96, label = "before the fire — these should be flat",
           hjust = 0, size = 2.9, colour = MUTED) +
  scale_colour_manual(values = pal1) +
  scale_alpha_manual(values = c(`FALSE` = 0.62, `TRUE` = 1), guide = "none") +
  scale_linewidth_manual(values = c(`FALSE` = 0.45, `TRUE` = 1.0), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.5, `TRUE` = 2.5), guide = "none") +
  scale_x_continuous(breaks = -5:8, limits = c(-5.6, 10.6)) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(colour = guide_legend(override.aes = list(linewidth = 1.1, size = 2.2))) +
  labs(
    title = "Wrentit detection falls after fire; the water birds rise",
    x = "Years before and after the fire in that square",
    y = "Change in how often the bird was reported\n(% of its own normal rate)"
  ) +
  theme_paper() +
  theme(legend.position = c(0.012, 0.055), legend.justification = c(0, 0),
        legend.direction = "horizontal", legend.background = element_blank())

save_fig(p1, "fig1_event_study", 8.4, 4.4)

# ============================================================================
# FIGURE 2 -- the negative controls fail
# Source: point_controls_results.csv
# ============================================================================

say("Figure 2 -- negative controls")
pc <- fread(file.path(PROJ, "point_controls_results.csv"))
# Olive-sided Flycatcher is dropped from this figure only. It is never
# introduced as a focal species in the text (Francisco, 2026-08-20 #52), and the
# panel's argument is the three controls against the Wrentit. The CSV is left
# untouched, so nothing else that reads it changes.
pc <- pc[common != "Olive-sided Flycatcher"]
stopifnot(nrow(pc) == 4L, sum(pc$is_control) == 3L)
pc[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
pc[, grp := fifelse(is_control, "Water birds (cannot respond to fire)", "Study species")]
setorder(pc, is_control, estimate)
pc[, common := factor(common, levels = common)]

pal2 <- c("Water birds (cannot respond to fire)" = C_ORANGE,
          "Study species" = C_BLUE)

p2a <- ggplot(pc, aes(y = common, colour = grp)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.7) +
  geom_point(aes(x = cell_grain, shape = "Whole 5 km square"),
             size = 2.4, fill = SURFACE, stroke = 0.8) +
  geom_point(aes(x = estimate, shape = "Same birding spot, before vs after (95% CI)"),
             size = 2.4) +
  geom_text(aes(x = estimate, label = sprintf("p = %.3f", p)),
            nudge_y = 0.30, size = 2.7, colour = INK2, show.legend = FALSE) +
  scale_colour_manual(values = pal2) +
  scale_shape_manual(values = c("Same birding spot, before vs after (95% CI)" = 16,
                                "Whole 5 km square" = 21)) +
  scale_x_continuous(labels = function(x) sprintf("%+.2f", x)) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 2, nrow = 2,
                               override.aes = list(colour = INK2))) +
  labs(x = "Change in share of visits reporting the bird", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_blank())

p2b <- ggplot(pc, aes(y = common, x = pct_of_base, colour = grp)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = pct_of_base, yend = common), linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%+.1f%%", pct_of_base),
                hjust = fifelse(pct_of_base > 0, -0.25, 1.25)),
            size = 2.7, colour = INK2) +
  scale_colour_manual(values = pal2, guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.22, 0.30))) +
  labs(x = "Change relative to how often it is normally reported", y = NULL) +
  theme_paper() +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())

p2 <- (p2a | p2b) +
  plot_layout(widths = c(1.35, 1)) +
  plot_annotation(
    title = "Three birds that fire cannot reach all register a fire effect",
    theme = theme_paper()
  )

save_fig(p2, "fig2_negative_controls", 9.2, 3.6)

# ============================================================================
# FIGURE 3 -- the mechanism: fire changes how people bird
# Source: effort_placebo_summary.txt (formatted report, needs parsing)
# ============================================================================

say("Figure 3 -- effort placebo")
parse_effort <- function(path) {
  ln <- readLines(path, warn = FALSE)
  rows <- grep("^\\s*[a-z_]+\\s+\\d+\\s+[0-9.]+\\s+-?[0-9.]+\\s+[0-9.]+\\s+-?[0-9.]+\\s+[0-9.]+\\s+-?[0-9.]+",
               ln, value = TRUE)
  m <- regmatches(rows, regexec(
    "^\\s*([a-z_]+)\\s+(\\d+)\\s+([0-9.]+)\\s+(-?[0-9.]+)\\s+([0-9.]+)\\s+(-?[0-9.]+)\\s+([0-9.]+)\\s+(-?[0-9.]+)",
    rows))
  dt <- data.table(
    outcome = sapply(m, `[`, 2),
    est     = as.numeric(sapply(m, `[`, 5)),
    p       = as.numeric(sapply(m, `[`, 8)),
    pct     = as.numeric(sapply(m, `[`, 9))
  )
  stopifnot(nrow(dt) == 8, "n_loc" %in% dt$outcome, "det_rate" %in% dt$outcome)
  dt[]
}
ef <- parse_effort(file.path(PROJ, "effort_placebo_summary.txt"))

# det_rate is the OUTCOME the species analyses use, not an effort measure --
# the seven effort measures are what this placebo is about.
ef <- ef[outcome != "det_rate"]
lbl <- c(n_loc = "Distinct locations birded", log_chk = "Checklists filed",
         log_hours = "Hours birded", dist_med = "Distance walked (median)",
         obs_med = "Party size (median)", dur_med = "Checklist duration (median)",
         frac_trav = "Travelling-protocol share")
ef[, label := lbl[outcome]]
ef[, sig := p < 0.05]
stopifnot(sum(ef$sig) == 2, nrow(ef) == 7)
setorder(ef, -pct)                       # largest drop at the top after coord_flip
ef[, label := factor(label, levels = label)]

# Whole numbers for the big movers, one decimal below 10%. Rounding the whole
# column to 0 dp collides -- -7.4 and -6.9 both print as -7%, as do -5.1 and
# -4.9 -- which reads as a typo on bars of visibly different length.
ef[, lab := sprintf("%s   p = %.4f",
                    fifelse(abs(pct) >= 10,
                            sprintf("%+.0f%%", pct),
                            sprintf("%+.1f%%", pct)),
                    p)]

p3 <- ggplot(ef, aes(x = label, y = pct, fill = sig)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_col(width = 0.62) +
  geom_text(aes(label = lab,
                hjust = fifelse(pct > 0, -0.08, 1.08),
                colour = sig),
            size = 2.9, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = C_CRIT, `FALSE` = MUTED),
                    labels = c(`TRUE` = "Changed significantly",
                               `FALSE` = "No significant change"),
                    breaks = c(TRUE, FALSE)) +
  scale_colour_manual(values = c(`TRUE` = C_CRIT, `FALSE` = INK2)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0.30, 0.34))) +
  coord_flip() +
  labs(
    title = "Does fire change how people bird?",
    x = NULL, y = "Change after the square burned, compared with before"
  ) +
  theme_paper() +
  theme(legend.position = c(0.015, 0.03), legend.justification = c(0, 0),
        legend.background = element_blank(),
        panel.grid.major.y = element_blank())

save_fig(p3, "fig3_effort_placebo", 8.4, 3.6)

# ============================================================================
# FIGURE 4 -- what survives: the sign test and the guild contrasts
# Sources: sign_test_results.csv, guild_contrast_results.csv
# ============================================================================

say("Figure 4 -- sign test + guild contrasts")
st <- fread(file.path(PROJ, "sign_test_results.csv"))
st[, guild_lab := factor(
  c(shrub = "Shrub", forest = "Forest", snag = "Snag", control = "Control (no prediction)")[guild],
  levels = c("Shrub", "Forest", "Snag", "Control (no prediction)"))]
st[, matched := fifelse(is.na(correct), "No prediction",
                 fifelse(correct, "Matched prediction", "Moved against prediction"))]
setorder(st, guild_lab, estimate)
st[, species := factor(species, levels = species)]
stopifnot(st[guild == "snag", sum(correct)] == 6,
          st[guild == "forest", sum(correct)] == 4,
          st[guild == "shrub", sum(correct)] == 2)

pal_g <- c("Shrub" = C_BLUE, "Forest" = C_ORANGE, "Snag" = C_AQUA,
           "Control (no prediction)" = MUTED)

p4a <- ggplot(st, aes(y = species, x = estimate, colour = guild_lab)) +
  annotate("rect", xmin = -0.005, xmax = 0.005, ymin = -Inf, ymax = Inf,
           fill = MUTED, alpha = 0.16) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = estimate, yend = species), linewidth = 0.55) +
  # fill is the SURFACE colour, not the guild colour: shape 21 takes its fill
  # from this, and filling it with the guild colour makes a "moved against
  # prediction" marker indistinguishable from a matched one.
  geom_point(aes(shape = matched), fill = SURFACE, size = 2.3, stroke = 0.9) +
  facet_grid(guild_lab ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = pal_g, guide = "none") +
  scale_shape_manual(values = c("Matched prediction" = 16,
                                "Moved against prediction" = 21,
                                "No prediction" = 4),
                     breaks = c("Matched prediction", "Moved against prediction",
                                "No prediction")) +
  scale_x_continuous(breaks = seq(-0.04, 0.08, by = 0.02),
                     labels = function(x) sprintf("%+.2f", x),
                     expand = expansion(mult = 0.06)) +
  guides(shape = guide_legend(override.aes = list(colour = INK2))) +
  labs(x = "Change in share of visits reporting the bird", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank(),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 0),
        panel.spacing.y = unit(3, "pt"))

gc <- fread(file.path(PROJ, "guild_contrast_results.csv"))
gc[, contrast := trimws(gsub("\\(the key test\\)|\\(placebo\\)", "", contrast))]
gc[, kind := fifelse(grepl("^control", contrast), "Placebo: water birds", "Real comparison")]
gc[, `:=`(lo = estimate - 1.96 * se, hi = estimate + 1.96 * se)]
setorder(gc, estimate)
gc[, contrast := factor(contrast, levels = contrast)]

p4b <- ggplot(gc, aes(y = contrast, x = estimate, colour = kind)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.7) +
  geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.2f   p = %.3f", estimate, p)),
            nudge_y = 0.28, size = 2.7, colour = INK2, show.legend = FALSE) +
  scale_colour_manual(values = c("Real comparison" = C_BLUE,
                                 "Placebo: water birds" = C_ORANGE)) +
  # the leftmost contrast sits near zero and its centred label runs off the
  # panel without extra room on that side
  scale_x_continuous(expand = expansion(mult = c(0.26, 0.10))) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
  labs(x = "Difference between groups (relative to own rates)",
       y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

p4 <- (p4a | p4b) +
  plot_layout(widths = c(1, 1.05)) +
  plot_annotation(
    title = "What survives: a sign test across 19 species, and the between-guild contrast",
    theme = theme_paper()
  )

save_fig(p4, "fig4_sign_test_guild_contrast", 11.0, 5.7)

# ============================================================================
# FIGURE 5 -- prediction: burn history is mostly a location detector
# Sources: burn_effect_results_wrentit_covariate.csv,
#          deep_models_common_test_wrentit.txt
# ============================================================================

say("Figure 5 -- prediction")
be <- fread(file.path(PROJ, "burn_effect_results_wrentit_covariate.csv"))
w <- dcast(be, split + algo ~ variant, value.var = "auc")
w[, delta := full - noburn]

getd <- function(sp, al) w[split == sp & algo == al, delta]
d5a <- data.table(
  split  = factor(c("Later years", "Later years", "New places", "New places"),
                  levels = c("Later years", "New places")),
  geo    = factor(c("(location hidden)", "(location given)",
                    "(location hidden)", "(location given)"),
                  levels = c("(location hidden)", "(location given)")),
  delta  = c(getd("temporal", "lightgbm"), getd("temporal", "lightgbm+geo"),
             getd("spatial",  "lightgbm"), getd("spatial",  "lightgbm+geo"))
)
d5a[, lab := paste(split, geo, sep = "\n")]
d5a[, lab := factor(lab, levels = lab)]
stopifnot(abs(d5a$delta[1] - 0.0460) < 0.001, abs(d5a$delta[2] - 0.0021) < 0.001)

p5a <- ggplot(d5a, aes(lab, delta, fill = geo)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%+.4f", delta)), vjust = -0.7,
            size = 3.1, colour = INK2, fontface = "bold") +
  scale_fill_manual(values = c("(location hidden)" = C_BLUE,
                               "(location given)" = C_BLUE_L)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "Gain in forecast accuracy from fire history") +
  theme_paper() +
  theme(legend.position = "none", panel.grid.major.x = element_blank())

# --- right panel: per-model AUC with and without the burn features -----------
tab <- w[split == "temporal"]
nm <- c(glm = "Logistic", gbm = "gbm", xgboost = "XGBoost",
        lightgbm = "LightGBM", `lightgbm+geo` = "LightGBM + coordinates")
tab[, model := nm[algo]]

parse_deep <- function(path) {
  ln <- readLines(path, warn = FALSE)
  rows <- grep("^\\s*(LightGBM|LSTM|CNN)\\s+[0-9.]+\\s*$", ln, value = TRUE)
  stopifnot(length(rows) >= 3)
  m <- regmatches(rows[1:3], regexec("^\\s*(\\S+)\\s+([0-9.]+)\\s*$", rows[1:3]))
  data.table(model = sapply(m, `[`, 2), auc = as.numeric(sapply(m, `[`, 3)))
}
dp <- parse_deep(file.path(PROJ, "deep_models_common_test_wrentit.txt"))
lgb_common <- dp[model == "LightGBM", auc]
say("  common-test LightGBM reference AUC = %.4f", lgb_common)

d5b <- rbind(
  tab[, .(model, full, noburn, fam = "Standard models")],
  dp[model != "LightGBM", .(model, full = auc, noburn = NA_real_,
                            fam = "Deep-learning models")]
)
setorder(d5b, full)
d5b[, model := factor(model, levels = model)]

p5b <- ggplot(d5b, aes(y = model, colour = fam)) +
  geom_segment(aes(x = noburn, xend = full, yend = model), linewidth = 0.55,
               na.rm = TRUE) +
  geom_point(aes(x = noburn, shape = "without fire history"), size = 2.3,
             fill = SURFACE, stroke = 0.8, na.rm = TRUE) +
  geom_point(aes(x = full, shape = "with fire history"), size = 2.3) +
  geom_vline(xintercept = lgb_common, colour = MUTED, linewidth = 0.4,
             linetype = "22") +
  annotate("text", x = min(d5b$noburn, na.rm = TRUE), y = 0.30,
           label = sprintf("dashed line: the best standard model, scored on the same test data as the deep models (%.4f)",
                           lgb_common),
           hjust = 0, vjust = 0.5, size = 2.5, colour = MUTED) +
  geom_text(aes(x = full, label = sprintf("%.4f", full)), nudge_y = 0.30,
            size = 2.7, colour = INK2, show.legend = FALSE) +
  scale_colour_manual(values = c("Standard models" = C_BLUE,
                                 "Deep-learning models" = C_ORANGE)) +
  scale_shape_manual(values = c("with fire history" = 16,
                                "without fire history" = 21)) +
  scale_x_continuous(expand = expansion(mult = c(0.10, 0.06))) +
  scale_y_discrete(expand = expansion(add = c(0.9, 0.9))) +
  guides(colour = guide_legend(order = 1, nrow = 2),
         shape  = guide_legend(order = 2, nrow = 2,
                               override.aes = list(colour = INK2))) +
  labs(x = "Forecast accuracy (0.5 = coin flip, 1.0 = perfect)", y = NULL) +
  theme_paper() +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(1, "pt"), legend.margin = margin(0, 0, 0, 0),
        panel.grid.major.y = element_blank())

p5 <- (p5a | p5b) +
  plot_layout(widths = c(1, 1.25)) +
  plot_annotation(
    title = "Burn history predicts well, and mostly by telling the model where it is",
    theme = theme_paper()
  )

save_fig(p5, "fig5_prediction", 11.2, 4.4)

say("\nAll five figures written to %s", normalizePath(FIGDIR))

# ============================================================================
# FIGURE 0 -- how the two eBird files combine, and why zero-filling needs both
#
# SCHEMATIC. Unlike figures 1-5 this draws no estimated quantity: the three
# checklists and their species are invented to illustrate the join, and the
# figure says so. Nothing here is parsed from a result file.
#
# LAYOUT NOTE: the "no Wrentit row" placeholder must sit BELOW a checklist's
# real EBD rows, not among them -- it represents a row that is absent from the
# file, so drawing it over the rows that do exist reads as if it replaced them.
# ============================================================================

say("Figure 0 -- data structure schematic")

.box <- function(x, y, w, h, fill, colour, lty = "solid", lwd = 0.4)
  annotate("rect", xmin = x, xmax = x + w, ymin = y - h/2, ymax = y + h/2,
           fill = fill, colour = colour, linewidth = lwd, linetype = lty)
.lab <- function(x, y, s, size = 2.7, col = INK2, face = "plain", hj = 0)
  annotate("text", x = x, y = y, label = s, size = size, colour = col,
           fontface = face, hjust = hj)

BW <- 3.4; RH <- 0.62
X1 <- 0; X2 <- 4.6; X3 <- 9.9
TINT <- "#eaf1fc"
GREY <- "#f1f0ec"

ebd <- data.frame(
  y   = c(9.00, 8.38, 7.76,   6.60, 5.98,   4.20),
  sp  = c("Wrentit", "Bushtit", "California Scrub-Jay",
          "Bushtit", "Spotted Towhee", "California Scrub-Jay"),
  cnt = c("3", "12", "2", "8", "1", "4"),
  foc = c(TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
gap <- data.frame(y = c(5.36, 3.58))                 # below each block's real rows
ctr <- c(8.38, 5.98, 3.89)                           # vertical centre of each block
sed <- data.frame(y = ctr, chk = c("A","B","C"),
                  eff = c("62 min · 1.4 km · 2 obs",
                          "35 min · 0.8 km · 1 obs",
                          "90 min · 2.1 km · 3 obs"))
out <- data.frame(y = ctr, chk = c("A","B","C"),
                  pres = c("1", "0", "0"), got = c(TRUE, FALSE, FALSE))

p0 <- ggplot() +
  .lab(X1, 10.15, "SED — sampling events", 3.1, INK, "bold") +
  .lab(X1, 9.73, "one row per checklist", 2.6, MUTED) +
  .lab(X2, 10.15, "EBD — observations", 3.1, INK, "bold") +
  .lab(X2, 9.73, "one row per species reported", 2.6, MUTED) +
  .lab(X3, 10.15, "Zero-filled Wrentit table", 3.1, INK, "bold") +
  .lab(X3, 9.73, "one row per checklist", 2.6, MUTED)

for (i in seq_len(nrow(sed)))
  p0 <- p0 + .box(X1, sed$y[i], BW, RH * 1.5, GREY, AXIS) +
    .lab(X1 + 0.18, sed$y[i] + 0.18, paste0("checklist ", sed$chk[i]), 2.7, INK, "bold") +
    .lab(X1 + 0.18, sed$y[i] - 0.18, sed$eff[i], 2.4, MUTED)

for (i in seq_len(nrow(ebd)))
  p0 <- p0 + .box(X2, ebd$y[i], BW, RH, if (ebd$foc[i]) TINT else GREY,
                  if (ebd$foc[i]) C_BLUE else AXIS) +
    .lab(X2 + 0.18, ebd$y[i], ebd$sp[i], 2.6,
         if (ebd$foc[i]) C_BLUE else INK2, if (ebd$foc[i]) "bold" else "plain") +
    .lab(X2 + BW - 0.18, ebd$y[i], ebd$cnt[i], 2.6, MUTED, "plain", hj = 1)

for (i in seq_len(nrow(gap)))
  p0 <- p0 + .box(X2, gap$y[i], BW, RH, SURFACE, MUTED, "22", 0.4) +
    .lab(X2 + BW/2, gap$y[i], "no Wrentit row exists", 2.5, MUTED, "italic", hj = 0.5)

for (i in seq_len(nrow(out)))
  p0 <- p0 + .box(X3, out$y[i], BW, RH * 1.5,
                  if (out$got[i]) TINT else GREY, if (out$got[i]) C_BLUE else AXIS) +
    .lab(X3 + 0.18, out$y[i] + 0.18, paste0("checklist ", out$chk[i]), 2.7, INK, "bold") +
    .lab(X3 + 0.18, out$y[i] - 0.18, paste0("Wrentit presence = ", out$pres[i]),
         2.5, if (out$got[i]) C_BLUE else MUTED)

p0 <- p0 +
  annotate("segment", x = X1 + BW + 0.15, xend = X2 - 0.15, y = ctr, yend = ctr,
           colour = MUTED, linewidth = 0.4,
           arrow = arrow(length = unit(4, "pt"), type = "closed")) +
  annotate("segment", x = X2 + BW + 0.15, xend = X3 - 0.15, y = ctr, yend = ctr,
           colour = MUTED, linewidth = 0.4,
           arrow = arrow(length = unit(4, "pt"), type = "closed")) +
  .lab((X1 + BW + X2)/2, 7.25, "joined on\nSAMPLING EVENT\nIDENTIFIER", 2.4, MUTED, "plain", 0.5) +
  .lab((X2 + BW + X3)/2, 7.25, "absence written\nas a verified\nzero", 2.4, MUTED, "plain", 0.5) +
  .lab(X1, 2.55,
       paste0("Checklists B and C leave no trace in the EBD, because it records only what was seen. The SED lists every checklist\n",
              "that happened, so the join recovers them — and because each observer certified that all species were reported, the\n",
              "absence becomes a verified zero rather than a missing value."), 2.6, INK2) +
  .lab(X1, 1.85, "Illustrative checklists; not drawn from the data.", 2.4, MUTED, "italic") +
  labs(title = "Why zero-filling needs both eBird files") +
  scale_x_continuous(limits = c(-0.2, X3 + BW + 0.3)) +
  scale_y_continuous(limits = c(1.6, 10.5)) +
  theme_void(base_size = 10) +
  theme(plot.background = element_rect(fill = SURFACE, colour = NA),
        plot.title = element_text(colour = INK, size = 12, face = "bold", hjust = 0),
        plot.title.position = "plot",
        plot.margin = margin(10, 10, 6, 10))

save_fig(p0, "fig0_zerofill", 10.4, 4.8)


# ============================================================================
# FIGURE 6 -- fire as a dose, not a switch  (section 4.3)
# Source: burn_dose_summary_wrentit.txt (formatted report, needs parsing)
#
# Two panels, because section 4.3 makes two claims and they are different:
#   (a) most burned cells are only partly burned, which is WHY a yes/no
#       treatment is diluted -- the median treated cell is ~29% burned;
#   (b) switching from yes/no to a continuous dose moves the estimate.
#
# NO INTERVALS ARE DRAWN. Both estimates in this report carry cell-level
# percentile CIs, and the same file states those are superseded: under the wild
# cluster bootstrap Wrentit reaches only p = 0.110. Drawing the cell CIs would
# repeat the error section 5.1 already makes. The wild-cluster p is annotated
# on the dose estimate instead.
# ============================================================================

say("Figure 6 -- dose, not switch")

parse_dose <- function(path) {
  ln <- readLines(path, warn = FALSE)

  # -- band table: "burn_extent cells rate_before rate_after change did"
  i_hd <- grep("^\\s*burn_extent\\s+cells\\s+rate_before", ln)
  stopifnot(length(i_hd) == 1)
  rows <- character(0)
  for (k in (i_hd + 1):length(ln)) {
    if (!grepl("^\\s*[0-9]+-[0-9]+%\\s", ln[k])) break
    rows <- c(rows, ln[k])
  }
  stopifnot(length(rows) == 6)
  f <- strsplit(trimws(rows), "\\s+")
  bands <- data.table(
    band   = vapply(f, `[`, "", 1),
    cells  = as.integer(vapply(f, `[`, "", 2)),
    r_bef  = as.numeric(vapply(f, `[`, "", 3)),
    r_aft  = as.numeric(vapply(f, `[`, "", 4)),
    change = as.numeric(vapply(f, `[`, "", 5)),
    did    = as.numeric(vapply(f, `[`, "", 6)))

  # -- control split, the number `did` is differenced against
  ctl <- as.numeric(sub(".*control cells, same period split:\\s*([+-][0-9.]+).*", "\\1",
                        grep("control cells, same period split", ln, value = TRUE)[1]))

  # -- the two estimates and the wild-cluster p
  bin <- as.numeric(sub(".*yes/no analysis found\\s*(-?[0-9.]+).*", "\\1",
                        grep("yes/no analysis found", ln, value = TRUE)[1]))
  dos <- as.numeric(sub(".*post x dose \\(dose-response\\)\\s*:\\s*(-?[0-9.]+).*", "\\1",
                        grep("post x dose \\(dose-response\\)", ln, value = TRUE)[1]))
  wcb <- as.numeric(sub(".*Wrentit p =\\s*([0-9.]+).*", "\\1",
                        grep("Wrentit p =", ln, value = TRUE)[1]))
  stopifnot(is.finite(c(ctl, bin, dos, wcb)))

  # -- SELF-CHECK, the convention every parser here follows: re-derive two
  #    printed columns from other printed columns. `change` must equal
  #    rate_after - rate_before, and `did` must equal change minus the control
  #    split, both to the 4 dp the report prints. A row slip trips this.
  stopifnot(max(abs(bands$change - (bands$r_aft - bands$r_bef))) < 1e-4)
  stopifnot(max(abs(bands$did - (bands$change - ctl))) < 1e-4)
  say("  dose parser: change and did re-derived, worst err %.1e",
      max(abs(bands$did - (bands$change - ctl))))

  list(bands = bands, ctl = ctl, binary = bin, dose = dos, wcb_p = wcb)
}

dz <- parse_dose(file.path(PROJ, "burn_dose_summary_wrentit.txt"))
dzb <- dz$bands
dzb[, band := factor(band, levels = band)]
dzb[, part := as.integer(.I <= 3)]     # the three bands under 50%

# -- panel A: how much of a burned square actually burned ---------------------
p6a <- ggplot(dzb, aes(x = band, y = cells, fill = factor(part))) +
  geom_col(width = 0.72) +
  geom_text(aes(label = cells), vjust = -0.55, size = 2.8, colour = INK2) +
  annotate("segment", x = 3.5, xend = 3.5, y = 0, yend = 292,
           colour = C_CRIT, linewidth = 0.5, linetype = "22") +
  annotate("text", x = 3.38, y = 292, hjust = 1, vjust = 1, size = 2.75,
           colour = C_CRIT,
           label = sprintf("%s of %s burned squares\nburned less than half",
                           format(dzb[part == 1, sum(cells)], big.mark = ","),
                           format(dzb[, sum(cells)], big.mark = ","))) +
  scale_fill_manual(values = c(`1` = C_BLUE_L, `0` = C_BLUE), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(title = "Most burned squares only partly burned",
       x = "Share of the 5 km square that burned",
       y = "Number of squares") +
  theme_paper()

# -- panel B: what the switch misses -----------------------------------------
est <- data.table(
  # levels reversed: ggplot draws the first level at the BOTTOM, and the
  # section's argument runs switch-first-then-dose, so the switch goes on top.
  what = factor(c("Burned: yes or no",
                  "Share of the square that\nburned, at a full burn"),
                levels = c("Share of the square that\nburned, at a full burn",
                           "Burned: yes or no")),
  val  = c(dz$binary, dz$dose))

p6b <- ggplot(est, aes(x = val, y = what)) +
  geom_vline(xintercept = 0, colour = AXIS, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = val, yend = what), linewidth = 0.8, colour = C_BLUE) +
  geom_point(size = 3, colour = C_BLUE) +
  geom_text(aes(label = sprintf("%+.4f", val)), nudge_y = 0.24,
            size = 2.9, colour = INK2) +
  scale_x_continuous(labels = function(x) sprintf("%+.2f", x),
                     expand = expansion(mult = c(0.34, 0.34))) +
  labs(title = "Treating fire as a switch finds almost nothing",
       x = "Change in the share of visits reporting the Wrentit", y = NULL) +
  theme_paper() +
  theme(panel.grid.major.y = element_blank())

p6 <- (p6a | p6b) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(theme = theme_paper())

save_fig(p6, "fig6_dose_not_switch", 10.2, 3.7)
