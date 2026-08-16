# =============================================================================
# DOES SEVERITY AFFECT DIFFERENT SPECIES DIFFERENTLY?  -- tested as a CONTRAST.
#
# WHY A CONTRAST AND NOT A TABLE OF SEPARATE ESTIMATES
#   Nineteen per-species estimates do not test whether species DIFFER. Reading a
#   table and observing "snag is positive, forest is negative" is eyeballing, not
#   inference -- especially when most individual estimates are not significant.
#
#   Worse, the per-species estimates are NOT independent. Every species is
#   measured on the SAME cells, the SAME fires and the SAME checklists; only the
#   bird being counted changes. Their errors are correlated, so neither a naive
#   Cochran's Q nor se1^2 + se2^2 is valid.
#
#   That correlation is the opportunity. The negative controls showed a shared
#   observer artifact -- after a fire, birders visit 43% fewer spots, which moves
#   EVERY species' detection the same way. In a DIFFERENCE between two guilds
#   measured on identical cell-years, that shared component cancels exactly.
#
# THE DESIGN
#   Outcome = (mean normalised detection of snag species)
#           - (mean normalised detection of forest species), per cell-year.
#   Each species is first divided by its own pre-fire treated mean, so guilds
#   with different baseline abundances are on a common scale and the contrast is
#   a proportional one. Then the usual within-cell model:
#
#     y_it = a_i + b_t + p*post + g*(post x extent) + d*(post x extent x mean-sev)
#
#   `d` IS THE ANSWER. It is the difference in the severity response between the
#   two guilds. Clustered by FIRE, wild cluster bootstrap.
#
#   The control guild is run the same way as a placebo: contrasting the controls
#   against themselves-minus-forest should show nothing beyond what the shared
#   artifact leaves.
# =============================================================================
suppressMessages({ library(data.table) })
setDTthreads(0)
say <- function(...) { message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"),
                                       paste0(...))); flush.console() }

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "panel_utils.R"))
source(file.path(dir_, "infer.R"))
B_REPS  <- as.integer(Sys.getenv("GC_B", "999"))
out_sum <- file.path(dir_, "guild_contrast_summary.txt")

sp  <- cfg_species(active_only = FALSE)
sp  <- sp[file.exists(file.path(dir_, sprintf("%s_cells_5km_16d.csv", prefix)))]
ref <- sp$prefix[1]        # computed BEFORE any drop, so the treatment
                           # assignment and severity cache never move

# GC_DROP: comma-separated prefixes to exclude, for robustness diagnostics.
# Unset (the default) reproduces the published run exactly; when set, output
# goes to *_<GC_TAG> files so a diagnostic can never overwrite a published one.
#
# Why this exists: control - snag came back at -2.19 (p 0.001), which is
# consistent BOTH with snag genuinely responding to severity and with the snag
# guild simply being noisier than the others. Each species is divided by its own
# pre-fire mean, and three snag species have very small ones -- Black-backed
# Woodpecker 0.0059, Mountain Bluebird 0.0065, Lazuli Bunting 0.0223 -- so their
# normalised series swing hard. Dropping those three tests which story holds.
drop_sp <- trimws(strsplit(Sys.getenv("GC_DROP", ""), ",")[[1]])
drop_sp <- drop_sp[nzchar(drop_sp)]
tag     <- Sys.getenv("GC_TAG", if (length(drop_sp)) "subset" else "")
if (length(drop_sp)) {
  miss <- setdiff(drop_sp, sp$prefix)
  if (length(miss)) stop("GC_DROP names species not present: ",
                         paste(miss, collapse = ", "))
  sp <- sp[!prefix %in% drop_sp]
  say("GC_DROP: excluded ", paste(drop_sp, collapse = ", "),
      " -> ", nrow(sp), " species remain")
}
suffix  <- if (nzchar(tag)) paste0("_", tag) else ""
out_sum <- file.path(dir_, paste0("guild_contrast_summary", suffix, ".txt"))

# severity composition + treatment, from the reference numbering. Geography is
# shared, so one treatment assignment serves every species; joining per species
# on (gx,gy) keeps it safe across numberings.
sev <- .cache_sev(ref)
ct  <- cell_treatment(ref, dose = "frac_max")
map <- unique(fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", ref)),
                    select = c("cell_id","gx","gy"), showProgress = FALSE))
base <- merge(merge(ct[, .(cell_id, treat, ev_year)], sev, by = "cell_id", all.x = TRUE),
              map, by = "cell_id")
for (cc in c("s1","s2","s3","s4","frac_tot","mean_sev"))
  set(base, which(is.na(base[[cc]])), cc, 0)

say("normalised detection per species")
L <- lapply(seq_len(nrow(sp)), function(i) {
  d <- fread(file.path(dir_, sprintf("%s_cells_5km_16d.csv", sp$prefix[i])),
             select = c("gx","gy","year","n_checklists","n_detections"), showProgress = FALSE)
  d <- d[, .(det = sum(n_detections), chk = sum(n_checklists)), by = .(gx, gy, year)][chk > 0]
  d <- merge(d, base[, .(gx, gy, treat, ev_year)], by = c("gx","gy"))
  d[, post := as.numeric(!is.na(ev_year) & year >= ev_year)]
  b <- d[treat == 1 & post == 0, sum(det)/sum(chk)]     # this species' own baseline
  if (!is.finite(b) || b <= 0) return(NULL)
  d[, .(gx, gy, year, chk, yn = (det/chk) / b, guild = sp$guild[i])]
})
A <- rbindlist(L[!vapply(L, is.null, logical(1))])
say("  ", uniqueN(A$guild), " guilds, ", format(nrow(A), big.mark=","), " species-cell-years")

# guild means per cell-year, then the contrasts
G <- dcast(A[, .(y = mean(yn), chk = chk[1]), by = .(gx, gy, year, guild)],
           gx + gy + year ~ guild, value.var = c("y","chk"))
setnames(G, sub("^y_", "", names(G)))

fit_contrast <- function(G, a, b, label) {
  if (!all(c(a, b) %in% names(G))) return(NULL)
  d <- G[is.finite(get(a)) & is.finite(get(b))]
  d[, yv := get(a) - get(b)]
  d <- merge(d, base[, .(gx, gy, cell_id, treat, ev_year, frac_tot, mean_sev)],
             by = c("gx","gy"))
  d[, post := as.numeric(!is.na(ev_year) & year >= ev_year)]
  d[, `:=`(rate = yv, w = 1,                       # unweighted: yv is already a ratio
           p_frac = post * frac_tot,
           p_msev = post * frac_tot * mean_sev)]
  d <- attach_fire(d[, .(cell_id, gx, gy, year, rate, w, post, p_frac, p_msev, treat)], ref)
  P  <- .prep(d, c("post","p_frac","p_msev"))
  bb <- .fit(P, P$y)
  V  <- crve(P, as.vector(P$y - P$X %*% bb))
  data.table(contrast = label, estimate = bb[3], se = sqrt(V[3,3]),
             t = bb[3]/sqrt(V[3,3]),
             p = wcb_p(P, 3L, bb[3], sqrt(V[3,3]), B = B_REPS),
             n_fires = uniqueN(d[treat == 1]$fire))
}

say("fitting contrasts")
res <- rbindlist(list(
  fit_contrast(G, "snag",   "forest",  "snag - forest      (the key test)"),
  fit_contrast(G, "snag",   "shrub",   "snag - shrub"),
  fit_contrast(G, "shrub",  "forest",  "shrub - forest"),
  fit_contrast(G, "control","forest",  "control - forest   (placebo)"),
  fit_contrast(G, "control","shrub",   "control - shrub    (placebo)"),
  # Added 2026-08-13 at Francisco's request. snag is the guild carrying the
  # surviving claim, so this is the placebo that most needed running: it asks
  # whether the contrast machinery invents a snag-sized difference between two
  # sets of species that severity cannot separate. APPENDED, not inserted, so
  # the row positions of the five published contrasts do not move -- a
  # positional check elsewhere would silently read the wrong row otherwise.
  fit_contrast(G, "control","snag",    "control - snag     (placebo)")
), use.names = TRUE)

fwrite(res, file.path(dir_, paste0("guild_contrast_results", suffix, ".csv")))

sink(out_sum)
cat("DOES SEVERITY AFFECT DIFFERENT GUILDS DIFFERENTLY?\n")
cat("==================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Outcome   : difference between two guilds' normalised detection, same\n")
cat("            cell-years. Each species divided by its own pre-fire mean.\n")
cat("Reported  : the post x extent x mean-severity coefficient = the DIFFERENCE\n")
cat("            in severity response between the two guilds.\n")
cat("Inference : cluster-robust by FIRE + wild cluster bootstrap, B = ", B_REPS, "\n\n", sep="")

cat("-- WHY A CONTRAST -----------------------------------------------------\n")
cat("  The negative controls exposed a shared observer artifact: after a fire,\n")
cat("  birders visit 43% fewer spots, which moves EVERY species the same way.\n")
cat("  In a difference between two guilds on identical cell-years, that shared\n")
cat("  component cancels. A contrast is therefore better identified than either\n")
cat("  guild's level -- which is exactly why the levels failed and this may not.\n\n")

print(as.data.frame(res[, .(contrast, est = round(estimate,4), se = round(se,4),
                            t = round(t,2), p = round(p,4),
                            sig = fifelse(p < 0.05, "*", ""), fires = n_fires)]),
      row.names = FALSE)

cat("\n-- HOW TO READ --------------------------------------------------------\n")
cat("  snag - forest significant  => severity moves these guilds DIFFERENTLY.\n")
cat("     That is the heterogeneity claim, tested directly rather than inferred\n")
cat("     from a table of separate estimates.\n")
cat("  The placebo rows involve the control guild, which should not respond to\n")
cat("  severity at all. If a placebo contrast is as large as snag - forest, the\n")
cat("  contrast is not cancelling what it is supposed to cancel.\n")
cat("\n-- LIMITS -------------------------------------------------------------\n")
cat("  * Cancellation is exact only if the artifact is common ACROSS GUILDS.\n")
cat("    If it scales with a species' abundance or detectability it will not\n")
cat("    fully cancel, and the normalisation only partly addresses that.\n")
cat("  * Still ", res$n_fires[1], " fires. This inherits the burned-cell ceiling;\n", sep="")
cat("    the distance treatment reaches 362-386 and has not been applied here.\n")
cat("  * Guild membership was assigned in advance, but the CONTRASTS to run were\n")
cat("    chosen after seeing the per-species table.\n")
sink()
cat(readLines(out_sum), sep = "\n")
say("DONE")
