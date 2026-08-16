# =============================================================================
# THE SIGN TEST -- the project's actual falsification design.
#
# Each species' predicted direction is committed to species_registry.csv BEFORE
# the species is run. This checks how many came out that way.
#
# WHY IT MATTERS MORE THAN ANY p-VALUE HERE
#   The binding constraint on every individual estimate is 225 independent
#   fires. The sign test does not depend on that: it scales with the number of
#   SPECIES, so it sidesteps the power problem instead of fighting it.
#
# WHY THE PATTERN LOOKS PARANOID
#   Two parsing bugs produced wrong counts before this was written:
#     1. grep("post x dose") matched the MODEL FORMULA line
#        ("rate ~ cell FE + year FE + post + post x dose"), which has no number,
#        so every estimate came back NA. verify_published.R already carries a
#        comment warning about exactly this. Require a colon AND a signed number.
#     2. The label CHANGED when dose_burn_effect.R became registry-driven:
#        pre-2026-07-30 summaries say "post x fraction", later ones "post x dose".
#        Matching only the new label silently dropped Brown Creeper and OSFL --
#        two species with results, missing from the count, no error.
#   Hence: match "dose-response" plus a signed number, which is unique to the
#   results line under either label and needs no escaped parentheses.
# =============================================================================
suppressMessages({ library(data.table) })

proj <- Sys.getenv("EBIRD_PROJ_ROOT", "..")
dir_ <- file.path(proj, "eBird")
source(file.path(dir_, "project_config.R"))
out_sum <- file.path(dir_, "sign_test_summary.txt")

s <- cfg_species(active_only = FALSE)
r <- rbindlist(lapply(seq_len(nrow(s)), function(i) {
  f <- file.path(dir_, sprintf("burn_dose_summary_%s.txt", s$prefix[i]))
  if (!file.exists(f)) return(NULL)
  ln <- grep("dose-response.*: *[-+][0-9]", readLines(f), value = TRUE)[1]
  if (is.na(ln)) return(NULL)
  v <- as.numeric(regmatches(ln, gregexpr("[-+][0-9]*[.][0-9]+", ln))[[1]][1])
  data.table(species = s$common_name[i], prefix = s$prefix[i], guild = s$guild[i],
             predicted = s$predicted[i], estimate = v)
}))
stopifnot(nrow(r) > 0, !any(is.na(r$estimate)))   # an NA here means the parse broke again

r[, observed := fifelse(estimate < 0, "negative", "positive")]
r[, correct  := fifelse(predicted == "none", NA, predicted == observed)]
fwrite(r, file.path(dir_, "sign_test_results.csv"))

t  <- r[!is.na(correct)]
bt <- binom.test(sum(t$correct), nrow(t), 0.5, alternative = "greater")

sink(out_sum)
cat("THE SIGN TEST -- predicted directions vs observed\n")
cat("================================================\n")
cat("Generated : ", format(Sys.time()), "\n", sep = "")
cat("Measure   : continuous burn-extent dose (frac_max), effect at 100% burned\n")
cat("Directions were recorded in species_registry.csv BEFORE each species ran.\n\n")

print(as.data.frame(r[order(guild, -estimate),
      .(species, guild, predicted, estimate = round(estimate, 4), observed,
        hit = fifelse(is.na(correct), "-- control --",
              fifelse(correct, "YES", "no")))]), row.names = FALSE)

cat(sprintf("\nOVERALL: %d of %d correct, one-tailed p = %.4f\n",
            sum(t$correct), nrow(t), bt$p.value))

cat("\n-- BY GUILD (this is where the story is) ----------------------------\n\n")
g <- t[, .(n = .N, correct = sum(correct)), by = guild]
g[, p := vapply(seq_len(.N), function(i)
       binom.test(correct[i], n[i], 0.5, alternative = "greater")$p.value, numeric(1))]
print(as.data.frame(g[order(p), .(guild, n, correct, p = round(p, 4))]), row.names = FALSE)

cat("\n-- NEGATIVE CONTROLS ------------------------------------------------\n")
cat("  Wetland and open-water birds. Fire should do NOTHING to them. A clear\n")
cat("  effect here is evidence the DESIGN is wrong, not that fire moves ducks.\n\n")
ctl <- r[predicted == "none"]
if (!nrow(ctl)) cat("  none have a result yet\n") else {
  print(as.data.frame(ctl[, .(species, estimate = round(estimate, 4))]), row.names = FALSE)
  cat(sprintf("\n  mean |effect| across controls : %.4f\n", mean(abs(ctl$estimate))))
  cat(sprintf("  mean |effect| across predicted: %.4f\n", mean(abs(t$estimate))))
  cat("  The controls should be the SMALLER of those two.\n")
}

cat("\n-- HOW TO READ ------------------------------------------------------\n")
cat("  This is a test of the PATTERN, not of any one species. It is the only\n")
cat("  claim in the project that does not depend on the 225-fire ceiling.\n")
cat("  Read the per-guild rows before the overall row: a guild that is right\n")
cat("  6 of 6 and one that is right 2 of 6 average out to something\n")
cat("  uninformative, and the average is the least interesting number here.\n")
sink()
cat(readLines(out_sum), sep = "\n")
