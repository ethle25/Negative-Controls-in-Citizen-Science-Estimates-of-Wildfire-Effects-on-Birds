#!/bin/bash
# Unattended overnight queue, 2026-08-01.
#
# Ordered by VALUE, not by dependency, so that if it overruns the most useful
# results are already on disk. Every phase logs to its own file and a failure
# never blocks a later phase -- the point of running unattended is to wake up to
# evidence, including evidence of what broke.
#
# NOTHING HERE EDITS SHARED CODE. Phases 0 and C are read-only computations;
# A and B use the existing verified build path. That is deliberate: no published
# number can move while nobody is watching. Phase D re-runs the 51 checks so any
# drift is visible immediately.
#
#   phase 0  fire-count gate            ~10 min   Step 1's decision
#   phase A  zero-fill 14 new species   ~1 hr     Step 3 prerequisite
#   phase B  join + regrid those 14     ~2.8 hr   Step 3 prerequisite
#   phase C  point-level distances      ~1 hr     Step 2 prerequisite
#   phase D  full verify                ~8 min    safety
cd "$(dirname "$0")"
LOG=overnight_run.log
: > $LOG

stamp () { date "+%H:%M:%S"; }
run () {                      # run <label> <logfile> <command...>
  local label=$1 lf=$2; shift 2
  echo "[$(stamp)] START  $label" | tee -a $LOG
  local t0=$SECONDS
  if "$@" > "$lf" 2>&1; then
    echo "[$(stamp)] OK     $label  ($(( (SECONDS-t0)/60 )) min)" | tee -a $LOG
  else
    echo "[$(stamp)] FAILED $label  ($(( (SECONDS-t0)/60 )) min) -- see $lf" | tee -a $LOG
    tail -5 "$lf" | sed 's/^/           /' | tee -a $LOG
  fi
}

echo "=== overnight queue started $(date) ===" | tee -a $LOG
df -h . | tail -1 | sed 's/^/disk: /' | tee -a $LOG
echo | tee -a $LOG

# --- phase 0: the highest-value unknown, first --------------------------------
run "0a fire_distance_lookup (cell grain)" ov_0a_distance.log \
    Rscript fire_distance_lookup.R
run "0b step1_gate (does distance add fires?)" ov_0b_gate.log \
    Rscript step1_gate.R

# --- phase A: zero-fill every active registry species --------------------------
# The species set changed, so the 52 GB filter pass re-runs. Species already
# built are re-written; harmless, and simpler than special-casing them.
run "A  zerofill_from_registry (14 new species)" ov_A_zerofill.log \
    Rscript zerofill_from_registry.R

# --- phase B: build each new species. run_species.sh already continues on
# --- failure per species, and rare birds (Black-backed Woodpecker, 5,924
# --- records) are the ones most likely to fall over.
NEW="spotow caltow bewwre buggna bkbwoo haiwoo wesblu moublu lazbun pacwre1 herwar mallar3 y00475 blkpho"
for P in $NEW; do
  run "B  build $P" "ov_B_$P.log" env STAGES=build ./run_species.sh "$P"
done

# --- phase C: point-level distances (Step 2's prerequisite) -------------------
run "C  point_fire_distance (snapped to 500 m)" ov_C_pointdist.log \
    Rscript point_fire_distance.R

# --- phase D: did anything move? ----------------------------------------------
run "D  verify_published (full tier)" ov_D_verify.log \
    env VERIFY_FULL=1 Rscript verify_published.R

echo | tee -a $LOG
echo "=== queue finished $(date) ===" | tee -a $LOG
echo "built species now:" | tee -a $LOG
Rscript -e 'suppressMessages(library(data.table)); r<-fread("species_registry.csv")
r[, built := file.exists(sprintf("%s_cells_5km_16d.csv", prefix))]
cat(sprintf("  %d of %d built\n", sum(r$built), nrow(r)))
print(r[, .(n = sum(built), of = .N), by = guild])' 2>&1 | tee -a $LOG
df -h . | tail -1 | sed 's/^/disk: /' | tee -a $LOG
