#!/bin/bash
# Run the pipeline for one species, or for every ACTIVE species in
# species_registry.csv.
#
#   ./run_species.sh                     every active species, build + analyse
#   ./run_species.sh wrentit             just that prefix
#   ./run_species.sh wrentit brncre      several prefixes
#   STAGES=analyse ./run_species.sh      skip join+regrid (they are the slow part)
#   STAGES=build   ./run_species.sh      rebuild data only, no analysis
#   DOSE_MEASURE=two_dose_sev STAGES=analyse ./run_species.sh    a different dose
#
# The old positional form still works and is unchanged in behaviour:
#   ./run_species.sh "Brown Creeper" brncre brncre
#
# WHAT CHANGED (2026-07-31)
#   This used to run join -> regrid -> dose, i.e. 3 of the 13 analysis scripts, so
#   "run everything for a new bird" was still a manual list. It now runs every
#   per-species analysis, then the three cross-species reports once at the end.
#
# COST. The two build stages are ~12 min per species. Of the analyses,
# fire_cluster_se and dose_burn_effect dominate because of their bootstrap
# subsets, so both subset flags default to OFF here -- set DOSE_SUBSETS=1 /
# FC_SUBSETS=1 deliberately when you want them (~45 min each, per species).
#
# Species metadata (common name, code, guild, predicted direction) is read from
# species_registry.csv, and the burn measure from treatment_registry.csv, so
# adding either means adding a row -- not editing this file.
set -e
cd "$(dirname "$0")"

DOSE_SUBSETS="${DOSE_SUBSETS:-0}"
FC_SUBSETS="${FC_SUBSETS:-0}"
STAGES="${STAGES:-all}"
export DOSE_SUBSETS FC_SUBSETS

case "$STAGES" in all|build|analyse) ;; *)
  echo "STAGES must be all, build or analyse (got '$STAGES')" >&2; exit 2 ;;
esac

# --- legacy positional form: "Common Name" code prefix ------------------------
if [ $# -eq 3 ] && [[ "$1" == *" "* ]]; then
  export EBIRD_SPECIES="$1" EBIRD_SPECIES_CODE="$2" EBIRD_PREFIX="$3"
  echo "=== $1 ($2) [legacy invocation] ==="
  Rscript join_species_all_years.R > "run_${3}_join.log"   2>&1; echo "  join   done"
  Rscript regrid_cells.R   > "run_${3}_regrid.log" 2>&1; echo "  regrid done"
  Rscript dose_burn_effect.R       > "run_${3}_dose.log"   2>&1; echo "  dose   done"
  exit 0
fi

# --- registry-driven ----------------------------------------------------------
if [ $# -ge 1 ]; then
  PREFIXES="$*"
else
  PREFIXES=$(Rscript -e 'source("project_config.R"); cat(cfg_species()$prefix)' 2>/dev/null)
fi

if [ -z "$PREFIXES" ]; then
  echo "no species to run (check for active=1 rows in species_registry.csv)" >&2
  exit 1
fi

# Per-species analyses, in dependency order. cell_fire_lookup must precede
# anything that clusters by fire; it writes one shared lookup, which is safe
# because attach_fire() joins it on (gx,gy) rather than cell_id.
BUILD_STAGES="join_species_all_years.R regrid_cells.R"
ANALYSE_STAGES="cell_fire_lookup.R dose_burn_effect.R did_burn_effect.R \
fe_decompose.R fire_cluster_se.R wild_cluster_boot.R cs_estimator.R effort_placebo.R"

echo "species to run : $PREFIXES"
echo "stages         : $STAGES"
echo "burn measure   : ${DOSE_MEASURE:-frac_max}"
echo "bootstrap subsets: DOSE_SUBSETS=$DOSE_SUBSETS FC_SUBSETS=$FC_SUBSETS"
echo

FAILED=""
for P in $PREFIXES; do
  META=$(Rscript -e "
    source('project_config.R')
    s <- cfg_species(active_only = FALSE)[prefix == '$P']
    if (!nrow(s)) quit(status = 3)
    cat(s\$common_name[1], '|', s\$species_code[1], '|', s\$guild[1], '|', s\$predicted[1])
  " 2>/dev/null) || { echo "!! '$P' not in species_registry.csv -- skipping" >&2; continue; }

  NAME=$(echo "$META" | cut -d'|' -f1 | sed 's/^ *//; s/ *$//')
  CODE=$(echo "$META" | cut -d'|' -f2 | sed 's/^ *//; s/ *$//')
  GUILD=$(echo "$META" | cut -d'|' -f3 | sed 's/^ *//; s/ *$//')
  PRED=$(echo "$META" | cut -d'|' -f4 | sed 's/^ *//; s/ *$//')

  echo "=== $NAME ($CODE) - guild: $GUILD, predicted: $PRED ==="
  export EBIRD_SPECIES="$NAME" EBIRD_SPECIES_CODE="$CODE" EBIRD_PREFIX="$P"

  # NB plain if-blocks, not `[ a ] || [ b ] && list=...`. That parses as
  # (a || b) && c, and when the condition is false the whole list returns
  # non-zero -- which under `set -e` exits the script instead of skipping a stage.
  STAGE_LIST=""
  if [ "$STAGES" = "all" ] || [ "$STAGES" = "build" ]; then
    STAGE_LIST="$BUILD_STAGES"
  fi
  if [ "$STAGES" = "all" ] || [ "$STAGES" = "analyse" ]; then
    STAGE_LIST="$STAGE_LIST $ANALYSE_STAGES"
  fi

  for S in $STAGE_LIST; do
    TAG=$(basename "$S" .R)
    # keep going on failure so one bad species does not abort a 12-species run,
    # but stop THIS species if a build stage fails -- everything after needs it
    if Rscript "$S" > "run_${P}_${TAG}.log" 2>&1; then
      printf "  %-26s done\n" "$TAG"
    else
      printf "  %-26s FAILED - see run_%s_%s.log\n" "$TAG" "$P" "$TAG"
      FAILED="$FAILED ${P}/${TAG}"
      case " $BUILD_STAGES " in
        *" $S "*)
          echo "  build stage failed -- skipping the rest of $P"
          break ;;
      esac
    fi
  done
  echo
done

# --- cross-species reports: these loop over the registry internally, so they run
# --- ONCE at the end rather than per species
if [ "$STAGES" != "build" ]; then
  echo "=== cross-species reports ==="
  for S in event_study.R severity_effect.R migration_phenology.R; do
    TAG=$(basename "$S" .R)
    unset EBIRD_SPECIES EBIRD_SPECIES_CODE EBIRD_PREFIX
    if Rscript "$S" > "run_${TAG}.log" 2>&1; then
      printf "  %-26s done\n" "$TAG"
    else
      printf "  %-26s FAILED - see run_%s.log\n" "$TAG" "$TAG"
      FAILED="$FAILED ${TAG}"
    fi
  done
  echo
fi

if [ -n "$FAILED" ]; then
  echo "FAILED stages:$FAILED"
  echo "Everything else completed. Inspect the logs above before quoting results."
else
  echo "all stages completed."
fi
echo "Predicted directions are recorded in species_registry.csv --"
echo "compare them against the results WITHOUT editing that file."
