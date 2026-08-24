#!/bin/bash
# Predictive task, for every ACTIVE species in species_registry.csv, under both
# effort treatments.
#
# The species used to be a literal shell string -- "Wrentit:wrenti:wrentit ..." --
# so adding a bird to the registry did not reach this script. Now the list comes
# from the registry, and only species whose cell frame has actually been built
# are run (a registry row without data would otherwise fail four models deep).
#
#   ./run_models.sh                  every active species with a built frame
#   ./run_models.sh wrentit olsfly   just these prefixes
#   BURN_MEASURE=frac_recent ./run_models.sh   swap the burn block (see
#                                    treatment_registry.csv; `frame` source only)
set -e
cd "$(dirname "$0")"

if [ $# -ge 1 ]; then
  PREFIXES="$*"
else
  PREFIXES=$(Rscript -e '
    source("project_config.R")
    s <- cfg_species()
    s <- s[file.exists(sprintf("%s_cells_5km_16d.csv", prefix))]
    cat(s$prefix)' 2>/dev/null)
fi

if [ -z "$PREFIXES" ]; then
  echo "no species to run: no active row in species_registry.csv has a built" >&2
  echo "  *_cells_5km_16d.csv frame. Build one with ./run_species.sh <prefix>." >&2
  exit 1
fi

echo "species: $PREFIXES"
echo "burn block: ${BURN_MEASURE:-block (all 22 registry burn columns)}"
echo

for PFX in $PREFIXES; do
  META=$(Rscript -e "
    source('project_config.R')
    s <- cfg_species(active_only = FALSE)[prefix == '$PFX']
    if (!nrow(s)) quit(status = 3)
    cat(s\$common_name[1], '|', s\$species_code[1])
  " 2>/dev/null) || { echo "!! '$PFX' not in species_registry.csv -- skipping" >&2; continue; }

  NAME=$(echo "$META" | cut -d'|' -f1 | sed 's/^ *//; s/ *$//')
  CODE=$(echo "$META" | cut -d'|' -f2 | sed 's/^ *//; s/ *$//')

  for mode in covariate fixed1; do
    export EBIRD_SPECIES="$NAME" EBIRD_SPECIES_CODE="$CODE" EBIRD_PREFIX="$PFX" EFFORT_MODE="$mode"
    echo "=== $NAME / $mode ==="
    # keep going on failure so one bad species does not abort the whole run
    Rscript model_burn_effect.R > "run_model_${PFX}_${mode}.log" 2>&1 \
      && echo "    done" \
      || echo "    FAILED - see run_model_${PFX}_${mode}.log"
  done
done
