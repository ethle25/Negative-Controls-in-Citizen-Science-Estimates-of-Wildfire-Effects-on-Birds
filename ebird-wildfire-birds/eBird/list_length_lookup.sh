#!/bin/bash
# =============================================================================
# list_length_lookup.sh -- number of species reported on each complete checklist
#
# WHY
#   The project's effort features measure how MUCH effort was spent (hours,
#   distance, observers, number of checklists). None measures how PRODUCTIVE the
#   trip was. "List length" -- the count of species on a checklist -- is the
#   standard proxy for observer skill, conditions and general detectability
#   rolled into one number. See LIST_LENGTH_PRESPEC.md for what it is for and
#   the pass/fail rule fixed in advance.
#
#   It is SPECIES-INDEPENDENT: computed once, reusable by every species and
#   every analysis. That is why it is worth a pass over the raw EBD.
#
# WHY NOT FROM THE JOINED TABLE
#   `<prefix>_2015_2026_joined.csv` is one row per checklist for ONE species
#   (zero-filled). It cannot say how many OTHER species were on that checklist.
#   Only the raw EBD has that, and it is 52 GB across two files.
#
# METHOD
#   The EBD holds one row per species per checklist. So for a given checklist,
#   the number of rows with CATEGORY == "species" IS the number of species.
#   No de-duplication hash is needed -- which matters, because a hash over
#   (checklist x species) would be ~130M keys.
#
#   Columns (verified identical in both files 2026-08-02):
#     4  CATEGORY        6  COMMON NAME      35 SAMPLING EVENT IDENTIFIER
#     45 ALL SPECIES REPORTED               46 GROUP IDENTIFIER
#
#   `cut` first, then awk. awk splits all ~50 fields per line otherwise, and at
#   52 GB that difference is large.
#
# KNOWN UNDERCOUNT (deliberate, and consistent across checklists)
#   Only CATEGORY == "species" is counted. Rows recorded solely as a
#   sub-specific group ("Canada Goose (Aleutian)", CATEGORY == "issf") without
#   the parent species are missed. This slightly understates a few lists. It is
#   consistent across checklists, and list length is used as a RELATIVE measure
#   of effort, not an absolute species count, so a small uniform undercount does
#   not affect its role. `spuh`, `slash` and `hybrid` are excluded on purpose --
#   "gull sp." is not a species detected.
#
# SHARED CHECKLISTS -- THE JOIN KEY TRAP (found 2026-08-02, before first run)
#   `auk_unique()` collapses a shared checklist to its GROUP IDENTIFIER, so
#   `checklist_id` in the joined table is a group id ("G...") for 626,512 of
#   4,445,934 rows -- 14.1%. Keying this lookup on SAMPLING EVENT IDENTIFIER
#   alone would silently return NA for every one of them, and NA list length on
#   a systematically different subset (group birding trips are longer and more
#   productive) is exactly the kind of bias this file exists to remove.
#
#   So both keys are emitted. For a group, list length is the MAX across its
#   member events: members of a shared checklist report the same list in
#   principle, but each may edit their own copy, and max is the least
#   destructive reconciliation. "S..." and "G..." keys cannot collide.
#
# OUTPUT  list_length_lookup.tsv  --  checklist_id <TAB> n_species
#         keyed by BOTH sampling event id and group id
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")"

A="ebd_US-CA_201501_202312_smp_relJun-2026/ebd_US-CA_201501_202312_smp_relJun-2026.txt"
B="ebd_US-CA_202401_202606_smp_relJun-2026/ebd_US-CA_202401_202606_smp_relJun-2026.txt"
OUT="list_length_lookup.tsv"

for f in "$A" "$B"; do
  [ -r "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
done

# The two date ranges share no checklists, so the parts concatenate with no
# merge step -- which avoids a second hash the size of the first.
pass () {
  local src="$1" dst="$2"
  echo "[$(date +%H:%M:%S)] $src"
  # cut order is file order: 4=CATEGORY 35=EVENT 45=ALL_SPECIES 46=GROUP
  cut -f4,35,45,46 "$src" \
    | awk -F'\t' '
        $1 == "species" && $3 == "1" {
          n[$2]++
          if ($4 != "") grp[$2] = $4          # remember this event is in a group
        }
        END {
          for (k in n) {
            print k "\t" n[k]                 # sampling-event key
            if (k in grp) {                   # group key = max over members
              g = grp[k]
              if (!(g in gm) || n[k] > gm[g]) gm[g] = n[k]
            }
          }
          for (g in gm) print g "\t" gm[g]
        }' > "$dst"
  echo "[$(date +%H:%M:%S)]   -> $(wc -l < "$dst") keys"
}

echo "start $(date +%H:%M:%S)"
pass "$A" .ll_part_a.tsv
pass "$B" .ll_part_b.tsv

cat .ll_part_a.tsv .ll_part_b.tsv > "$OUT"
rm -f .ll_part_a.tsv .ll_part_b.tsv

echo "[$(date +%H:%M:%S)] TOTAL $(wc -l < "$OUT") checklists -> $OUT"
awk -F'\t' '{s+=$2; if($2>m) m=$2; n++}
            END { printf "  mean list length %.2f   max %d   n %d\n", s/n, m, n }' "$OUT"
echo "done $(date +%H:%M:%S)"
