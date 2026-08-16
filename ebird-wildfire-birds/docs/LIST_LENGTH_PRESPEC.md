# Pre-specification — list-length outcome test

**Written 2026-08-02, BEFORE any list-length data existed.** The point of this
file is that the pass/fail rule is fixed in advance. If the rule is changed
later, the change and its reason must be appended below, dated, not edited in
place.

## The problem this is meant to fix

Per-species occurrence effects are currently unusable. Three negative controls
(Mallard, American Coot, Black Phoebe) show larger apparent fire effects than the
focal species. The channel is observer behaviour: after a fire, birders visit
**43% fewer distinct locations** per cell (p = 0.0001).

Current state, from `point_controls_results.csv` / `RESULTS_WRITEUP.md` §2:

```
                  cell grain   point grain   p (point)
Mallard             +0.0265      +0.0148      0.032
American Coot       +0.0370      +0.0231      0.001
Black Phoebe        +0.0431      +0.0186      0.043
Wrentit             -0.0298      -0.0102      0.090

mean |effect|  CONTROLS : cell 0.0355 -> point 0.0188
mean |effect|  PREDICTED: cell 0.0290 -> point 0.0076
```

## The idea being tested

If fewer/shorter/less productive trips after a fire depress the recorded rate for
**every** species by a roughly common factor, then a measure that divides a
species out by total recording activity should cancel that factor. **List
length** — the number of species reported on a checklist — is the standard proxy
for "how hard did this person look and how good were conditions."

This is the same cancellation logic that makes the guild contrast work
(`guild_contrast.R`), applied at the species level instead of the guild level.

## What gets built

1. `list_length_lookup.R` — one pass over the raw EBD, counting distinct true
   species per complete checklist. Species-independent, computed once.
2. List length joined to the checklist table, aggregated to 5 km x 16-day cells.
3. Two new outcomes alongside `det_rate`:
   - **covariate form** — `det_rate` with mean list length as a control
   - **share form** — the species' detections normalised by total recording
     activity in the same cell-time
4. The three control species rerun on both.

## THE VERDICT RULE — fixed in advance

The new outcome is declared **CLEAN** only if **both** hold:

1. **All three controls** (Mallard, American Coot, Black Phoebe) return
   p >= 0.05, and
2. **Mean |effect| across the three controls drops below 0.0100**, i.e. a
   material fall from the current point-grain 0.0188 — not a rounding change.

If either condition fails, the outcome is declared **NOT CLEAN**, per-species
occurrence claims stay retired, and the negative result is reported in
`RESULTS_WRITEUP.md` as a fourth pre-specified check that failed.

Partial outcomes (e.g. two of three controls quiet) are **NOT CLEAN**. There is
no partial credit and no post-hoc subsetting to a favourable pair.

## What a CLEAN result would and would not establish

**Would:** that the shared component of the observer artifact is removable, and
that per-species occurrence work can resume on the corrected outcome.

**Would NOT:** vindicate the guild contrast. This cancels effects common to all
species. It does nothing about **differential detectability** — burning off dense
chaparral plausibly makes a skulking Wrentit easier to see while barely changing
a bird that already perches in the open. That is species-specific, would survive
this correction, and would mimic a guild pattern with no change in bird numbers.
It remains an open threat either way.

## Reporting commitment

The result is written into `RESULTS_WRITEUP.md` whichever way it comes out. A
failed repair is evidence that the standard fix does not work, which strengthens
the methods argument rather than weakening it.

---

## Amendments

*(none)*
