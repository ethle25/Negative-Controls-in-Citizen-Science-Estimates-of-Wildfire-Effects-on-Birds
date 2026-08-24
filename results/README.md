# results/

**Empty by design.** Everything here is generated and `.gitignore`d.

`main.R` writes `sessionInfo.txt` here at the end of each run. Include it when
reporting a result that does not reproduce.

The analysis scripts write their own outputs into `eBird/` alongside the code,
and figures land in `eBird/figures/`. That is a consequence of the scripts
resolving paths relative to a single directory — see README section 4.

| Output | Written by |
|---|---|
| `eBird/*_summary.txt` | Each analysis script — the formatted reports |
| `eBird/*_results.csv` | Each analysis script — the machine-readable results |
| `eBird/figures/fig[0-6]_*.png` and `.pdf` | `make_figures.R` |
| `results/sessionInfo.txt` | `main.R` |
