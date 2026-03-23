# Code — clinical-trial-fragility-robustness

**Language:** R  
**R version tested:** 4.3+  
**Dependencies:** base R (most scripts); `ggplot2`, `dplyr`, `tidyr`, `scales` (figure generation only — auto-installed if missing)

---

## Scripts

### 1. `sim_SFW.R` — Main simulation pipeline (canonical)

Runs the full Monte Carlo simulation for Chapter 5 (§5.5–§5.6). Computes SFI, PFI, MFQ, and RQ across a factorial design grid and characterizes Pattern (1,1,0) — the SFW pattern.

**Inputs:** None (generates synthetic 2×2 tables)

**Key parameters (configurable via environment variables):**

| Variable | Default | Description |
|----------|---------|-------------|
| `FRAGILITY_SEED` | 20251203 | Random seed |
| `FRAGILITY_ALPHA` | 0.05 | Significance threshold |
| `FRAGILITY_REPS` | 2000 | Replicates per scenario cell |
| `FRAGILITY_OUT` | `out/` | Output directory |

**Design grid:**
- Sample sizes (N): 60, 100, 200, 400, 800
- Allocation ratios: 1:1, 2:1, 3:2
- Control event rates: 0.05, 0.10, 0.20, 0.40
- True RR: 0.60–1.40 (13 values)

**Outputs (written to `FRAGILITY_OUT`):**

| File | Contents |
|------|----------|
| `replicates_raw.csv` | All individual simulation replicates |
| `operating_characteristics.csv` | Aggregated metric summaries by scenario |
| `table_5_9_property_checks.csv` | Formal property verification (Table 5.9) |
| `correlations_by_scenario.csv` | Spearman correlations among metrics by scenario |
| `thresholds_used.csv` | RQ and MFQ cutoffs derived from null distribution |
| `fpr_overall_pattern110.csv` | False positive rate for Pattern (1,1,0) overall |
| `fpr_by_design_pattern110.csv` | False positive rate by design cell |
| `perf_mfq_*.csv` | MFQ sensitivity/specificity at varying thresholds |
| `perf_rq_*.csv` | RQ sensitivity/specificity at varying thresholds |
| `roc_curve_mfq_sigReq.csv` | ROC data for Figure 5.1 |
| `sfw_unconditional_by_rr.csv` | Pattern (1,1,0) rate by RR_true |
| `sfw_conditional_by_rr.csv` | Pattern (1,1,0) rate conditional on significance, by RR_true |
| `sfw_conditional_by_scenario.csv` | Pattern (1,1,0) rate conditional on significance, by full scenario |

**Run:**
```r
Rscript sim_SFW.R
```

---

### 2. `rq_cutoff_sim_2026-02-01_final.R` — RQ empirical cutoff simulation

Establishes empirical percentile cutoffs for RQ under the null hypothesis (RR = 1.0) using 50,000 simulated 2×2 tables drawn from realistic trial parameter distributions.

**Inputs:** None

**Outputs:**
- `rq_simulation_results.csv` — raw simulation results
- `rq_distribution.png` — RQ distribution histogram with 33rd/67th percentile lines

**Run:**
```r
Rscript rq_cutoff_sim_2026-02-01_final.R
```

---

### 3. `mfq_vs_allocation_sim_2026-02-01.R` — MFQ allocation independence test

Demonstrates that MFQ eliminates the allocation bias present in FQ. Simulates trials across five allocation ratios and tests whether ρ = MFQ/(2×FQ) scales linearly with allocation ratio — confirming FQ is allocation-biased while MFQ is not.

**Inputs:** None

**Outputs:**
- `dissertation_mfq_allocation_test.csv` — raw results
- `dissertation_rho_by_allocation.png` — boxplot of ρ by allocation ratio

**Run:**
```r
Rscript mfq_vs_allocation_sim_2026-02-01.R
```

---

### 4. `generate_figures_ch5_2026-02-01.R` — Chapter 5 figure generation

Produces Figures 5.1, 5.2, and 5.3 from simulation output CSVs.

**Inputs (must exist in same directory or `FRAGILITY_OUT`):**
- `roc_curve_mfq_sigReq.csv` (from `sim_SFW.R`)
- `correlations_by_scenario.csv` (from `sim_SFW.R`)

**Outputs (written to `out_figures/`):**

| File | Figure |
|------|--------|
| `Figure_5_1_MFQ_ROC_Curve_Zoomed.png/.pdf` | MFQ ROC curve |
| `Figure_5_2_SFI_RQ_Correlation_Heatmap.png/.pdf` | SFI–RQ correlation heatmap |
| `Figure_5_3_Pattern_110_Component_Analysis.png/.pdf` | Pattern (1,1,0) flow diagram |

**Run:**
```r
Rscript generate_figures_ch5_2026-02-01.R
```

---

## Recommended Execution Order

```
1. sim_SFW.R                          # generates all simulation data
2. rq_cutoff_sim_2026-02-01_final.R   # generates RQ cutoffs independently
3. mfq_vs_allocation_sim_2026-02-01.R # generates allocation test independently
4. generate_figures_ch5_2026-02-01.R  # requires output from step 1
```

---

## Metric Implementations

All metric computations are self-contained within `sim_SFW.R`:

| Function | Metric |
|----------|--------|
| `compute_fi_fq_mfq()` | FI, FQ, MFQ (Walsh et al. 2014 toggle rule) |
| `compute_sfi()` | SFI (larger arm, BFU = 1/n_large) |
| `compute_pfi_along_path()` | PFI (Pearson χ² fixed-margin path) |
| `RQ_metric()` | RQ = \|ad − bc\| / (N²/4) |

---

## Citation

Heston TF. *A Four-Metric Fragility-Robustness Framework for Binary Clinical Trials.* PhD Dissertation, Global Humanistic University, 2026 [in progress].

## License

CC BY 4.0 — https://creativecommons.org/licenses/by/4.0/
