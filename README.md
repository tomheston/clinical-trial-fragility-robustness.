# clinical-trial-fragility-robustness

**A four-metric fragility-robustness framework for binary clinical trial evidence quality assessment.**

**Author:** Thomas F. Heston, MD  
**ORCID:** [0000-0002-5655-2512](https://orcid.org/0000-0002-5655-2512)  
**Affiliation:** University of Washington; Washington State University  
**Dissertation:** PhD Business Administration, Global Humanistic University, 2026 [in progress]  
**License:** CC BY 4.0  

---

## What This Is

This repository contains the empirical data and simulation code supporting a PhD dissertation that proposes four metrics for evaluating the quality of statistical evidence in binary clinical trials — going beyond the p-value alone.

The framework addresses a fundamental problem: a p-value of 0.049 and a p-value of 0.0001 are both reported as "statistically significant," but they represent vastly different levels of evidence quality. These four metrics quantify the missing dimensions.

---

## The Four Metrics

| Metric | Full Name | Measures | Range |
|--------|-----------|----------|-------|
| **SFI** | Standardized Fragility Index | Minimum outcome toggles (larger arm) to flip significance | 0–∞ (integer) |
| **PFI** | Percent Fragility Index | Fixed-margin continuous shift to flip significance | [0, 1] |
| **MFQ** | Modified-arm Fragility Quotient | Proportion of toggled arm required to flip significance | [0, 1] |
| **RQ** | Risk Quotient | Geometric distance from therapeutic neutrality | [0, 1] |

SFI, PFI, and MFQ measure **fragility** — how stable is the significance classification?  
RQ measures **robustness** — how far is the result from no effect?

Together they form the **p–fr–nb triplet**: significance, fragility, and robustness — complete statistical evidence.

---

## Repository Structure

```
clinical-trial-fragility-robustness/
├── data/
│   ├── corpus_empirical_trials_v2_2026-03-23.xlsx   # 100 binary-outcome clinical trials
│   └── README_empirical.md                           # Data dictionary
└── code/
    ├── sim_SFW.R                                     # Main Monte Carlo simulation (canonical)
    ├── rq_cutoff_sim_2026-02-01_final.R              # RQ empirical cutoff simulation
    ├── mfq_vs_allocation_sim_2026-02-01.R            # MFQ allocation independence test
    ├── generate_figures_ch5_2026-02-01.R             # Chapter 5 figure generation
    └── README_code.md                                 # Code documentation
```

---

## Empirical Corpus

100 binary-outcome clinical trials (Phase I–IV) with published 2×2 contingency tables. Includes raw cell counts, computed fragility and robustness metrics, and the Pattern (1,1,0) flag — the clinically problematic pattern of significant + fragile + weak robustness.

See `data/README_empirical.md` for full data dictionary.

---

## Simulations

Monte Carlo simulations across 720,000 synthetic trials characterize metric operating characteristics under null and alternative conditions. Key finding: 18% of pharmaceutical trials in the empirical corpus exhibit Pattern (1,1,0) versus 1.34% expected under the null — a 13.4× enrichment.

See `code/README_code.md` for execution instructions.

---

## Key Definitions

**Fragility:** The minimum perturbation to observed data required to reverse the statistical significance classification. High fragility (low metric value) means a small change flips the conclusion.

**Robustness:** Geometric distance from therapeutic neutrality (RR = 1, no effect). High RQ means the result is far from "no difference."

**Pattern (1,1,0):** A trial is statistically significant (p ≤ 0.05), fragile (MFQ ≤ 0.10), and weakly robust (RQ < 0.075). Evidence consistent with a trivial or unreliable effect.

**MFQ formula:** MFQ = FI / n_mod, where n_mod = sample size of the arm subjected to toggling.

**RQ formula:** RQ = |ad − bc| / (N²/4), where a, b, c, d are the 2×2 table cells and N is total sample size.

---

## Citation

Heston TF. *A Four-Metric Fragility-Robustness Framework for Binary Clinical Trials.* PhD Dissertation, Global Humanistic University, 2026 [in progress].

Heston TF. Fragility Metrics Toolkit. Zenodo. 2025. https://doi.org/10.5281/zenodo.17254763

---

## License

CC BY 4.0 — https://creativecommons.org/licenses/by/4.0/
