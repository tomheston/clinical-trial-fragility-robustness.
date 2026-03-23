# Data Dictionary — corpus_empirical_trials_v2_2026-03-23.xlsx

**Version:** 1.0  
**Date:** 2026-02-01  
**Trials:** 100  
**Columns:** 20  

## Overview

Empirical corpus of binary-outcome randomized controlled trials used to validate the four-metric fragility-robustness framework (SFI, PFI, MFQ, RQ) in Heston (2026). Each row is one trial. All trials have clearly reported 2×2 cell counts from parallel-group designs with fixed follow-up periods.

**Inclusion criteria:** Parallel-group RCT; binary primary outcome; all four cell counts (a, b, c, d) explicitly reported; pharmaceutical or diagnostic intervention.

**Exclusion criteria:** Adaptive designs without stable final counts; multi-arm trials not reducible to a pairwise contrast; time-to-event outcomes without binary reduction.

---

## Columns

### Study Identification

| Column | Type | Description |
|--------|------|-------------|
| DOI | string | Digital Object Identifier — primary source reference |
| YEAR | integer | Publication year |

### Study Design

| Column | Type | Description |
|--------|------|-------------|
| ALLOCATION | float | Treatment-to-control allocation ratio (e.g., 1.0 = 1:1; 2.0 = 2:1). Range: 1.0–23.1 |
| SPECIALTY | string | Medical specialty (see codes below) |
| INTERVENTION_TYPE | string | Intervention category (see codes below) |
| PHASE | integer | Trial phase: 1, 2, 3, or 4 |
| EVENT | string | Name of the primary binary outcome as reported in the trial |

### Raw 2×2 Table

Arm A = treatment/intervention. Arm B = control/comparator.

| Column | Type | Description |
|--------|------|-------------|
| AEvent | integer | Arm A event count (cell a) |
| ANonEvent | integer | Arm A non-event count (cell b) |
| BEvent | integer | Arm B event count (cell c) |
| BNonEvent | integer | Arm B non-event count (cell d) |
| N | integer | Total sample size = a + b + c + d |
| pvalue | float | Observed p-value (Fisher's exact, two-sided, α = 0.05) |

### Fragility Metrics

| Column | Type | Description |
|--------|------|-------------|
| FI | integer | Classic Fragility Index — minimum toggle count to flip significance (Walsh et al., 2014) |
| FQ | float | Fragility Quotient = FI / N ∈ [0, 1] |
| SFI | integer | Standardized Fragility Index — toggle count using larger arm; allocation-invariant |
| MFQFI | float | Modified-arm Fragility Quotient = FI / n_mod ∈ [0, 1], where n_mod = arm subjected to toggling |
| PFIPearson | float | Percent Fragility Index (Pearson χ² fixed-margin path) ∈ [0, 1] |

### Robustness Metric

| Column | Type | Description |
|--------|------|-------------|
| RQ | float | Risk Quotient = \|ad − bc\| / (N²/4) ∈ [0, 1] — geometric distance from therapeutic neutrality |

### Pattern Flag

| Column | Type | Description |
|--------|------|-------------|
| PATTERN110 | float | Pattern (1,1,0) flag: 1.0 = trial is statistically significant + fragile (MFQFI ≤ 0.10) + weak robustness (RQ < 0.075). Null = pattern absent. 18 of 100 trials flagged. |

---

## Specialty Codes

| Code | Description |
|------|-------------|
| M-cardiology | Cardiology |
| M-oncology | Oncology |
| M-neurology | Neurology |
| M-infection_vaccines | Infectious disease and vaccines |
| M-pulmonary | Pulmonology |
| M-endocrinology | Endocrinology |
| M-autoimmunity | Autoimmune disease |
| M-hepatology | Hepatology |
| M-nephrology | Nephrology |
| M-OBGYN | Obstetrics and gynecology |
| M-Peds | Pediatrics |
| S-cardiovascular | Cardiovascular surgery |

## Intervention Type Codes

| Code | Description |
|------|-------------|
| P-Small Molecule Drug | Small molecule pharmaceutical |
| P-Biologic | Biologic agent |
| P-Vaccine | Vaccine |
| P-Hormone | Hormone therapy |
| P-peptide | Peptide therapy |
| P-Novel-Other | Novel pharmaceutical, other |
| D-Molecular | Molecular diagnostic |
| D-MRI_CT_XR_US | Imaging diagnostic |

---

## Notes

- **Event polarity:** EVENT records the primary outcome name as reported. Event direction (favorable vs. unfavorable) was not harmonized in this archive; consult the source trial via DOI.
- **Zero-cell policy:** Haldane–Anscombe 0.5 continuity correction applied only where required by the specific estimator. Not applied to RQ or raw toggle counts.
- **MFQFI denominator:** n_mod is the arm with fewer events; if tied, the smaller arm — consistent with the Walsh et al. (2014) toggle rule.
- **PATTERN110:** The clinically problematic pattern of significant + fragile + weak robustness per the p–fr–nb framework (Heston, 2026). Null displayed as #NULL! in Excel indicates the pattern is absent.

---

## Citation

Heston TF. *A Four-Metric Fragility-Robustness Framework for Binary Clinical Trials.* PhD Dissertation, Global Humanistic University, 2026.

## License

CC BY 4.0 — https://creativecommons.org/licenses/by/4.0/
