# Replication table (all cohorts analysed, October 2026)

Discovery: GSE76895 (islet), GSE18732 (muscle), GSE15653 (liver), GSE27951 (adipose).
29 datasets received, 23 analysed. Statistics as defined in `python/` modules; n per group in brackets.

## R1. A systemic transcriptional state shared across organs and specific to the individual

| Evidence | Cohort | n | Statistic | Result |
|---|---|---|---|---|
| Cross-tissue network coherence (sheaf) | discovery (3 tissues, separate cohorts) | 83/118/33 | z vs gene-correspondence null | −10.6 / −7.5 / −7.0 (healthy / IGT / T2D) |
| Same, **paired donors** | GTEx v11 (muscle, SAT, pancreas) | 253 donors | z | **−26** (RIN, ischaemic time, batch, age, sex, Hardy regressed) |
| Individual-level coherence | GTEx (muscle–SAT, muscle–pancreas, SAT–pancreas) | 253 | canonical ρ vs donor-shuffle | **0.93 / 0.91 / 0.89 vs 0.23** (p = 0.002) |
| Same, **in vivo**, two depots | GSE20950 (SAT–omental), bariatric patients | 19 | ρ vs shuffle | **0.93 vs 0.59** (p = 0.0005) |
| Blood reflects tissue state (same donor) | GTEx (blood → muscle/SAT/pancreas) | 224 | 5-fold R² | 0.30 / 0.23 / 0.23 |

**Status: replicated (post mortem and in vivo). Caveat: unmeasured donor-level factors in GTEx; disease label unavailable.**

## R2. No second attractor; continuous drift

| Cohort | Tissue | n | Stable basins | ΔBIC(2 vs 1) | Cramér's V (component × stage) |
|---|---|---|---|---|---|
| GSE76895 | islet | 83 | 1 | +37 | 0.12 |
| GSE18732 | muscle | 118 | 1 | +95 | 0.04 |
| GSE27951 | adipose | 33 | 2* | +62 | 0.25 (second component = 1 outlier) |
| GSE164416 (batch-balanced) | islet, living donors | 73 | 1 | +62 | 0.32 |
| GSE50244 | islet, RNA-seq | 77 | 1 | −27 | 0.17 |
| GSE50398 (same donors, array) | islet | 77 | 1 | −28 | 0.28 |
| GSE25462 | muscle | 50 | 1 | −17 | — |
| METSIM GSE135134 | adipose, by BMI | 434 | 1 | −31 | 0.09 |

**Status: replicated in 8 cohorts, 3 tissues, 2 platforms. The GMM mixtures favoured by BIC never align with stage.** *Adipose discovery not supported on inspection.

## R3. Fastest transcriptional change inside the prediabetic range

| Cohort | Control variable | Fisher–Rao peak | Stage at peak |
|---|---|---|---|
| GSE18732 muscle | HbA1c | **6.27 %** | T2D-side of prediabetes |
| GSE76895 islet | fasting glucose | 5.81 mmol/L | healthy side |
| GSE27951 adipose | fasting glucose | 5.51 mmol/L | intermediate |
| GSE50244 islet (RNA-seq) | HbA1c | **5.86 %** | intermediate |
| GSE50398 islet (array, same donors) | HbA1c | **5.78 %** | intermediate |

**Status: replicated across cohorts and platforms (diagnostic thresholds: HbA1c 6.5 %, glucose 7.0 mmol/L).**

## R4. Loss of the coordinated response to insulin (within-person perturbation, batch-immune)

Direction coherence = mean cosine of individual responses to the group mean response; cos = cosine of group direction to healthy direction.

| Organ, state | Coherence | cos to healthy | Magnitude | Cohort (n) |
|---|---|---|---|---|
| Muscle, healthy, 4 h clamp | **0.80** | 1 | ref | GSE22309 (20) |
| Muscle, healthy, 4 h | **0.77** | — | — | GSE9105 (12); at 30 min: 0.29 |
| Muscle, healthy, 2 h | 0.49 | — | — | GSE7146 (6) |
| Muscle, insulin-resistant | 0.45 | 0.58 | = | GSE22309 (20) |
| Muscle, prediabetes | 0.43 / 0.27 | — | — | GSE157988 placebo / NMN arms (12/11) |
| Muscle, T2D | 0.55 | **0.25** | ↓ (net −50 %, p = 0.03) | GSE22309 (15) — **single cohort** |
| Adipose, non-obese | 0.53 | 1 | ref | Rydén CAGE (23) |
| Adipose, sensitive | 0.89 | 1 | ref | GSE26637 (5) |
| Adipose, obese | 0.35 (p = 0.07) | 0.88 | **½** | Rydén (23) |
| Adipose, resistant | 0.77 | 0.95 | **½** | GSE26637 (5) |
| Adipose, 2 y after bariatric surgery | 0.33 | 0.86 | **restored** (7.9 vs 7.8) | Rydén (23, same women) |
| Muscle, prediabetes after 10 wk NMN | 0.18 (vs 0.27, p = 0.57) | — | — | GSE157988 |

**Status: replicated for healthy coordination (2 muscle cohorts), for loss of coordination in the insulin-resistant / prediabetic state (3 muscle cohorts), and for the half-magnitude same-direction adipose response (2 cohorts). The orthogonal direction of diabetic muscle rests on one cohort.**

## Withdrawn after replication

| Claim | Why |
|---|---|
| State-space contraction healthy → T2D (discovery islet, muscle, adipose) | Reversed in two islet cohorts (GSE164416, GSE50244/50398), flat in muscle (GSE25462), expands with BMI in METSIM; in longitudinal designs (GSE59034, Rydén, GSE157988) group dispersion is confounded with time-point batch in both directions. |
| IGT as a critical state (DNB-type signature) | I_c monotone or flat in discovery and replication islets; weak non-significant IGT maximum only in GSE164416. |
| Depth score readable in blood | No separation of prediabetes/T2D from controls in GSE156993, GSE21321 (n = 6–9/group); only poorly controlled T2D (HbA1c 11 %) separates. |
| Diurnal response as a third perturbation | GSE104674: n = 6/6, opposite direction, uninformative. |

## Still open and what decides it

| Question | Dataset | Access |
|---|---|---|
| Orthogonal insulin response in diabetic muscle, second cohort | mixed-meal muscle study (Physiol Genomics 2023) or any clamp with T2D | GEO, accession in paper |
| Systemic state in vivo by glycaemic stage, paired tissues | FUSION (331, muscle + SAT, NGT/IFG/IGT/T2D) | dbGaP phs001048 |
| Blood readout with disease label | GSE26168 (whole blood, IFG/T2D); KORA, FHS, MESA for prospective | GEO / controlled |
| Disease label for GTEx paired tissues | GTEx MHT2D | dbGaP phs000424 |
