# Response to the editorial and methodological review

Status of each of the 13 checks. "Done" means the code was changed and re-run; the number in the
manuscript must be taken from the regenerated tables, not from an earlier version of this document.

| # | Check | Status | What changed |
|---|---|---|---|
| 1 | Synthetic null in Fig. 1f | **Done** | `03b` now writes every permutation draw to `results/gtex/sheaf_null_draws.tsv` and a summary to `sheaf_summary.tsv`; the figure plots those draws. Permutations raised from 100 to 200 (`T2D_SHEAF_PERM`), which changes **z from −26 to −29.7** (E = 0.953 against 1.911 ± 0.032). The earlier z was computed on 100 permutations while the legend said 200. |
| 2 | One version of the results | **In progress** | The analyses re-run so far write their own tables; the remaining step is a single clean execution of `run_replication.sh` from which text, legends and figures are regenerated. Until then, AUDIT.md and the manuscript may disagree and the manuscript is not final. |
| 3 | Pair identity in GSE22309 | **Checked, holds** | Sample titles carry consecutive numbering with basal before insulin (907/908, 911/912, 915/916 …), so `arange(n)//2` reproduces donor identity. Now documented in Methods. The check also confirms that hybridisation run differs within some pairs, which is why the same-run restriction is reported. |
| 4 | Source Data | **Partly done** | New exports: `sheaf_null_draws.tsv`, `sheaf_summary.tsv`, `canonical_summary.tsv`, `canonical_null_draws.tsv`, `blood_to_tissue.tsv`, `exercise_group_tests.tsv`, `power_and_confounding.tsv`. Still missing: per-panel source tables for Figs. 2–4 and the network metrics. |
| 5 | Methods must match the code | **Partly done** | Welch is now actually used (`equal_var=False` in `05` and `09`). Still to correct in the text: the Python spectral embedding uses `|bicor|^β`, i.e. unsigned; and the implemented energy is a normalised Frobenius difference between aligned subspaces, which is a Grassmann-type distance but should be named as implemented. |
| 6 | P values and permutation counts | **Done** | Defaults set to 2,000 for group comparisons and 3,000 for paired swaps, matching the legends. One-sided for coherence, two-sided for paired and per-gene tests, as stated. |
| 7 | Blood-to-tissue leakage | **Done** | Residualisation, feature selection, scaling and both PCAs are now fitted inside each training fold only. R² changes little: 0.296 → **0.279 ± 0.062** (muscle), 0.233 → **0.229 ± 0.037** (adipose), 0.228 → **0.217 ± 0.062** (pancreas), so the earlier estimate was not driven by leakage. |
| 8 | Resting association model | **Done** | Now a partial correlation: the M-value is residualised on the same covariates as expression. This **changes a number quoted in the manuscript**: genes at *P* < 0.01 go from 14 to **231, against 215 expected by chance**. No gene reaches FDR < 0.1 either way, but the claim becomes "as expected by chance" rather than "fewer than expected", which was an artefact of residualising only one side. Aggregate depth versus M-value: ρ = −0.37, *P* = 0.008. M-value units are taken as deposited and still need checking against the source publication. |
| 9 | GSE9105 timing metadata | **Resolved** | The series summary states biopsies basal and after 30 and 240 min of insulin infusion, and the sample titles agree (`Baseline_0min`, `Insulin_30min`, `Insulin_240min`). The "180-min" in `overall_design` is an inconsistency in the deposit. Our use of the 240-min labels is correct. Participants are 12 Mexican-American adults without a family history of diabetes, now recorded in Supplementary Table 1. |
| 10 | Exercise contrasts not traceable | **Done** | `11` now exports `exercise_group_tests.tsv` with Δcoherence, cosine and both P values. Re-run values: muscle at recovery Δ = −0.079, **P = 0.70**, cosine **0.84** (null 0.90); adipose at recovery Δ = −0.133, P = 0.76, cosine **0.92**. The direction of the difference favours T2D in both tissues, so the specificity argument stands, but it rests on a non-significant difference and must be phrased as such. |
| 11 | Robustness and power | **Partly done** | The null for the power simulation is now isotropic (`isotropic_null`), not a sign flip. The sign-flip null was indeed too lenient at the 95th percentile (e.g. n = 10, d = 5: 0.46 versus 0.30 isotropic), although both centre on zero. A TOST was added for magnitude: with a margin of 0.5 log2 units neither group is equivalent to the sensitive group (P = 0.26 and 0.60), so "magnitude preserved" must become "no reduction detected". Still to do: vary gene number, dimension and scaling. |
| 12 | Cohort inventory | **Done** | `13_cohort_inventory.py` builds Supplementary Table 1 and counts the three quantities separately: **24 independent studies, 29 accessions, 1,766 participants summed over accessions and 1,689 after removing donors shared between accessions of the same study** (GSE50398 overlaps GSE50244; GSE182117/182120 belong to the superseries GSE182121). Seventeen cohorts have paired sampling around a perturbation, totalling 646 participants. Platform and sample number are read from the series matrices where available (23 of 29). The manuscript's "23 cohorts" and "25 cohorts" are both replaced by these figures. |
| 13 | Reproducible archive | **Open** | Needs a tagged release with a lockfile; `/home/claude` paths removed from `03b`. |

## Items where the review changes what we can claim

- "Magnitude is preserved" is not supported by the equivalence test and becomes "no reduction in magnitude was detected".
- The specificity of the loss to insulin rests on a comparison between separate cohorts and protocols and on a non-significant group difference within the exercise cohorts; it is evidence of stimulus dependence, not proof of specificity.
- The two prediabetes arms are two arms of one trial, not independent cohorts.
- "Clock output" includes both output factors and regulatory components; "clock-related" is the accurate term for the genes examined.

## Characterisation of the measure (added after the review)

The review asked what the statistic does in cases where no coordinated response is expected. Six
perturbations were analysed with the same pipeline:

| Perturbation | Time scale | Tissue | *n* | Coherence |
|---|---|---|---|---|
| Insulin, hyperinsulinaemic clamp, 4 h (sensitive) | hours | muscle | 20 | 0.77 |
| Acute exercise, before training (GSE224310) | hours | muscle | 25 | **0.80** (*P* = 0.0005) |
| Acute exercise, after training (GSE224310) | hours | muscle | 23 | 0.70 (*P* = 0.0005) |
| Acute exercise, 3 h recovery (GSE202295) | hours | muscle | 17–20 | 0.45–0.53 |
| Acute exercise (GSE224310) | hours | adipose | 12–14 | 0.31–0.38 |
| Training, weeks (GSE224310) | weeks | muscle | 24 | 0.31 |
| Training, weeks (GSE224310) | weeks | adipose | 13 | −0.09 (*P* = 0.72) |
| Diurnal change, aligned (GSE106800) | 12 h | muscle | 12 | 0.10 |
| Diurnal change, misaligned (GSE106800) | 12 h | muscle | 11 | −0.07 |
| Cold acclimation, 10 days (GSE67297) | days | adipose | 7 | −0.01 (*P* = 0.54) |
| Palmitate ex vivo, 4 days (GSE159984) | days | islet | 5 | 0.42 |

Three consequences. First, coherence near 0.8 is what a healthy acute response looks like in muscle,
which places the 0.77 observed under insulin in context and shows it is not an artefact of the 2007
platform or its batch structure, since GSE224310 is an independent cohort on another platform.
Second, the contrast between acute and chronic is demonstrated **within the same participants, tissue
and platform** (0.80 versus 0.31 in muscle; 0.31–0.38 versus −0.09 in adipose), so slow perturbations
produce individual rather than shared trajectories. Third, adipose tissue shows lower coherence than
muscle even in healthy people and even for an acute stimulus, which qualifies the two-failure-modes
statement in Fig. 6.

Circadian misalignment does **not** provide the mechanistic link we had hoped for: the contrast that
design measures is evening versus morning twelve hours apart, which is not a coordinated response even
in the aligned condition, so the comparison is uninformative rather than negative. Any argument that
misalignment reduces coordination has been removed from the plan.
