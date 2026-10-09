# Audit of results and methodology (October 2026)

Each claim is re-tested with the stricter version of its statistic and judged on three
levels: **survives**, **weakened** (direction holds, significance does not), **fails**.
Numbers are from re-analyses run for this audit.

## A. Claims that survive

### A1. A systemic transcriptional state shared across organs and specific to the individual
- Sheaf coherence in paired GTEx tissues z = −26 after RIN, ischaemic time, batch, age, sex, Hardy. The null (gene-correspondence shuffle) is the right one for "same architecture"; the effect size is not an artefact of gene properties alone because the shuffle preserves every per-tissue property.
- Individual-level coherence: canonical ρ 0.89–0.93 vs 0.23 under donor shuffle (p = 0.002); robust to the same covariates. **In vivo replicate** GSE20950: 0.93 vs 0.59 (p = 0.0005).
- Residual risk: an unmeasured donor-level factor (terminal illness, medication) in GTEx; the genes carrying the axis (FKBP5, MT1X/2A, SERPINA3, XBP1, IL1RL1) are consistent with agonal stress. The GSE20950 replicate (living surgical patients) argues the phenomenon is not only agonal, but its n is 19.
- **Verdict: survives.** The decay of coherence with stage (z −10.6 → −7.5 → −7.0; E differences p = 0.15–0.53) is a trend, not a result; write it as such.

### A2. One basin, continuous drift
- 8 cohorts, 3 tissues, 2 platforms, no stable second basin; GMM components unrelated to stage (Cramér's V 0.04–0.32; adipose "second component" = 1 outlier).
- Limitation: tested on the first 3 PCs of 800 HV genes; a bimodality confined to a minor direction would be missed. Bridge flux result (J ≈ gradient) rests on the ergodicity assumption.
- **Verdict: survives**, as a transcriptomic contribution to the continuum-vs-subtypes debate, with the dimensionality caveat stated.

### A3. The basal muscle transcriptome does not encode insulin sensitivity
- GSE182120 (49 clamped biopsies, two labs regressed, age/BMI regressed): 0 genes at FDR < 0.1 vs M-value; 14 at p < 0.01 where ~215 are expected by chance; clock-output genes flat (all p > 0.2). Positive controls present (PPARGC1A ρ = 0.33, PDK4 ↑ T2D).
- **Verdict: survives, and it is the cleanest result in the paper.** It reframes everything else: the information is in the response.

### A4. Healthy insulin response is coordinated across individuals and assembled over hours
- Leave-one-out coherence (removes the self-inclusion bias of the naive statistic): IS 0.77 (n = 20, GSE22309); 0.77 at 4 h vs 0.29 at 30 min (GSE9105, n = 12); 0.49 at 2 h (GSE7146, n = 6); 0.51 at 1 h meal (GSE231509).
- Caveat: repeat-biopsy artefact. A second biopsy near the first induces FOS/JUN/EGR1 regardless of insulin; the "immediate-early" part of the healthy programme is partly this. It does not affect between-group comparisons (all groups had repeat biopsies) but it does mean the gene list of the "healthy programme" overstates insulin specificity. Tuvia et al. 2021 had a saline arm in adipose; no public muscle dataset does.
- **Verdict: survives** for the time-course and for coordination; the gene-level "programme" needs the repeat-biopsy caveat.

## B. Claims that are weakened

### B1. Loss of coordination in insulin resistance and T2D
- Full cohort, LOO coherence: IS 0.77 vs IR 0.35 (p = 0.013) vs T2D 0.45 (p = 0.087). Direction of T2D response: cos 0.20 to healthy vs 0.69 under permutation (p = 0.034).
- **Batch problem found in audit:** sequencing "Run" is confounded with group, and only 64 % of basal/insulin pairs were hybridised on the same run. Restricting to same-run pairs (IS 11, IR 13, T2D 11): IS 0.76 vs IR 0.45 (p = 0.11) vs T2D 0.46 (p = 0.15); cos(IS, T2D) = 0.14. Effect sizes hold; significance does not at n ≈ 11.
- Replicates in prediabetes: GSE157988 (coherence 0.27 / 0.43, n = 11/12, RNA-seq, one lab) — these need the LOO recomputation but are far below 0.77 either way. Adipose: obese/IR respond at half magnitude, same direction (Rydén n = 23, GSE26637 n = 5).
- **Verdict: weakened.** The insulin-resistant loss of coordination is supported by three cohorts in direction and by one in significance; the diabetic orthogonal direction rests on 11–15 subjects with a batch caveat. Must be written as "consistent across three cohorts, significant in one" and the T2D direction as suggestive.

### B2. Loss of the insulin → clock-output coupling
- Interaction tests (IS vs IR, permutation, FDR over 9 clock genes): DBP p = 0.003, PER2 0.005, NR1D2 0.012 (FDR 0.016–0.024); BHLHE40 FDR 0.08. IS vs T2D: DBP, NR1D2, HLF, TEF, BHLHE40 all p 0.02–0.05, FDR 0.06.
- Same-run pairs: DBP response IS −0.13 / IR −0.03 / T2D +0.02; PER2 +0.31 / +0.17 / +0.31; NR1D2 −0.19 / −0.49 / −0.43: DBP loss holds; PER2 does not in T2D; NR1D2 shows a *stronger* repression in IR/T2D (a change of pattern, not a simple loss).
- In vitro convergence (GSE182117): chronic insulin represses DBP in NGT (t −3.3) more than T2D (−1.7); DBP/TEF/HLF/PER2 amplitude ~30 % lower in T2D myotubes (n = 7/5, p ≈ 0.15). Basal clock output unrelated to M-value (A3).
- Prior art: insulin resets PER2 in human adipose (Tuvia 2021); clock disrupted in T2D myotubes (Gabriel 2021). New here: the *coupling* is lost in vivo in muscle of insulin-resistant people, with DBP as the most robust marker.
- **Verdict: weakened but real for IR (3 genes survive FDR in the full cohort; DBP survives the same-run restriction); borderline for T2D.** It is the paper's most novel biological statement and its most fragile; it needs the Zierath cohort or a new clamp.

### B3. Interventions restore magnitude, not coordination
- Bariatric surgery (Rydén, n = 23): magnitude 5.0 → 7.9 (non-obese 7.8); coherence 0.35 → 0.33 (p = 0.55). NMN (n = 11): 0.27 → 0.18 (p = 0.57).
- "Not restored" is absence of evidence at n = 11–23; power to detect a return to 0.5–0.7 is modest.
- **Verdict: weakened.** Report as "no detectable restoration" with the power statement.

## C. Claims that fail

### C1. Fastest transcriptional change inside the prediabetic range (Fisher–Rao peak)
- Bootstrap CI of the peak: GSE50244 [5.45, 5.99], GSE50398 [5.36, 6.17]. **Under θ-shuffle the peak also lands at HbA1c 5.69 (IQR 5.50–5.86)**, because the windowed estimator peaks where sampling density is highest (median HbA1c 5.60). The observed peaks (5.86, 5.78) sit inside the null IQR.
- **Verdict: fails.** The "prediabetic window" is a density artefact of the estimator. Withdrawn. (A density-corrected estimator could be built, but there is no evidence to rescue at present.)

### C2–C4. Already withdrawn: state-space contraction (reversed in replication; batch-confounded longitudinally); IGT as a critical state; blood readout; immune decoupling (no FDR).

## D. Methodological corrections adopted

1. Coherence is reported leave-one-out, never naive.
2. Group differences in coherence and in clock responses are tested by permutation with FDR over the gene set; t-values are never compared directly.
3. GSE22309 is reported both in full and restricted to same-run pairs; the restricted analysis is primary.
4. Fisher–Rao along a control is reported only against a θ-shuffle null; at present it does not pass and is removed from the results.
5. Repeat-biopsy artefact is stated for every within-person perturbation gene list.
6. Any longitudinal before/after comparison of dispersion is not reported (batch-confounded by design).

## E. What survives, in order of strength

1. Basal muscle transcriptome is blind to insulin sensitivity (A3).
2. The organs share one architecture and the individual one position (A1).
3. The disease does not switch basins (A2).
4. The healthy insulin response is coordinated and assembled over hours (A4).
5. Insulin resistance loses that coordination, including the clock coupling (B1, B2): consistent in three cohorts, significant in one, with a batch caveat; novel.
6. Interventions that restore magnitude do not detectably restore coordination (B3).

Lost: the prediabetic window, contraction, IGT criticality, blood, immune decoupling.
