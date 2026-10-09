# Contraction and decoupling: what the data say about the progression to type 2 diabetes

> **Replication status (islet, added after the first replication round).**
> F2 (contraction) does **not** replicate in two independent islet RNA-seq
> cohorts: GSE164416 (living donors; batch-balanced subset ND 14 / IGT 25 /
> T2D 34) gives dispersion 258 → 431 → 385 and GSE50244 (organ donors, HbA1c
> strata ND 39 / PreD 27 / T2D 11) gives 306 → 370 → 392; the mean pairwise
> distance in gene space also increases in both. The sign is reversed relative
> to the discovery islet cohort; confidence intervals include zero in all
> three. The contraction claim is therefore **withdrawn for islet** and held
> as *untested* for muscle and adipose until their replication cohorts are
> analysed. What replicates in islet: no second attractor (F5: GSE50244
> unimodal, ΔBIC −27; GSE164416 batch-balanced, one basin), and the Fisher
> peak inside the prediabetic range (GSE50244: HbA1c 5.86 %). GSE164416 shows
> a non-significant IGT maximum of I_c (2.06 vs 1.82 / 1.61, p = 0.19) and of
> dispersion, i.e. a weak critical-transition-like signature that the
> discovery cohort did not show; the full GSE164416 (with the IGT-only
> library batch) must not be used, as batch and stage are confounded. The
> theory below is kept as the discovery-cohort formulation; §2 consequence 1
> and §4 prediction 1 are the parts that now lack independent support.
> First remission test (GSE66306, monocytes, 19 paired, 3 months after
> surgery): dispersion 252 → 223, I_c 2.57 → 2.03, paired-permutation
> p = 0.63 — direction consistent with re-expansion of coherence (lower I_c)
> but underpowered and in blood cells, not tissue.


This document formulates the theory that the results support. It replaces the
pre-registered hypothesis of a saddle-type transition (`THEORY.md` §2, P1–P3),
which the data refute, and keeps what survives. Every statement below points
to a number in `results/`.

## 1. The five empirical facts

**F1. The organs share one transcriptional architecture, and the sharing decays
with progression.** Sheaf energy is 7–11 standard deviations below the
gene-correspondence null in every stage (z = −10.6 healthy, −7.5 IGT, −7.0
T2D; p = 0.005, B = 200) for islet, muscle and adipose tissue. The shared
component is strongest in health and weakens monotonically; stage-to-stage
differences are not individually significant (p = 0.15–0.53 at β = 4, 6, 8),
but the direction is the same in all three analyses and in the n → ∞
extrapolation (`results/sheaf*/main_*`).

**F2. The occupied region of state space contracts along its organised
directions.** The dispersion of patients in the covariate-adjusted embedding
(trace of the covariance over the first three principal components) falls
from healthy to T2D in all three tissues with power: islet 415 → 246 → 254,
muscle 156 → 133 → 94, adipose 514 → 422 → 175. The Gaussian differential
entropy falls accordingly (islet 11.0 → 10.7, muscle 10.2 → 9.4, adipose
11.3 → 9.7). The equal-n bootstrap of (dispersion_healthy − dispersion_T2D)
excludes zero in muscle (63, CI 10–118) and is positive with CI touching zero
in islet (154, CI −51–338) and adipose (321, CI −38–722). The contraction is
unchanged without covariate adjustment and without stage re-weighting.

In the full 800-gene space the mean pairwise distance also decreases, but
modestly: muscle 39.4 → 38.5 (−2 %), islet 38.6 → 37.0 (−4 %, minimum in IGT),
adipose 40.3 → 33.7 (−16 %). The difference is expected and is itself part of
the claim: in 800 dimensions the pairwise distance is dominated by isotropic
noise directions that are identical across stages, whereas the contraction
lives in the structured directions of variation. The funnel narrows along
the organised axes, not along noise. Liver (n = 5 Lean vs 9 Obese_T2D, an
obesity-stratified design) shows the opposite sign (−189, CI −354 to −5) and
is excluded from the conclusion (`results/landscape*/landscape_summary.tsv`,
columns `dispersion_*`, `entropy_*`, `meanpairdist_genes_*`, `contraction_*`).

**F3. Within-stage heterogeneity of gene–gene coupling rises.** The critical
transition index I_c (mean |gene–gene correlation| / mean |patient–patient
correlation|) increases monotonically: islet 1.55 → 1.74 → 2.18, muscle 2.50 →
2.63 → 2.93. It does not peak in IGT (p = 0.85 and 0.71 for "IGT maximal").
Together with F2 this means T2D patients are closer to each other along the
dominant axes of variation but less similar in their full gene profiles:
the disease state is contracted on the main axes and fragmented in the
residual ones.

**F4. Network organisation degrades at equal sample size.** With β fixed per
tissue and n equalised, global efficiency falls and effective resistance
rises from healthy to T2D in islet (EG 0.067 → 0.053; R̄ 0.69 → 0.77) and
adipose (EG 0.163 → 0.156; R̄ 0.068 → 0.071), with muscle flat. Raw metrics
without n equalisation would have shown the opposite (islet EG 0.038 vs
0.016 vs 0.022), driven entirely by group size (`results/summary/all_metrics_combined.tsv`).

**F5. There is no evidence of a second attractor, a saddle, or non-equilibrium
flux.** Persistent basins stable across bandwidths: one, in every tissue
with power. The two-component Gaussian mixtures that BIC favours do not
separate stages: Cramér's V between mixture component and stage is 0.12
(islet), 0.04 (muscle) and 0.25 in adipose where the second component is a
single healthy outlier (`gmm_min_component_n` = 1); components are spread
uniformly over batches.
IGT patients are not located at a pass: on the healthy → T2D axis they sit
on the healthy side (islet, mean −0.27) or in the middle (adipose, 0.67)
with a within-group spread larger than the healthy–T2D distance (SD 1.1–3.3).
The drift inferred by the Schrödinger bridge is at least as gradient-like as
under label permutation in all four tissues (p_gt_null = 0.43–0.90): no
detectable non-conservative component.

**F6 (gene level, exploratory).** The genes whose network position
reorganises most tissue-specifically in T2D are complement and innate-immune
genes: S100A4, C1QA, C1QB, C3AR1, ITGB2, LSP1, IL1R1, HLA-DOA, HPGDS. None
reaches FDR < 0.1 at the present n; the pattern is hypothesis-generating.
The IGT-specific list is unrelated (correlation of z-scores 0.15) and has no
such coherence.

## 2. The theory

**The healthy metabolic state is a broad, high-entropy, inter-organ-coherent
region of transcriptional state space; progression to type 2 diabetes is a
monotone contraction of that region accompanied by decoupling between organs
and stiffening of the coexpression networks. There is no intermediate
critical state and no second basin separated by a barrier.**

In landscape terms: not two valleys and a pass, but a single funnel whose
mouth narrows. The quasi-potential U = −ln P becomes steeper and its basin
smaller as glycaemia rises; the dynamics, as far as cross-sectional data can
tell, remain gradient (J ≈ 0). The entropy of the state distribution is the
order parameter, and it decreases.

Three consequences follow and each is supported:

1. *Loss of variability is loss of adaptability.* A broad healthy region
   means that many transcriptional configurations are compatible with
   normoglycaemia; the system has slack. Contraction removes that slack: the
   diabetic transcriptome is canalised into a narrow region it cannot easily
   leave. This is the transcriptional version of the loss-of-complexity
   paradigm of physiological ageing and disease (Lipsitz & Goldberger 1992),
   where reduced variability of heart rate, gait or glucose dynamics marks
   fragility. The "return barrier" of the original hypothesis is not a saddle
   to be crossed back but the small volume of the diabetic region: a
   perturbation large enough to re-expand it (bariatric surgery, intensive
   caloric restriction) is what remission requires.

2. *Decoupling precedes diagnosis.* Inter-organ coherence is already lower in
   IGT than in health (F1), and the Fisher–Rao information along a continuous
   glycaemic control peaks inside the prediabetic range (muscle, HbA1c 6.27 %;
   islet, fasting glucose 5.8 mmol/L; adipose 5.5 mmol/L): the distribution of
   transcriptional states changes fastest *before* the diagnostic thresholds
   (HbA1c 6.5 %, glucose 7.0 mmol/L). The IGT label, however, does not capture
   this: IGT patients are a heterogeneous, mostly healthy-side group. The
   transition is a property of the continuous control, not of the categorical
   stage.

3. *The diabetic state is one region with many residents.* Contraction on the
   main axes (F2) with fragmentation on the residual ones (F3) and
   tissue-specific immune reorganisation (F6) is what a single wide-but-rugged
   basin looks like from cohort data. It is compatible with the clinical
   subtypes of Ahlqvist et al. (2018) without requiring separate attractors
   for them.

## 3. Relation to the physiological bistability of T2D

The models of Topp et al. (2000) and Ha, Satin & Sherman (2016) are bistable
in the fast variables (glucose, insulin, β-cell mass). Our result does not
contradict them: the tissue transcriptome is a slow variable that integrates
over the glycaemic history. A fast bistable switch driving a slow canalising
response produces exactly a monotone, contracting slow variable with the
fastest change near the switch threshold (F2, consequence 2) and no
bimodality in the slow variable itself. The transcriptome records that the
switch has been thrown; it does not reproduce the switch. Testing the switch
requires fast variables or longitudinal sampling across the threshold.

## 4. What the theory predicts (testable, not yet tested)

- **Remission re-expands.** After bariatric surgery or diet-induced remission,
  tissue dispersion and inter-organ coherence should increase towards the
  healthy values. Paired pre/post adipose and muscle transcriptomes exist
  (GEO) and this is the cleanest next test.
- **Contraction is monotone in the control.** Along HbA1c (or fasting glucose)
  within a single cohort, dispersion in sliding windows should decrease
  monotonically, with the steepest change in the 5.7–6.4 % range.
- **No reversal of I_c.** With larger n per stage (≥ 60) the critical
  transition index should remain monotone; a peak in IGT would falsify the
  theory and revive the saddle hypothesis.
- **Decoupling is driven by immune genes.** The tissue-specific
  reorganisation in T2D should be enriched in complement/innate immunity in
  independent cohorts; removing those genes should reduce the loss of
  coherence more than removing a random gene set of equal size.
- **Order of tissues.** Because the sheaf localises incoherence per gene and
  per tissue, the theory predicts which organ decouples first for each gene
  module; with longitudinal data this becomes a causal ordering.

## 5. What is new relative to the literature

Dynamical-network-biomarker approaches (Chen et al. 2012) and population
optimal transport (GGOT, 2025) look for a critical state in a single tissue;
applied here, their statistic (I_c) finds none. The two measurable
quantities that do change monotonically — the volume of the occupied state
region and the inter-organ coherence — are not part of those frameworks.
The cellular-sheaf construction is what makes the second one measurable at
all, and the equal-n, configuration-null design is what prevents the first
one from being an artefact of cohort sizes (as the raw network metrics were).

## 6. Scope and limits, unchanged

Cross-sectional cohorts, one per tissue; correlation networks; n = 33–118;
liver underpowered. Nothing here is mechanistic. The contraction is measured
on the first three principal components of 800 high-variance genes; it
should be re-checked on the full gene space with a dispersion measure that
does not depend on a fixed embedding (e.g. mean pairwise distance on the
residualised expression), which `landscape.py` now also reports (`meanpairdist_genes_*`).
