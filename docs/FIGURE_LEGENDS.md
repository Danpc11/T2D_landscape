# Figure legends

Panels carry no titles; everything is described here. Every number below is written by
`run_replication.sh` into `results/` and read by `python/figures/make_figures.py`; regenerate both
before submission so that text and figures cannot drift apart.

## Statistical conventions used throughout

*Coherence* is a leave-one-out statistic: the mean cosine between each individual's displacement
(expression during insulin minus expression at baseline, in the principal subspace of the 800 most
variable genes, r = 5 components) and the mean displacement of the same group computed **without** that
individual. It is 1 when everyone moves the same way and 0 when directions are random, and it is
independent of how large the displacements are.

Group comparisons of coherence and of response direction are permutation tests that relabel individuals
between the two groups compared (2,000 permutations unless stated); the resulting null distribution is
drawn in the panel rather than summarised. Paired comparisons (before/after an intervention) permute the
two time points within each participant. Per-gene group × insulin interactions are permutation tests of
the difference in mean response between groups (5,000 permutations) with Benjamini–Hochberg FDR over the
nine clock-output genes tested. Differential expression at rest is a two-sided Welch *t*-test per gene on
expression residualised for the covariates named in the panel, with Benjamini–Hochberg FDR over all genes
tested. Gene-set statistics are the mean |statistic| of the set compared with 2,000 random sets of the
same size drawn from the same data; the grey band is the central 95% of that null and *P* is the
proportion of random sets that reach the observed value. Error bars are stated in each panel as s.d.,
s.e.m. or a 95% confidence interval. Exact *P* values are given to three decimals and reported as
*P* < 0.001 below that. No correction is applied across panels; each family of tests is corrected within
itself as stated.

---

## Fig. 1 | The metabolic organs share one transcriptional state, each person occupies one position in it, and that state drifts continuously with disease

**a**, A coexpression network of the kind the sheaf is built on: the 150 most variable genes of healthy
skeletal muscle (GSE25462, *n* = 15), with edges drawn for the strongest 3% of the biweight-midcorrelation
adjacency raised to the soft-thresholding power β = 6, node size proportional to connectivity and colour
given by hierarchical clustering into four modules. The sheaf does not use the drawing: it uses the
spectral subspace of this network's normalised Laplacian.

**b**,**c**, Network organisation by glycaemic stage in the four discovery cohorts, computed on the 800
genes shared by all organs at equal *n* per stage (subsampling, 20 repetitions; mean ± 95% CI). Global
efficiency $E_G$ (**b**) falls and effective resistance $\\bar{R}$ (**c**) rises from healthy to T2D in
islet (0.067 → 0.053 and 0.69 → 0.77) and adipose tissue (0.163 → 0.156 and 0.068 → 0.071); muscle is
flat and liver, with four samples per group, is uninformative. These metrics describe each organ
separately and motivate the cross-organ statistic that follows.

**d**, The cellular sheaf used throughout. Each organ's coexpression network contributes a stalk, the
spectral subspace of its normalised Laplacian (r = 8 components); edges carry orthogonal alignments
(restriction maps) fitted by Procrustes to a neutral reference; the sheaf energy *E* is the mean squared
Grassmann distance between aligned stalks, and is zero when the organs share one architecture.

**e**, Sheaf energy by glycaemic stage in the three discovery cohorts (islet GSE76895, muscle GSE18732,
adipose GSE27951; 800 shared high-variance genes, equal *n* per stage by subsampling, 20 repetitions).
Circles, observed; squares, mean of the gene-correspondence null, in which the gene labels of each tissue
are permuted independently so that the correspondence between organs is destroyed while every
within-organ property is preserved (*B* = 200 draws); error bars, 95% CI of the null. z = −10.6, −7.5 and
−7.0 for healthy, intermediate and T2D (*P* = 0.005 for each, the smallest value attainable with
*B* = 200). Stage-to-stage differences are not significant (*P* = 0.15–0.53).

**f**, The same statistic in GTEx v11 donors with skeletal muscle, subcutaneous adipose tissue and
pancreas from the same individual (*n* = 253). Expression was residualised for age, sex, Hardy
death-classification, RNA integrity number, ischaemic time and sequencing batch before the networks were
built. Grey, gene-correspondence null (*B* = 200); blue line, observed (E = 0.95 versus 1.91 ± 0.04 s.d.;
z = −26).

**g**, Position of each donor along the canonical variate shared by muscle and adipose tissue, expressed
as a rank and traced across the three organs (one line per donor, coloured by rank in muscle). Left,
observed; right, the same after shuffling donors between organs. First canonical correlation 0.93
(muscle–adipose), 0.91 (muscle–pancreas) and 0.89 (adipose–pancreas) against 0.23 ± 0.04 under donor
shuffling (*P* = 0.002, 500 permutations). In 19 living bariatric-surgery patients with paired
subcutaneous and omental biopsies (GSE20950) the same statistic is 0.93 against 0.59 ± 0.13
(*P* = 0.0005).

**h**, Effect of technical adjustment in GTEx: |z| of the sheaf energy (blue, left axis) and the
muscle–adipose canonical correlation (orange, right axis) before and after adding RNA integrity,
ischaemic time and sequencing batch to a model that already contained age, sex and Hardy classification
(z from −31.2 to −26.1; ρ 0.929 to 0.932).

**i**, Cross-validated R² of whole-blood expression predicting the same donor's organ state (ten blood
components predicting five organ components, five-fold cross-validation, *n* = 224 donors): 0.30
(muscle), 0.23 (adipose) and 0.23 (pancreas).

**j**, Number of basins of the quasi-potential *U* = −ln *p* above the bootstrap stability threshold τ in
each cohort, at three kernel bandwidths (0.7×, 1.0× and 1.4× the Silverman rule). Basins are counted by
sublevel-set persistence; τ is twice the standard deviation of the persistence of the leading basin over
200 bootstrap resamples. No cohort retains two basins across bandwidths.

**k**, Covariate-adjusted embedding of 77 islet donors (GSE50244) coloured by HbA1c stratum (ND < 5.7%,
prediabetes 5.7–6.4%, T2D ≥ 6.5%). Cramér's V between the component assigned by a two-component Gaussian
mixture and the glycaemic stratum is 0.17: the mixture favoured by BIC is unrelated to stage.

**l**, The genes that carry the shared individual position: Pearson *r* between each gene's residual
expression and the canonical variate shared by muscle and adipose tissue, for the 28 genes with the
largest effect of consistent sign in all three organs. These are the drivers identified by the present
approach; the canonical differential-expression drivers are in Fig. 2g,h.

---

## Fig. 2 | The disease is not visible at rest in any organ; it is visible in the response, and muscle is where that response fragments

**a**, Cohorts analysed per organ, at rest and with biopsies taken before and during a
hyperinsulinaemic–euglycaemic clamp. Ranges are *n* per stage (at rest) or *n* per group (clamp).

**b**, Sheaf energy between the three organs by stage, as in Fig. 1b. Coherence weakens monotonically
with stage (E = 1.75, 1.80, 1.81) but no stage-to-stage difference reaches significance (*P* = 0.15, 0.27
and 0.53 for intermediate versus healthy, T2D versus healthy and T2D versus intermediate; permutation of
stage labels, *B* = 500).

**c**, Critical-transition index I_c (mean |gene–gene correlation| divided by mean |donor–donor
correlation|, computed at equal *n* by subsampling) by stage, one line per cohort, coloured by organ. A
critical state would appear as a peak at the intermediate stage; none is present in any organ
(*P* ≥ 0.19 for "intermediate is the maximum" in every cohort, permutation of stage labels).

**d**, Ratio of state dispersion (trace of the covariance over the first three principal components of
the covariate-adjusted embedding) in T2D relative to healthy, with the 95% confidence interval of 2,000
bootstrap resamples at equal *n*; one point per cohort. Every interval crosses 1.

**e**, Under insulin the organs differ. Filled squares, ratio of mean response magnitude (resistant or
diabetic group divided by sensitive group); open circles, 1 + the difference in coherence between the
same two groups. Adipose tissue loses magnitude and keeps direction; muscle keeps magnitude and loses
direction. *n* per group is given beside each contrast; the underlying statistics are in Figs. 5 and 6.

**f**, Summary of what each organ contributed. Muscle is the only organ with paired insulin biopsies in
all three states (sensitive, resistant and diabetic).

**g**, The canonical read-out for comparison: genes differentially expressed at rest at FDR < 0.1
(two-sided Welch *t*-test per gene on expression residualised for the covariates of each study,
Benjamini–Hochberg over all genes tested). T2D versus normal glucose tolerance in islet (GSE164416,
*n* = 14/34: 1,079 of 14,077 genes; GSE50244, *n* = 39/11: 451 of 15,428) and in muscle (GSE25462,
*n* = 15/10: 2 of 48,496 probes); obese versus normal weight in adipose tissue (METSIM, *n* = 142/68:
7,332 of 13,718).

**h**, Set enrichment of the same contrasts: z of the mean |*t*| of each curated set relative to 2,000
random sets of equal size; asterisks mark *P* < 0.05. The canonical analysis recovers the expected
signatures (inflammation and complement, lipid and adipogenesis, ECM) in islet and adipose tissue and
almost nothing in muscle; none of these sets concerns the response to insulin, which is where the muscle
signal lies (Figs. 4 and 5). Set membership is listed in Methods and in
`python/analyses/09_classical_de.py`.

---

## Fig. 3 | The resting muscle transcriptome carries no information about insulin sensitivity

**a**, Design (GSE182120): resting vastus lateralis biopsies from 49 individuals profiled in two
laboratories, 24 with normal glucose tolerance and 25 with type 2 diabetes, all with insulin sensitivity
measured by hyperinsulinaemic–euglycaemic clamp (M-value 34.0 ± 15.8 versus 15.3 ± 9.3 mg kg⁻¹ min⁻¹,
mean ± s.d.). Expression was residualised on laboratory, age and BMI before correlating with the M-value.

**b**, Quantile–quantile plot of the genome-wide association between residual expression and the M-value
(Spearman's ρ, 21,595 genes). No gene reaches FDR < 0.1 (Benjamini–Hochberg); 14 genes reach *P* < 0.01
where ~215 are expected by chance. Known markers are highlighted (PPARGC1A ρ = 0.33, *P* = 0.021; PDK4
and TXNIP at *P* ≈ 0.05).

**c**, Volcano plot of T2D versus normal glucose tolerance at rest (two-sided Welch *t*-test per gene on
the same residuals). No gene reaches FDR < 0.1; the dashed line marks nominal *P* = 0.05.

**d**, Mean |ρ| with the M-value for curated gene sets (set size in parentheses) against 2,000 random
sets of equal size (grey bars, central 95% of the null). Only the oxidative programme separates from the
null (*P* = 0.012); every other set is indistinguishable from random (*P* > 0.1).

**e**, Spearman's ρ with the M-value for the clock-output genes examined in Figs. 4 and 5; all
*P* > 0.2. At rest these genes carry no information about insulin sensitivity, in contrast with their
behaviour during insulin (Fig. 5f).

Aggregate geometry is equally uninformative: dispersion 243 (NGT) versus 308 (T2D) and I_c 1.89 versus
1.59, neither significant; the distance of each individual from the healthy reference correlates weakly
with the M-value (ρ = −0.36, *P* = 0.010), a diffuse signal that no individual gene captures.

---

## Fig. 4 | The healthy response to insulin is coordinated across individuals, assembled over hours, and includes a reset of the clock output

**a**, Definition of the measure. Each person's displacement is the difference between the insulin and
the basal biopsy in the principal subspace of the 800 most variable genes; coherence is the mean cosine
of each displacement to the leave-one-out group mean.

**b**, Displacements of 20 insulin-sensitive individuals after 4 h of insulin (GSE22309) projected onto
the group mean direction (abscissa) and its principal orthogonal direction (ordinate). Coherence 0.77;
mean magnitude 16.3; net group response 13.1.

**c**, Coherence against time after the stimulus in four independent healthy datasets: 0.20 at 30 min
(GSE9105, *n* = 12), 0.24 at 1 h after a mixed meal (GSE231509, *n* = 7), 0.19 at 2 h (GSE7146, *n* = 6),
0.73 at 4 h (GSE9105, *n* = 12) and 0.77 at 4 h (GSE22309, *n* = 20). Coordination is absent early and is
assembled over hours.

**d**, Paired *t* per gene at 30 min versus 4 h in 12 healthy men (GSE9105). Blue, the 55 genes of the
replicated healthy programme, defined as |*t*| > 4 in the insulin-sensitive group of GSE22309 and the
same sign with |*t*| > 3 at 4 h in GSE9105.

**e**, Gene sets in the healthy 4-h response (GSE22309, insulin-sensitive group): mean |paired *t*| of
each set against 2,000 random sets of equal size (grey, central 95%). Immediate-early transcription
factors, chaperones, metallothioneins, canonical insulin targets and the clock output all separate from
the null (*P* < 0.05); ECM and the oxidative programme do not.

**f**, Paired *t* for clock-output genes at 30 min and 4 h in GSE9105 and at 4 h in the independent
GSE22309 cohort. Insulin represses DBP and NR1D2 and induces PER2 and BHLHE40, and does so only at 4 h.

**g**, Control for the repeat-biopsy artefact. Mean paired *t* of immediate-early genes, which respond to
the biopsy itself, and of clock-output genes at each time point: immediate-early genes respond at every
time point (0.5, 4.9 and 5.8), whereas the clock output responds only at 4 h (−0.1, 1.8 and 3.2).
Comparisons between groups are unaffected because every group underwent the same two-biopsy protocol.

---

## Fig. 5 | Insulin resistance preserves the magnitude of the response but fragments its direction and uncouples the clock

**a**, Magnitude of each individual's displacement by group (GSE22309; log scale; bars, mean ± s.d.):
16.3 (insulin-sensitive, *n* = 20), 18.5 (insulin-resistant, *n* = 20) and 13.5 (T2D, *n* = 15); no group
differs from the sensitive group (*P* > 0.4, permutation of group labels).

**b**, Loss of coherence relative to the insulin-sensitive group (circles: 0.77 − 0.35 = 0.42 for
insulin-resistant, 0.77 − 0.45 = 0.32 for T2D) against the distribution obtained by permuting group
labels between the two groups compared (violins, 2,000 permutations). *P* = 0.007 and *P* = 0.037.

**c**, The same comparison restricted to the 35 individuals whose two biopsies were hybridised in the
same batch (filled; 11, 13 and 11 per group) alongside all pairs (open). Effect sizes are preserved
(0.76, 0.45 and 0.46) but at this sample size the permutation tests give *P* = 0.10 and *P* = 0.19; the
direction contrast remains significant (cosine 0.14 against a null of 0.62, *P* = 0.023).

**d**, Direction of each person's response relative to the healthy mean direction, which points north, as
a rose of angles (bins of 20°; one count per individual, symmetrised about the vertical axis because only
the angle to the healthy direction is meaningful). Insulin-sensitive individuals cluster at north;
insulin-resistant and diabetic individuals spread around the circle. The cosine between the group mean
directions of the insulin-sensitive and T2D groups is 0.20, against 0.85 under permutation of group
labels (*P* = 0.001).

**e**, Coherence across all muscle datasets analysed, including two independent arms of a prediabetes
trial at baseline (GSE157988): 0.77 and 0.73 in healthy or insulin-sensitive groups, 0.35 in
insulin-resistant muscle, 0.42 and 0.02 in the two prediabetic arms and 0.45 in T2D; *n* beside each bar.

**f**, Response to insulin of each clock-output gene by group, shown as the absolute response on the
radius; dark sectors are repressed genes, light sectors induced. Asterisks, FDR < 0.05 for the
group × insulin interaction versus the insulin-sensitive group; daggers, FDR < 0.1 for the T2D
comparison (permutation of group labels, 5,000 permutations, Benjamini–Hochberg over nine clock genes).
Insulin-resistant versus sensitive: DBP *P* = 0.003, FDR = 0.013; PER2 *P* = 0.004, FDR = 0.013; NR1D2
*P* = 0.012, FDR = 0.024; BHLHE40 *P* = 0.047, FDR = 0.070. T2D versus sensitive: DBP *P* = 0.024, TEF
*P* = 0.045, HLF *P* = 0.041, NR1D2 *P* = 0.022 and BHLHE40 *P* = 0.045, all at FDR = 0.054. The
response of BHLHE40 persists in all three groups.

**g**, Effect of chronic insulin on the same genes in primary myotubes from donors with normal glucose
tolerance (*n* = 7) and with type 2 diabetes (*n* = 5) (GSE182117); points are the paired *t* across
donors of the difference between treated and control arms averaged over sampling times. DBP is repressed
in NGT cells (*t* = −3.3) more than in T2D cells (*t* = −1.7); BHLHE40 responds in both.

**h**, Amplitude of the 24-h oscillation of clock-output genes in the same myotubes, from a cosinor fit
per donor in the control arm; mean ± s.d. Amplitudes are ~30% lower in T2D cells (for example DBP
1.03 ± 0.31 versus 0.71 ± 0.18), a difference that does not reach significance at this sample size
(*P* ≈ 0.15, two-sided Mann–Whitney).

**i**, Fate of the 55-gene healthy programme in insulin-resistant and diabetic muscle: kept (same sign
and |*t*| > 2), lost (|*t*| < 2) or inverted (opposite sign and |*t*| > 2). A third of the programme is
lost (38% in insulin-resistant, 42% in T2D muscle) and essentially none of it is inverted.

---

## Fig. 6 | Adipose tissue fails by magnitude, and interventions restore magnitude without restoring coordination

**a**, Magnitude (grey, left axis) and coherence (blue, right axis) of the adipose response to a clamp in
23 non-obese individuals, 23 obese individuals, and the same 23 obese women two years after bariatric
surgery (Rydén et al. 2016; CAGE sequencing). Magnitude 7.4, 5.3 and 6.8; coherence 0.37, 0.19 and 0.14.
The loss of coherence with obesity does not reach significance at this sample size (*P* = 0.11,
permutation of group labels), and surgery does not restore it (*P* = 0.71, paired permutation of the two
time points within each woman).

**b**, The same quantities in an independent adipose clamp cohort (GSE26637; 5 insulin-sensitive and 5
insulin-resistant individuals): magnitude 28.4 versus 15.6, coherence 0.83 versus 0.63. The resistant
group keeps the direction of the sensitive group (cosine 0.95).

**c**, Magnitude ratio (resistant or obese divided by sensitive or non-obese) in the two adipose cohorts
and in muscle, on one scale: 0.71 and 0.55 in adipose tissue, 1.14 in muscle.

**d**, Coherence before and after 10 weeks of nicotinamide mononucleotide or placebo in prediabetic women
(GSE157988; 11 and 12 participants). Neither arm changes (*P* = 0.77 and *P* = 0.82, paired permutation
of the two time points within each participant, 3,000 permutations), although the trial improved
clamp-measured insulin sensitivity in the NMN arm.

**e**, Change in coherence detectable with 80% power as a function of group size, for a two-sided test at
α = 0.05 assuming a between-individual s.d. of 0.3 (the value observed across these cohorts); dashed
lines mark the two intervention cohorts available. At *n* = 11–23 a restoration smaller than ~0.25 would
be missed, so the absence of restoration in **a** and **d** is an absence of evidence rather than
evidence of absence.

**f**, Summary: healthy organs respond along one shared direction; insulin-resistant muscle responds in
scattered directions, and insulin-resistant adipose tissue responds in the same direction with half the
magnitude.
