# Insulin resistance reduces the directional concentration of the human muscle transcriptional response to insulin

**Authors.** [to be completed]

**Affiliations.** [to be completed]

**Correspondence.** [to be completed]

**Article type.** Article (*Nature Metabolism*)

**Word count.** Main text ~4,700; six main figures; three Extended Data figures.

\newpage

# Abstract

Type 2 diabetes is defined by what insulin fails to do, yet most molecular studies profile tissues at rest. We reanalysed 24 public human transcriptomic studies across 29 accessions and 1,689 participants, spanning pancreatic islet, skeletal muscle, adipose tissue and blood, including 253 GTEx donors with three organs sampled from the same individual and 140 participants biopsied before and during a defined perturbation. In resting muscle from 49 clamp-phenotyped individuals, no gene was associated with insulin sensitivity at a false discovery rate below 0.1, and the previously reported association in this cohort did not survive adjustment for age or adiposity. The response to insulin was informative. Treating each person's response as a direction on a sphere, healthy individuals concentrated around a shared direction whereas insulin-resistant and diabetic individuals did not: per-person alignment fell by 0.49 (P = 0.006 within hybridisation batch) and the von Mises–Fisher concentration κ fell from 8.2 to 3.6. In a hierarchical model across four cohorts and 140 participants, metabolic impairment reduced κ by 44% under insulin, and this reduction was abolished under acute exercise (interaction P = 0.005). Insulin resistance therefore reduces the directional concentration of the muscle transcriptional response to insulin specifically, a property invisible at rest and measurable in any study that biopsies before and during a perturbation.

\newpage

# Introduction

Skeletal muscle disposes of most postprandial glucose and its resistance to insulin is among the earliest defects on the path to type 2 diabetes. The molecular description of that defect has been sought mainly in two places: in insulin signalling, where impaired Akt phosphorylation is reproducible but accompanied by heterogeneous downstream effects, and in the resting transcriptome, where two decades of profiling have produced a short list of consistent findings, essentially the reduced oxidative programme marked by PGC-1α, alongside many signatures that have not replicated.

A third place has been explored less. Insulin regulates several hundred muscle genes within hours, and specific insulin-responsive genes are blunted in obesity and diabetes. Insulin is also a timing signal: it resets the clock of adipose tissue and of isolated adipocytes and of peripheral tissues in rodents. The muscle clock in turn governs insulin sensitivity, circadian misalignment lowers insulin sensitivity in humans, and muscle cells from people with type 2 diabetes show intrinsically disrupted circadian oscillations. Recent proteomic work has shown that the molecular signatures of insulin resistance are highly individual, raising the possibility that what differs between insulin-resistant people is not a shared defect but the loss of a shared response.

Two obstacles have kept this question open. The first is measurement. Comparing responses between groups requires separating *how far* a tissue moves from *where* it moves to, because a group whose members all respond strongly but in different directions produces the same average as a group that barely responds; conventional differential expression conflates the two. Recent work in single-cell CRISPR screens has introduced a directional coherence statistic and concluded that it is largely a restatement of effect magnitude, with correlations of 0.84 to 0.98 between cells. Whether the same holds between people, where genetic background, age and adiposity differ, is unknown and is not a question single-cell data can answer.

The second obstacle is design. The informative experiment is a perturbation with paired biopsies in the same person, which exists in public repositories but in small, scattered cohorts that have not been analysed together.

Here we ask whether the information about insulin sensitivity in human muscle lies in how the tissue *is* or in how it *responds*. We first establish what the resting state does and does not show (Figs. 1–3). We then treat each person's response as a direction, model those directions explicitly, and test whether what we find is specific to insulin by applying the same model to a different physiological perturbation in the same tissue (Figs. 4–6).

# Results

## Tissues of the same person share transcriptional architecture, and this is not specific to metabolic organs

We represented each tissue's coexpression network by the spectral subspace of its normalised Laplacian and compared tissues after orthogonal alignment, so that organs expressing different genes can be compared through the organisation of their networks rather than their expression levels (Fig. 1a,d; Methods). A null that permutes the gene labels of each tissue independently destroys the correspondence between organs while preserving every within-organ property.

In GTEx v11 donors with skeletal muscle, subcutaneous adipose tissue and pancreas sampled from the same individual (n = 253), and after adjusting for age, sex, death classification, RNA integrity, ischaemic time and batch, the residual distance was 0.953 against a null of 1.911 ± 0.032 (z = −29.7, 200 permutations; Fig. 1f). The first canonical correlation between a donor's scores in muscle and adipose tissue was 0.93, and 0.91 and 0.89 for the two other pairs, against 0.23 under donor shuffling (Fig. 1g). The same held in 19 living bariatric-surgery patients with paired adipose depots (0.93 versus 0.59), indicating that the effect is not an artefact of post-mortem sampling, and whole blood predicted 22 to 28% of the principal variation of each organ in leakage-free cross-validation (Fig. 1i).

This shared individual position is, however, **not specific to the organs of metabolism**. Extending the analysis to eleven GTEx tissues, all 30 available tissue pairs showed canonical correlations far above the shuffled null, including muscle–ileum (0.87), adipose–kidney (0.87) and liver–adrenal (0.89). Metabolic pairs showed it more strongly (excess over null 0.62 versus 0.47; Mann–Whitney P = 0.001; +0.126 ± 0.045 adjusting for the number of donors), and replacing pancreas with stomach in the three-organ comparison raised the energy from 0.95 to 1.21 and lowered z from −29 to −22. We therefore report this as a general property of tissues from the same person, graded rather than exclusive, and treat it as context for what follows rather than as a metabolic finding.

Nor does the resting organisation separate disease stages. In eight cohorts, three tissues and two platforms, the quasi-potential had one basin, with no cohort retaining a second across kernel bandwidths (Fig. 1j), and where a mixture model favoured two components they were unrelated to glycaemic stage (Fig. 1k). Recomputed at equal sample size per stage, the sheaf energy was 0.861, 0.848 and 0.790 for healthy, intermediate and diabetic (z = −6.9, −6.6, −7.4), showing shared architecture at every stage with **no monotonic change** (Figs. 1e, 2b). The critical-transition index shows no peak at the intermediate stage in any organ (Fig. 2c) and the dispersion ratio crosses one in every cohort (Fig. 2d).

## In resting muscle, no gene is associated with insulin sensitivity once age and adiposity are accounted for

We used the largest public dataset with both a resting transcriptome and a clamp-measured phenotype: 49 individuals profiled in two laboratories, all with an M-value (Fig. 3a).

Correlating basal expression with the M-value without covariate adjustment, as previously reported for this cohort, returns 739 genes at FDR < 0.1 in our pipeline, with 85 of 211 mitochondrial genes among them, reproducing the published mitochondrial enrichment. The association does not survive adjustment. In this cohort the M-value is strongly correlated with body-mass index (ρ = −0.58, P = 1 × 10⁻⁵) and with age (ρ = −0.39, P = 0.005), and adjusting expression for either one alone removes every gene at FDR < 0.1 (Fig. 3b and Extended Data Fig. 1). Notably, body-mass index by itself has no genes at that threshold, so the signal arises from variance shared between adiposity, age and insulin sensitivity rather than from a separate adiposity signature.

In the fully adjusted analysis, using partial correlation so that the M-value is residualised on the same covariates as expression, no gene among 21,595 reaches FDR < 0.1 and 231 reach P < 0.01 where 215 are expected by chance. The comparison of diabetic with normoglycaemic individuals at rest gives no gene either, and with Welch's test the independent muscle cohort GSE25462 gives none (Fig. 2g). Of the curated gene sets, only the oxidative programme separates from size-matched random sets (Fig. 3d), and the clock-output genes that become central below show no association at rest (Fig. 3e).

Because age and adiposity are themselves determinants of insulin sensitivity, a cross-sectional design cannot separate a transcriptional correlate of sensitivity from a correlate of its determinants. This is a limitation of the resting design, not of this cohort, and it motivates the shift from state to response.

## The healthy response to insulin is concentrated around a shared direction and is established over hours

To compare responses between groups without choosing genes in advance, we represented each biopsy by its position in the principal subspace of the 800 most variable genes and each person's response as the displacement between their paired biopsies (Fig. 4a). Two properties are then independent of magnitude: the *direction* each person moves in, and how concentrated those directions are.

We quantify concentration in three complementary ways. The *coherence* of a group is the mean cosine between each displacement and the leave-one-out group mean; the leave-one-out construction is essential, since the naive statistic returns 0.41 at n = 4 and 0.19 at n = 20 when no direction is shared at all, whereas the leave-one-out version returns zero at every group size (Extended Data Fig. 3a). Because coherence collapses n people into a single number per group, our primary analysis instead uses the *alignment of each participant* with the reference direction, which preserves one observation per person. Finally, we fit the generative model these statistics correspond to: a von Mises–Fisher distribution over the sphere with mean direction μ and concentration κ.

In 20 insulin-sensitive individuals, coherence after 4 h of insulin was 0.77 and κ was 8.2 (Fig. 4b). Coherence was lower in cohorts sampled earlier, 0.20 at 30 min, 0.24 at 1 h after a mixed meal and 0.19 at 2 h, and higher at 4 h in two independent cohorts (Fig. 4c); these are separate cohorts rather than a within-person time course, but the same ordering appears gene by gene within the single cohort sampled at two times (Fig. 4d) and, as shown below, in an entirely different perturbation.

The genes carrying the shared 4-h direction, replicated across two healthy cohorts, are not the canonical insulin-signalling genes alone but a broader programme: immediate-early and stress transcription factors, metallothioneins and HMOX1, chaperones, the amino-acid transporter SLC7A5, and BHLHE40 (Fig. 4e). Insulin repressed the clock-output genes DBP and NR1D2 and induced PER2, and did so only at 4 h (Fig. 4f). Because a second biopsy taken near the first induces immediate-early genes irrespective of the stimulus, we verified that this cannot explain the clock finding: immediate-early genes respond at every sampling time whereas clock-output genes respond only at 4 h (Fig. 4g).

## Insulin resistance preserves magnitude but reduces directional concentration

In the same cohort, insulin-resistant and diabetic individuals responded with magnitudes indistinguishable from the sensitive group; a two-one-sided-test with a margin of half a log2 unit does not establish equivalence either (P = 0.26 and 0.60), so we report this as no detected reduction rather than as preserved magnitude (Fig. 5a).

What differed was the direction. Per-person alignment with the healthy response direction fell from 0.75 to 0.31 in insulin-resistant and 0.26 in diabetic muscle (difference 0.44, P = 0.0017; and 0.49, P = 0.0010; permutation with full re-estimation of the reference direction in each permutation). Hybridisation batch is unevenly distributed across groups in this 2007 dataset, and we show that batch alone can produce coherences of 0.75 to 0.87 between people who share nothing but a hybridisation run (Extended Data Fig. 3b). We therefore treat the analysis restricted to the 35 individuals whose two biopsies were processed in the same batch as primary. There the effect is undiminished: differences of 0.51 (P = 0.013) and 0.47 (P = 0.006), and the decisive control, pseudo-groups defined by hybridisation run alone, gives −0.10 (P = 0.57) (Fig. 5b).

The von Mises–Fisher fit puts a parameter on this: κ = 8.2 in the sensitive group against 3.6 in the insulin-resistant and 4.0 in the diabetic group, a 2.3-fold loss of concentration (likelihood-ratio P = 0.0016 and 0.0086; 0.0054 and 0.061 within batch; Fig. 5c). The bootstrap intervals for κ are wide and partly overlapping (4.7–21.0 versus 2.3–6.7), so κ itself is imprecisely estimated at this sample size even though the paired comparison is significant.

Crucially, concentration is not a restatement of magnitude in these data. Across 29 groups the correlation between magnitude and coherence was 0.33 (P = 0.077), and within the key cohort it is inverted: the insulin-resistant group has the largest mean displacement (18.8) and nearly the lowest alignment (0.31). This differs from single-cell CRISPR screens, where the two are tightly coupled (ρ = 0.84–0.98), and the difference is biologically sensible, since variation between people has sources that isogenic cells do not have.

The canonical insulin programme was intact in these individuals: TXNIP, KLF15 and IRS2 were repressed and PPP1R3B and HES1 induced in every participant of an independent prediabetic cohort. Of the 55-gene healthy programme, a third was lost in insulin-resistant and diabetic muscle and essentially none was inverted (Fig. 5i). Among clock-output genes, tested as a pre-specified set of nine, the group-by-insulin interaction was significant for DBP, PER2 and NR1D2 in the insulin-resistant comparison and borderline for five genes in the diabetic comparison, whereas BHLHE40 responded in all three groups (Fig. 5f). In primary myotubes from donors with and without diabetes exposed to high glucose plus insulin, DBP was repressed in cells from normoglycaemic donors more than in diabetic cells, and the amplitude of clock-output oscillations was lower in diabetic cells by a factor of 0.71 on author-independent processing (Wilcoxon across clock genes, P = 0.008) and 0.82 on the authors' batch-corrected matrix (P = 0.055) (Fig. 5g,h). Resting expression of these genes carries no information about insulin sensitivity (Fig. 3e), so this is a property of the response.

No single gene survives genome-wide correction for the group-by-insulin interaction: two genes at FDR < 0.05 for the insulin-resistant comparison and one for the diabetic comparison, with 95 genes at P < 0.001 where 8 are expected. This excess without individual significance is the expected signature of a distributed effect and is what the directional analysis is designed to capture.

## The loss of concentration is specific to insulin

A reduction in concentration could mean that diabetic muscle has lost the ability to respond in concert to anything. We tested this with a different physiological perturbation in the same tissue and in participants recruited by the same group: an acute bout of exercise with biopsies at rest and three hours into recovery, in 20 men with type 2 diabetes and 17 with normal glucose tolerance, and the same protocol sampled in adipose tissue in 23 and 16 participants.

Taken as separate comparisons these are reassuring but not significant: per-person alignment differs by −0.12 (P = 0.41) in muscle and −0.18 (P = 0.30) in adipose, in both cases favouring the diabetic group. The question is better answered jointly. We fitted a hierarchical von Mises–Fisher model to all 140 participants from four cohorts, letting each cohort keep its own mean direction, with log κ = b₀ + b_cohort + b₁·impairment + b₂·(impairment × exercise).

Metabolic impairment reduced κ by a factor of 0.56 under insulin (b₁ = −0.585, 95% CI −1.31 to +0.01; likelihood ratio P = 0.014), and the interaction with exercise was positive and significant (b₂ = +0.864, 95% CI +0.13 to +2.15; likelihood ratio **P = 0.005**), so that under exercise the effect of impairment is abolished (κ × 1.32). The specificity that is not significant as a pair of separate tests is significant as an interaction estimated on 140 people.

Two further observations follow. Immediately after exercise the response is idiosyncratic in both groups and becomes concentrated only at three hours of recovery, reproducing in a different stimulus the time dependence seen with insulin. And because the exercise response is intact, the reduced concentration of the insulin response cannot be attributed to a general transcriptional failure of diabetic muscle, to biopsy trauma, or to the clinical state of the participants.

Concentration is also not a generic property of any perturbation. Across eleven contrasts, values near 0.8 occur for acute stimuli in muscle (insulin 4 h, 0.77; acute exercise, 0.80 in an independent cohort at P = 0.0005) whereas slow perturbations produce individual rather than shared trajectories: weeks of training give 0.31 in muscle and −0.09 in adipose tissue **in the same participants**, a diurnal contrast gives 0.10, and ten days of cold acclimation give −0.01 (Supplementary Note 3).

## Adipose tissue fails differently, and interventions restore magnitude without demonstrably restoring concentration

Adipose tissue sampled before and during a clamp showed the opposite pattern to muscle. Obese and insulin-resistant adipose tissue responded along the healthy direction, with cosines of 0.88 and 0.95, at approximately half the magnitude (Fig. 6a–c). Two years after bariatric surgery in the same 23 women, the magnitude of the response had returned to the level of non-obese controls while coherence had not (Fig. 6a), and ten weeks of nicotinamide mononucleotide, which improved clamp-measured insulin sensitivity, did not change it either (Fig. 6d). Simulation on these data shows that with 11 to 23 participants the power to detect a restoration affecting a third of individuals is 0.3 to 0.4 (Extended Data Fig. 3c), so these are absences of evidence. Adipose tissue also shows lower concentration than muscle in healthy people and for an acute stimulus, which qualifies the two-failure-modes summary.

# Discussion

Four statements summarise this work. The resting transcriptome of human muscle is uninformative about insulin sensitivity once age and adiposity are accounted for. The healthy response to insulin is concentrated around a shared direction, is established over hours, and includes late modulation of clock-output genes. Insulin resistance preserves the canonical signalling programme and shows no detectable reduction in magnitude, but reduces the directional concentration of the response. And that reduction is specific to insulin: in a joint model of 140 participants, the same impairment does not reduce concentration under exercise.

The first statement helps explain an old frustration. Signatures of insulin resistance in resting muscle have rarely replicated beyond the oxidative programme. Our reanalysis adds a specific reason: in the one cohort with a gold-standard phenotype, the association with insulin sensitivity is not separable from age and adiposity, and vanishes when either is adjusted for. Any cross-sectional resting design faces this, because the covariates are themselves determinants of the phenotype. If the disease lives in the response, resting profiles of any size will keep finding little.

The measurement framework matters as much as the result. Directional coherence has been introduced for single-cell CRISPR screens, where it proved to be largely a restatement of effect magnitude. Between people it is not: the correlation is 0.33 and inverts within the key cohort. Two methodological points follow for anyone applying such a statistic to human cohorts. The naive form is badly biased at the sample sizes typical of biopsy studies, returning 0.41 at n = 4 when no direction is shared, so the leave-one-out form is required. And a group-level statistic discards the individual observations that make the comparison powerful: using per-person alignment instead of group coherence changed P from 0.10 to 0.013 on identical data.

The clock results connect literatures that have developed separately. We observe the corresponding step in vivo: modulation of clock output is part of the normal 4-h response of human muscle and is altered in insulin resistance, while BHLHE40 keeps responding. DBP, TEF and HLF are direct BMAL1 targets and the PAR-bZip output through which the clock times metabolism; their altered response, with resting expression unchanged, points to the coupling rather than to the oscillator. We stop short of causal language, and we note that this was a pre-specified test of nine genes, not a discovery: genome-wide, these genes do not stand out.

That no single gene survives correction is a prediction of the model rather than a shortcoming of the data. If what is lost is the agreement between people about which direction to move in, the effect is distributed by construction, and the 95 genes at P < 0.001 where 8 are expected is what that looks like. It also means the data do not nominate a drug target, and we do not propose one. What they support is a statement about strategy rather than molecule: if what fails is a response that takes hours to assemble and involves clock output, interventions acting on the timing of metabolic signals are the class worth testing, and the concentration of the response is an endpoint current trials do not measure.

## Limitations of the study

The insulin comparison rests on a single clamp cohort from 2007. We show that hybridisation batch can by itself generate apparent coordination as high as that we attribute to health, which is why the same-batch analysis is primary; there the effect is undiminished and the batch-only control is null, but the cohort remains one cohort and replication in a modern sequenced clamp study with paired biopsies is the single most valuable next experiment. The deposited data of the largest recent clamp cohort contain only fasted biopsies, so that replication requires new data. Two sequential biopsies induce an immediate-early response independent of the stimulus; this does not affect comparisons between groups but inflates the apparent insulin specificity of the healthy gene list. The time course in Fig. 4c combines independent cohorts. The hierarchical model estimates the interaction well but the main effect of impairment only imprecisely (95% CI −1.31 to +0.01), and it assumes a common concentration parameter structure across tissues and stimuli that differ in many ways besides the one modelled. The adipose clamp cohorts are small and the myotube comparison has seven and five donors. GTEx donors are post mortem. Pancreatic islet could not be included in the main response analyses because no in vivo perturbation with paired sampling exists for that organ in humans. All cohorts are cross-sectional at the group level, and no claim here is causal.

# Methods

Methods are written purpose first: each subsection states what question the procedure was designed to answer before describing how it was implemented. Code implementing every step, with fixed seeds, is available as described under Code availability; each numbered script corresponds to one subsection.

## Cohorts and preprocessing

*Purpose.* To assemble every public human dataset that could address either the resting state or the response, without selecting cohorts on their results.

Twenty-four independent studies across 29 accessions and 1,689 participants were used, after removing donors shared between accessions of the same study; cohort-by-cohort details, platforms and group sizes are in Supplementary Table 1. Cohorts were included if they provided human transcriptomes with a glycaemic or insulin-sensitivity phenotype, or paired sampling around a perturbation; none was excluded after analysis. Processed expression was used as deposited where available; counts were converted to log2 counts per million with an expression filter; probes were mapped to symbols with the platform annotation. Where a technical factor was unevenly distributed across groups it was either regressed out or used to define a restricted analysis, as stated per analysis.

## Networks, the cellular sheaf, and the specificity of the shared individual position

*Purpose.* To ask whether organs expressing different genes are organised alike and whether an individual occupies a consistent position across their tissues, and then to test whether that is a property of metabolic organs or of tissues in general.

Each organ contributes the leading spectral subspace of its normalised Laplacian, built on biweight-midcorrelation adjacency at a soft-thresholding power fixed once in the healthy state. Organs are compared after orthogonal alignment by Procrustes to a neutral reference, because the eigenbasis is defined only up to rotation; the energy is the mean squared distance between aligned subspaces. The gene-correspondence null permutes gene labels within each organ, destroying only the correspondence between organs; the donor-shuffle null permutes individuals between organs. All permutation draws are exported. The specificity analysis repeats the same computation for five metabolic and four non-metabolic GTEx tissues, and for all 30 tissue pairs with sufficient donors, comparing the excess over the null between metabolic and other pairs with and without adjustment for the number of donors.

Blood-to-tissue prediction fits residualisation, feature selection, scaling and both principal-component decompositions **inside each training fold only**.

## Quasi-potential landscape

*Purpose.* To test whether health and diabetes are two stable states separated by a barrier.

The density of covariate-adjusted, stage-balanced positions is converted to a quasi-potential by negative logarithm; basins are counted by sublevel-set persistence with a bootstrap stability threshold and must survive three kernel bandwidths. Gaussian mixture models and Cramér's V between component and stage provide an independent check.

## Geometry of the within-person response

*Purpose.* To compare responses between groups in a way that distinguishes a group that responds weakly from one that responds strongly but inconsistently.

Samples are projected into the principal subspace of the 800 most variable genes (five components) and each person's response is the displacement between their paired biopsies. Three quantities are reported. *Coherence* is the mean cosine between each displacement and the mean of the person's own group computed without that individual. *Alignment* is the same cosine evaluated per person against the reference group's leave-one-out mean, and is the primary outcome, because coherence collapses n individuals into one number and so compares two observations rather than n. *Magnitude* is the Euclidean norm.

Group comparisons of alignment use permutation of individual labels with **full re-estimation of the reference direction in each permutation**, since that direction is defined by the reference group; 3,000 permutations. Comparisons before and after an intervention permute the two time points within each participant. Per-gene group-by-stimulus interactions use permutation with Benjamini–Hochberg correction over the nine pre-specified clock-output genes, and are reported as a targeted test, not a screen; the genome-wide interaction test is reported separately.

## Von Mises–Fisher model and the hierarchical pooled analysis

*Purpose.* To replace a descriptive statistic with a generative model whose parameter is interpretable, comparable across cohorts of different size, and testable by likelihood ratio; and to combine cohorts using participants rather than summary statistics.

Each person's response direction is taken as a unit vector and modelled as a draw from a von Mises–Fisher distribution with mean direction μ and concentration κ, fitted by maximum likelihood with the normalising constant evaluated through the exponentially scaled modified Bessel function. Groups are compared by a likelihood-ratio test against the restricted model in which both share κ but keep separate directions; confidence intervals for κ come from 400 bootstrap resamples.

The hierarchical model pools all paired cohorts. Each cohort keeps its own mean direction, because tissue, stimulus and platform differ, and the concentration is modelled as log κ = b₀ + b_cohort + b₁·impairment + b₂·(impairment × exercise). The parameter of interest is b₁, the change in concentration associated with metabolic impairment, and b₂, which tests whether that change is specific to insulin. Significance is by likelihood ratio against the nested models; confidence intervals come from bootstrap resampling of participants stratified by cohort.

## Confounding, bias and power

*Purpose.* To establish how much apparent concentration the measures can produce in the absence of a shared direction, how much a technical factor can generate, and what the available cohorts could detect.

Three simulations run on the GSE22309 data. First, displacements are replaced by isotropic directions of the same length, which destroys any shared direction without preserving the original axes as a sign flip would, and the coherence is recomputed at every group size for both the naive and the leave-one-out statistic. Second, pseudo-groups are defined by hybridisation run alone, ignoring clinical status; and the clinical comparison is repeated restricted to individuals whose two biopsies were processed in the same run. Third, pairs of groups are drawn from the sensitive displacements with a fixed fraction of one group randomised in direction, and the proportion of significant permutation tests gives the power.

## Differential expression and gene sets

*Purpose.* To provide the canonical read-out against which the response-based analysis is compared.

Differential expression at rest uses a two-sided Welch *t*-test per gene on expression residualised for the covariates of each study, with Benjamini–Hochberg correction. The contrast of each cohort is declared explicitly rather than inferred. Association with insulin sensitivity uses partial Spearman correlation, residualising the M-value on the same covariates as expression; the ladder of adjustments from unadjusted to fully adjusted is reported. Gene sets are curated rather than taken from a database and enrichment compares the mean absolute statistic with 2,000 random sets of the same size drawn from the same data. Transcription-factor activity is inferred from DoRothEA regulons of confidence A–C as the mean target response weighted by mode of regulation, z-scored against size-matched random target sets.

## Reporting

Exact sample sizes, statistics and *P* values are given in the figure legends. Random seeds are fixed. Analyses were run with Python 3.10 and R 4.3; versions are pinned in the repository.

# Data availability

All data analysed are public. GEO accessions are listed in Supplementary Table 1. GTEx v11 expression and open-access annotations are from the GTEx portal. The adipose clamp data of Rydén et al. are from the export accompanying that publication. No new data were generated.

# Code availability

All code is at https://github.com/Danpc11/t2d_landscape, including cohort preparation, the network, sheaf and landscape implementations, the response-geometry statistics, the von Mises–Fisher and hierarchical models, the specificity, confounding and power analyses, and the scripts that generate every figure from the resulting tables.
