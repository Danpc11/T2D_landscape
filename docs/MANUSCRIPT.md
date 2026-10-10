# A coordinated transcriptional response to insulin is a tissue-level property of human skeletal muscle and is lost in insulin resistance

**Authors.** [to be completed]

**Affiliations.** [to be completed]

**Correspondence.** [to be completed]

**Article type.** Article (*Nature Metabolism*)

\newpage

# Abstract

Type 2 diabetes is defined by what insulin fails to do, yet molecular studies of skeletal muscle almost always profile the tissue at rest. We reanalysed 24 public human studies across 29 accessions and 1,689 participants, including 140 individuals biopsied before and during a defined perturbation. In resting muscle from 49 clamp-phenotyped individuals no gene was associated with insulin sensitivity, and the association previously reported in that cohort disappeared once age or adiposity was accounted for. The response to insulin was informative. Treating each person's response as a direction in expression space, healthy individuals concentrated around a shared direction whereas insulin-resistant and diabetic individuals did not: per-person alignment fell by 0.49 and the von Mises–Fisher concentration of response directions fell from 8.2 to 3.6. In a hierarchical model across four cohorts, metabolic impairment reduced concentration by 44% under insulin and not at all under acute exercise. The coordinated response is a property of the organised tissue rather than of the myocyte: it is accounted for by the cellular composition of the biopsy, and in primary myotubes from 24 donors stimulated with insulin, coordination was weak in everyone and did not differ between healthy and diabetic donors. We propose that insulin elicits a coordinated transcriptional response only in intact muscle, and that insulin resistance degrades that tissue-level coordination.

\newpage

# Introduction

Skeletal muscle disposes of most postprandial glucose and its resistance to insulin is among the earliest defects on the path to type 2 diabetes. Two decades of transcriptional profiling of resting muscle have produced a short list of reproducible findings, essentially the reduced oxidative programme marked by PGC-1α, alongside many signatures that have not replicated. Proteomic work has recently shown that the molecular signatures of insulin resistance are highly individual, raising a possibility that resting profiles cannot address: that what differs between insulin-resistant people is not a shared defect but the loss of a shared response.

Testing this requires a different measurement and a different design. The measurement must separate *how far* a tissue moves from *where* it moves to, because a group whose members all respond strongly but in different directions produces the same average as a group that barely responds; conventional differential expression conflates the two. The design must be a perturbation with paired biopsies in the same person. Such cohorts exist in public repositories but are small and have never been analysed together.

Here we assemble them. We first establish what the resting state does and does not show, then describe the geometry of the within-person response to insulin, show that its coordination is reduced in insulin resistance and that this reduction is specific to insulin, and finally localise the phenomenon: the coordinated response is a property of intact tissue that the isolated myocyte does not display.

# Results

## The resting muscle transcriptome carries no information about insulin sensitivity that is independent of age and adiposity

We used the largest public dataset with both a resting transcriptome and a gold-standard phenotype: 49 individuals profiled in two laboratories, all with insulin sensitivity measured by hyperinsulinaemic–euglycaemic clamp (Fig. 3a).

Correlating basal expression with the M-value without covariate adjustment returns 739 genes at a false discovery rate below 0.1, with 85 of 211 mitochondrial genes among them, reproducing the mitochondrial enrichment previously reported for this cohort. The association is not independent of the covariates. In these participants the M-value is strongly correlated with body-mass index (ρ = −0.58, P = 1 × 10⁻⁵) and with age (ρ = −0.39, P = 0.005), and adjusting expression for either one alone removes every gene at that threshold (Fig. 3b, Extended Data Fig. 1). Body-mass index by itself has no genes at FDR < 0.1, so the signal arises from variance shared between adiposity, age and insulin sensitivity rather than from a separate adiposity signature.

In the fully adjusted analysis, with the M-value residualised on the same covariates as expression, no gene among 21,595 reaches FDR < 0.1 and 231 reach P < 0.01 where 215 are expected by chance (Fig. 3b). Diabetic and normoglycaemic individuals do not differ at rest either, here or in an independent muscle cohort (Fig. 2g). Of ten curated gene sets only the oxidative programme separates from size-matched random sets (Fig. 3d). Because age and adiposity are themselves determinants of insulin sensitivity, a cross-sectional resting design cannot distinguish a transcriptional correlate of sensitivity from a correlate of its determinants. This is a limitation of the design rather than of the cohort, and it motivates the move from state to response.

## Healthy muscle answers insulin along a shared direction, established over hours

We represented each biopsy by its position in the principal subspace of the 800 most variable genes and each person's response as the displacement between their paired biopsies (Fig. 4a). Two properties of a group's responses are then independent of magnitude: the direction each person moves in, and how concentrated those directions are. We quantify concentration per person, as the alignment of each participant's displacement with the leave-one-out mean direction of the reference group, and at the group level by fitting a von Mises–Fisher distribution over the sphere with mean direction μ and concentration κ (Methods).

In 20 insulin-sensitive individuals, coherence 4 h into a clamp was 0.77 and κ was 8.2 (Fig. 4b): healthy people answer insulin along one direction. Coherence was low in cohorts sampled earlier, 0.20 at 30 min, 0.24 at 1 h after a mixed meal and 0.19 at 2 h, and high at 4 h in two independent cohorts (Fig. 4c), and the same ordering appears gene by gene within the single cohort sampled at two times (Fig. 4d). A shared direction is therefore not an immediate consequence of the stimulus; it is assembled over hours.

The genes carrying that direction, replicated across two healthy cohorts, extend well beyond canonical insulin signalling: immediate-early and stress transcription factors, metallothioneins and HMOX1, chaperones, the amino-acid transporter SLC7A5, and BHLHE40 (Fig. 4e). Insulin repressed the clock-output genes DBP and NR1D2 and induced PER2, and did so only at 4 h (Fig. 4f). Because a second biopsy taken near the first induces immediate-early genes irrespective of the stimulus, we verified that this cannot explain the clock result: immediate-early genes respond at every sampling time, clock-output genes only at 4 h (Fig. 4g).

## Insulin resistance preserves magnitude but reduces the concentration of response directions

In the same cohort, insulin-resistant and diabetic individuals moved as far as insulin-sensitive ones; a two-one-sided test does not establish equivalence at this sample size, so we report no detected reduction in magnitude rather than preserved magnitude (Fig. 5a).

What differed was direction. Per-person alignment with the healthy response direction fell from 0.75 to 0.31 in insulin-resistant and 0.26 in diabetic muscle (differences of 0.44, P = 0.0017, and 0.49, P = 0.0010; permutation with full re-estimation of the reference direction). Hybridisation batch is unevenly distributed across groups in this cohort and can by itself generate apparent coordination, so we take as primary the analysis restricted to the 35 individuals whose two biopsies were processed in the same batch: there the differences are 0.51 (P = 0.013) and 0.47 (P = 0.006), and groups defined by hybridisation run alone show no difference (−0.10, P = 0.57) (Fig. 5b, Extended Data Fig. 3). The von Mises–Fisher fit puts a parameter on this: κ falls from 8.2 to 3.6 and 4.0, a 2.3-fold loss of concentration (likelihood ratio P = 0.0016 and 0.0086; Fig. 5c).

Concentration is not a restatement of magnitude in these data. Across 29 groups the correlation between the two is 0.33, and within this cohort it is inverted: the insulin-resistant group has the largest mean displacement and nearly the lowest alignment. This differs from single-cell CRISPR screens, where directional coherence tracks effect magnitude with correlations of 0.84 to 0.98, and the difference is expected, since variation between people has sources that isogenic cells do not have.

The canonical insulin programme was intact: TXNIP, KLF15 and IRS2 were repressed and PPP1R3B and HES1 induced in every participant of an independent prediabetic cohort. Of the 55-gene healthy programme, a third was lost in insulin-resistant and diabetic muscle and essentially none was inverted (Fig. 5i). Among nine pre-specified clock-output genes, the group-by-insulin interaction was significant for DBP, PER2 and NR1D2 in the insulin-resistant comparison and borderline for five genes in the diabetic comparison, whereas BHLHE40 responded in all three groups (Fig. 5f); the same contrast appeared in primary myotubes exposed to high glucose plus insulin, where clock-output oscillations were lower in amplitude in diabetic cells (Fig. 5g,h). No single gene survives genome-wide correction for the interaction, with 95 genes at P < 0.001 where 8 are expected: the signature of a distributed effect, which is what the directional analysis is built to capture.

## The reduction is specific to insulin

Diabetic muscle might simply have lost the ability to respond in concert to anything. We tested this with a different physiological perturbation in the same tissue, in participants recruited by the same group: an acute bout of exercise with biopsies at rest and three hours into recovery, in 20 men with type 2 diabetes and 17 with normal glucose tolerance, and the same protocol sampled in adipose tissue.

Diabetic muscle responds to exercise as coherently as healthy muscle (0.53 versus 0.45; Fig. 5j), and the adipose arm agrees. To estimate the contrast jointly we fitted a hierarchical von Mises–Fisher model to all 140 participants from four cohorts, each cohort keeping its own mean direction, with log κ = b₀ + b_cohort + b₁·impairment + b₂·(impairment × exercise). Metabolic impairment reduced κ by a factor of 0.56 under insulin (b₁ = −0.585, 95% CI −1.31 to +0.01) and the interaction with exercise was positive and significant (b₂ = +0.864, 95% CI +0.13 to +2.15; likelihood ratio **P = 0.005**), abolishing the effect under exercise.

Two observations follow. Immediately after exercise the response is idiosyncratic in both groups and becomes concentrated only at three hours of recovery, reproducing in a different stimulus the time dependence seen with insulin. And because the exercise response is intact, the reduced concentration of the insulin response is not a general transcriptional failure of diabetic muscle, nor an artefact of biopsy trauma or of the clinical state of the participants.

High concentration is, moreover, specific to acute stimuli. Values near 0.8 occur for insulin at 4 h and for acute exercise in muscle, whereas slow perturbations produce individual rather than shared trajectories: weeks of training give 0.31 in muscle and −0.09 in adipose tissue in the same participants, a twelve-hour diurnal contrast gives 0.10, and ten days of cold acclimation give −0.01 (Supplementary Note 3).

## The coordinated response is a property of the tissue, not of the myocyte

Where in the tissue does this coordination reside? Muscle biopsies contain myofibres, vasculature, immune cells and progenitors in proportions that differ between people, and insulin acts on several of them.

Fibre type does not account for the result. Using fibre-type signatures validated against myosin heavy-chain determination, the slow-to-fast axis does not differ between groups (P = 0.83), does not predict individual alignment (ρ = −0.07), and adjusting for it leaves the effect intact (0.46 to 0.40, P = 0.0013). The mononuclear compartment does. Estimating the proportions of fibro-adipogenic progenitors, endothelium, pericytes, satellite cells, macrophages and lymphocytes by deconvolution against a single-cell reference of human vastus lateralis, and adjusting each person's displacement for those proportions, reduces the group difference from 0.46 to 0.08 (P = 0.50), whereas adjusting for the same number of random covariates leaves it unchanged (0.44, P = 0.0007) (Fig. 6).

If coordination depends on the composition of the tissue, it should be absent from the isolated myocyte. We tested this directly in primary myotubes from 24 donors, 12 with normal glucose tolerance and 12 with type 2 diabetes, cultured under identical conditions and stimulated with 100 nM insulin, with samples at 0, 0.5, 1 and 2 h. Coordination was weak in everyone: the highest coherence at any time point was 0.27 in healthy donors, against 0.77 in intact healthy muscle, and the mean displacement was 2.2 against 16. It did not differ between groups at any time point (differences of −0.27, +0.09 and +0.21; P = 0.32, 0.69 and 0.41), and κ at 2 h was 1.90 and 1.96 (Fig. 6).

The myocyte therefore retains the donor phenotype, as these cultures are known to do and as the resting expression of these same cells shows, but it does not produce a coordinated response to insulin, and consequently it cannot be the site at which coordination is lost.

# Discussion

Insulin elicits from human skeletal muscle a transcriptional response that is shared between individuals, is assembled over hours, extends well beyond canonical insulin signalling, and includes modulation of clock output. In insulin resistance that response retains its size and its canonical programme but loses its common direction. The loss is specific to insulin, since the same muscle answers acute exercise in concert, and it is a property of the organised tissue, since it is accounted for by cellular composition and is absent from myocytes cultured from the same people.

This reframes what insulin resistance looks like at the transcriptional level. The question has generally been posed as which pathway is defective in the myocyte. Our results suggest the unit of failure is larger than the cell: a coordinated response requires the assembled tissue, where insulin acts on endothelium, pericytes, immune cells and progenitors as well as on the myofibre, and where its delivery to the interstitium is itself regulated. A cultured myocyte has the donor's genome and epigenome but not that architecture, and it does not produce the response at all. That the healthy response takes hours rather than minutes fits the same picture: a direction that must be assembled is more plausibly a tissue-level process than a cell-autonomous one.

It also explains an old frustration. Signatures of insulin resistance in resting muscle have rarely replicated beyond the oxidative programme, and in the one cohort with a gold-standard phenotype we find that the reported association is not separable from age and adiposity. If the disease is visible in a coordinated tissue response rather than in a resting state, resting profiles of any size will keep finding little, and the informative design is two biopsies around a perturbation.

Two methodological points follow for anyone measuring response geometry in human cohorts. A directional coherence statistic has been introduced for single-cell screens, where it proved to be largely a restatement of effect magnitude; between people the two are largely dissociated, and in our key cohort they are inverted. And a group-level statistic discards the individual observations that make the comparison powerful: using per-person alignment rather than group coherence changed P from 0.10 to 0.013 on identical data.

The immediate question our results raise is where within the tissue the coordination is generated and which compartment degrades it, and this requires an experiment we could not assemble from public data: single-nucleus profiling of paired biopsies around a clamp. That design separates three possibilities that bulk tissue cannot. If myonuclei diverge between people, the defect is in the muscle cell after all, in a form the cultured myocyte does not retain. If myonuclei within each person are concentrated and the divergence lies in the vascular and immune compartments, coordination is set by the tissue environment, as our results suggest. And if myonuclei within a single person already diverge, the heterogeneity is intracellular and a different phenomenon from the one measured here. The measure used in this work requires only what such a study would already collect, so the question is answerable with existing technology.

## Limitations of the study

The insulin contrast comes from a single clamp cohort in which hybridisation batch is unevenly distributed; the primary analysis is restricted to same-batch pairs, where the effect is undiminished and the batch-only control is null, but replication in a modern sequenced clamp cohort remains the most valuable next step and no such deposited dataset currently exists. The hierarchical model estimates the interaction well and the main effect of impairment only imprecisely. The myotube experiment reaches 2 h, whereas coordination in tissue is established nearer 4 h, so it bounds rather than excludes a late cell-autonomous component. Composition and coordination cannot be fully separated in bulk data: adjusting for composition may remove part of the phenomenon if altered composition is itself part of the pathophysiology. The time course in Fig. 4c combines independent cohorts. All cohorts are cross-sectional at the group level and no claim here is causal.

# Methods

Each subsection states the question the procedure was designed to answer before describing how it was implemented. Code implementing every step, with fixed seeds, is available as described under Code availability.

## Cohorts and preprocessing

Twenty-four independent studies across 29 accessions and 1,689 participants were used, after removing donors shared between accessions of the same study; details, platforms and group sizes are in Supplementary Table 1. Cohorts were included if they provided human transcriptomes with a glycaemic or insulin-sensitivity phenotype, or paired sampling around a perturbation. Processed expression was used as deposited where available; counts were converted to log2 counts per million with an expression filter; probes were mapped to symbols with the platform annotation. Where a technical factor was unevenly distributed across groups it was either regressed out or used to define a restricted analysis, as stated per analysis.

## Resting associations

Differential expression at rest used a two-sided Welch *t*-test per gene on expression residualised for the covariates of each study, with Benjamini–Hochberg correction; the contrast of each cohort was declared explicitly rather than inferred. Association with insulin sensitivity used partial Spearman correlation, with the M-value residualised on the same covariates as expression. The ladder from unadjusted to fully adjusted models is reported because it is the substance of the result.

## Geometry of the within-person response

Samples were projected into the principal subspace of the 800 most variable genes (five components) and each person's response taken as the displacement between their paired biopsies. *Alignment* is the cosine between a participant's displacement and the reference group's mean direction computed without that participant, and is the primary outcome: group coherence collapses n individuals into one number and so compares two observations rather than n. *Coherence* is the mean of those cosines within a group and is reported for comparability with other work. *Magnitude* is the Euclidean norm.

Group comparisons of alignment use permutation of individual labels with full re-estimation of the reference direction in each permutation, since that direction is defined by the reference group; 3,000 permutations. Intervention comparisons permute the two time points within each participant. Per-gene group-by-stimulus interactions use permutation with Benjamini–Hochberg correction over the nine pre-specified clock-output genes and are reported as a targeted test; the genome-wide interaction test is reported separately.

Two properties of the statistics required explicit control. The naive coherence, which includes each individual in the mean it is compared with, returns 0.41 at n = 4 and 0.19 at n = 20 when no direction is shared; the leave-one-out form returns zero at every group size and is used throughout. Nulls that destroy a shared direction were generated by replacing displacements with isotropic directions of the same length, since a sign flip preserves the original axes and gives a lenient null. Because alignment is computed from the expression matrix, any signature derived from the same matrix correlates with it by construction; such correlations were therefore evaluated against size-matched random gene sets rather than against zero.

## Von Mises–Fisher and hierarchical models

Each person's response direction was taken as a unit vector and modelled as a draw from a von Mises–Fisher distribution with mean direction μ and concentration κ, fitted by maximum likelihood with the normalising constant evaluated through the exponentially scaled modified Bessel function. Groups were compared by likelihood ratio against the restricted model in which both share κ but keep separate directions; confidence intervals for κ come from 400 bootstrap resamples.

The hierarchical model pools all paired cohorts, each keeping its own mean direction because tissue, stimulus and platform differ, with log κ = b₀ + b_cohort + b₁·impairment + b₂·(impairment × exercise). Significance is by likelihood ratio against the nested models; confidence intervals come from bootstrap resampling of participants stratified by cohort.

## Cellular composition

Fibre-type proportions were estimated from the twenty markers per type validated against myosin heavy-chain determination in single human fibres, summarised as the contrast between the slow and fast marker sets. Mononuclear proportions were estimated by non-negative least squares against a signature matrix built from single-cell RNA-seq of human vastus lateralis mononuclear cells (2,876 cells, four samples), clustered and annotated by canonical markers into fibro-adipogenic progenitors, endothelium, pericytes and smooth muscle, satellite cells, macrophages and lymphocytes. Composition was adjusted for by regressing centred, scaled covariates out of each person's displacement, which preserves the group mean direction; proportions were transformed to centred log-ratios because they are compositional. The same procedure applied to matched numbers of random covariates provides the control.

## Myotube experiment

Primary myotubes from 24 donors (12 normal glucose tolerant, 12 with type 2 diabetes, balanced for obesity and sex) were stimulated with 100 nM insulin with sampling at 0, 0.5, 1 and 2 h. Displacements were computed per donor between 0 h and each later time point and analysed exactly as the tissue data.

## Reporting

Exact sample sizes, statistics and *P* values are in the figure legends. Random seeds are fixed. Analyses used Python 3.10 and R 4.3; versions are pinned in the repository.

# Data availability

All data analysed are public. Accessions are listed in Supplementary Table 1. No new data were generated.

# Code availability

All code is at https://github.com/Danpc11/t2d_landscape.
