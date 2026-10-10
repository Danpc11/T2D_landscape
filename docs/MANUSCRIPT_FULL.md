---
title: "Insulin resistance disrupts transcriptional coordination in skeletal muscle"
geometry: margin=2.4cm
fontsize: 10pt
linestretch: 1.25
colorlinks: true
---

**Authors.** [to be completed]

**Affiliations.** [to be completed]

**Correspondence.** [to be completed]

**Article type.** Analysis (*Nature Metabolism*)

**Word count.** Main text ~4,600; six main figures; Supplementary Information with one supplementary figure and five supplementary tables.

\newpage

# Abstract

Type 2 diabetes is defined by what insulin fails to do, yet most molecular studies profile tissues at rest. Whether the transcriptional signature of insulin resistance lies in the resting state of a tissue or in the way it responds has not been resolved. We reanalysed 23 public human transcriptomic cohorts spanning pancreatic islet, skeletal muscle, adipose tissue and blood, including 253 GTEx donors with three organs sampled from the same individual, 19 surgical patients with paired adipose depots, 49 individuals with resting muscle biopsies and clamp-measured insulin sensitivity, and 123 individuals biopsied before and during hyperinsulinaemic clamps. Across organs, coexpression architecture is shared far beyond chance and each person occupies a consistent position within it, whereas resting geometry separates glycaemic stages poorly and no cohort shows two stable states. In resting muscle, no gene is associated with insulin sensitivity at a false discovery rate below 0.1. The response to insulin is informative: healthy muscle moves along one direction shared across individuals, established over hours, that includes late modulation of clock-output genes, whereas insulin-resistant and diabetic muscle respond with comparable magnitude but without a shared direction and with altered clock-output responses. Adipose tissue instead loses response magnitude while preserving direction. Insulin resistance therefore disrupts the coordination of the muscle transcriptional response rather than its size, a property invisible at rest and measurable in any study that biopsies before and during a clamp.

\newpage

# Introduction

Skeletal muscle disposes of most postprandial glucose and its resistance to insulin is among the earliest defects on the path to type 2 diabetes^1,2^. The molecular description of that defect has been sought mainly in two places. In insulin signalling, where impaired Akt phosphorylation is reproducible but accompanied by heterogeneous downstream effects^3,4^. And in the resting transcriptome, where two decades of profiling have produced a short list of consistent findings, essentially the reduced oxidative programme marked by PGC-1α^5,6^, alongside many signatures that have not replicated across cohorts.

A third place has been explored less. Insulin regulates several hundred muscle genes within hours^7,8^, and specific insulin-responsive genes are blunted in obesity and diabetes^9^. Insulin is also a timing signal: it resets the clock of adipose tissue and of isolated adipocytes^10^ and of peripheral tissues in rodents^11^. The muscle clock in turn governs insulin sensitivity^12^, circadian misalignment lowers insulin sensitivity in humans^13^, and muscle cells from people with type 2 diabetes show intrinsically disrupted circadian oscillations^14^. These literatures have not been joined: whether insulin's modulation of clock output is part of the normal response of human muscle in vivo, and whether it is altered in insulin resistance, has not been examined. Recent proteomic and phosphoproteomic work has shown that the molecular signatures of insulin resistance are highly individual^15,16^, raising the possibility that what differs between insulin-resistant people is not a shared defect but the loss of a shared response.

Two obstacles have kept this question open. The first is measurement. Comparing responses between groups requires a statistic that separates *how far* a tissue moves from *where* it moves to, because a group whose members all respond strongly but in different directions produces the same average as a group that barely responds at all; conventional differential expression conflates the two. The second is design. The informative experiment is a perturbation with paired biopsies in the same person, which exists in public repositories but in small, scattered cohorts that have not been analysed together.

Here we ask whether the information about insulin sensitivity in human muscle lies in how the tissue *is* or in how it *responds*. We first establish, across organs and cohorts, what the resting state does and does not show (Figs. 1–3), then define a measure of the geometry of a within-person response and apply it to every public clamp cohort we could obtain (Figs. 4–6).

# Results

## The metabolic organs share one transcriptional architecture and each person occupies one position in it

We first needed to know whether the organs of metabolism can be treated as one system at all, because if each organ drifts independently then an individual's muscle tells us nothing about that individual beyond muscle. To ask this we needed a comparison that does not depend on which genes are expressed in which tissue, since a hepatocyte and a myofibre share few transcripts at comparable levels but may still be organised alike. We therefore compared organs through the *organisation* of their coexpression networks rather than through their expression levels: each organ contributes the leading spectral subspace of its network, the organs are aligned by orthogonal maps, and the residual distance after alignment measures how differently they are organised (Fig. 1a,d; Methods). A null that permutes the gene labels of each tissue independently destroys the correspondence between organs while preserving every within-organ property, which isolates exactly the quantity of interest.

Within each organ, standard network descriptors behave along the glycaemic axis as expected at equal sample size: global efficiency falls and effective resistance rises from healthy to diabetic in islet and adipose tissue, while muscle is flat (Fig. 1b,c). These descriptors treat each organ separately and say nothing about whether the organs are organised alike.

They are. In three discovery cohorts of islet, muscle and adipose tissue from different individuals, the residual distance was 7 to 11 standard deviations below the correspondence null (Fig. 1e). Because those cohorts contain different people, we repeated the analysis in GTEx v11 donors with skeletal muscle, subcutaneous adipose tissue and pancreas sampled from the same individual: after adjusting for age, sex, death classification, RNA integrity, ischaemic time and batch, the distance was 26 standard deviations below the null (Fig. 1f).

Paired tissues allow a stronger question: is the shared architecture also a property of the person? The first canonical correlation between a donor's scores in muscle and in adipose tissue was 0.93, and 0.91 and 0.89 for the two other pairs, against 0.23 when donors were shuffled between organs (Fig. 1g). The same held in 19 living bariatric-surgery patients with paired subcutaneous and omental biopsies (0.93 versus 0.59), indicating that the effect is not an artefact of post-mortem sampling. It survived the addition of technical covariates (Fig. 1h). Whole blood from the same GTEx donors predicted 23 to 30% of the principal variation of muscle, adipose tissue and pancreas in cross-validation (Fig. 1i), showing that part of this shared state is visible in an accessible tissue. The genes that carry the shared position form a stress, glucocorticoid and inflammatory module together with ribosome-biogenesis genes and PTPN1 (Fig. 1l).

## Resting geometry does not separate glycaemic stages, and no cohort shows two stable states

If health and diabetes were two attractors of this shared state, the distribution of individuals should be bimodal with a saddle between the basins, and the intermediate stage should show the signatures of a system approaching a transition. This matters practically: a bistable disease would justify searching for a tipping point and for interventions timed to it, whereas a continuous drift would not. We estimated the quasi-potential on a covariate-adjusted, stage-balanced embedding and counted basins by sublevel-set persistence, requiring stability across kernel bandwidths rather than accepting the modes of a single smoothing choice. In eight cohorts, three tissues, two platforms and 73 to 434 individuals each, there was one basin (Fig. 1j). Where a mixture model favoured two components, the components were unrelated to glycaemic stage (Fig. 1k).

Nor is the intermediate stage distinguished. The critical-transition index, which would peak at an intermediate state in a system approaching a transition, shows no such peak in any organ (Fig. 2c). The dispersion of individuals in diabetes relative to health crosses one in every cohort (Fig. 2d). Between organs, coherence weakens monotonically with stage but no stage-to-stage difference reaches significance (Fig. 2b). The resting organisation therefore drifts with disease rather than switching, and the drift is too small, at these sample sizes, to be resolved stage by stage.

## In resting muscle, no gene is associated with insulin sensitivity

We then asked how much the resting transcriptome of a single organ knows about the physiological phenotype, using the largest public dataset with both: 49 individuals profiled in two laboratories, all with insulin sensitivity measured by hyperinsulinaemic–euglycaemic clamp (Fig. 3a). The point of using clamp data rather than diagnosis is that the diagnosis is a threshold on a continuum, whereas the M-value is the quantity the tissue is supposed to determine.

After adjusting for laboratory, age and body-mass index, no gene among 21,595 was associated with the M-value at a false discovery rate below 0.1, and 14 genes reached *P* < 0.01 where approximately 215 are expected by chance (Fig. 3b). The comparison of diabetic with normoglycaemic individuals at rest gave no gene below that threshold either (Fig. 3c). Positive controls behaved as expected: PPARGC1A tracked the M-value nominally and PDK4 was higher in diabetes, reproducing the oxidative signature that has been the field's most consistent resting finding^5,6^. At the level of gene sets, only the oxidative programme separated from size-matched random sets (Fig. 3d). The clock-output genes that become central below show no association with insulin sensitivity at rest (Fig. 3e). Aggregate geometry is similarly uninformative, although the distance of each individual from the healthy reference correlates weakly with the M-value, a diffuse signal that no individual gene captures.

This is a negative result with the design and sample size to be informative: two laboratories, a gold-standard phenotype in every participant, and the known markers recovered. For comparison, the same analysis applied to islet and adipose tissue at rest returns hundreds to thousands of genes and the expected biology, inflammation and complement, lipid handling and extracellular matrix, while muscle returns two (Fig. 2g,h). The contrast is not that the canonical approach fails in general; it is that in muscle the resting state carries little, which motivates the shift from state to response.

## The healthy response to insulin is coordinated across individuals and is established over hours

To compare responses between groups we needed a description of a response that does not depend on choosing genes in advance and that separates size from direction. We represented each biopsy by its position in the principal subspace of the most variable genes and each individual's response as the displacement between the paired biopsies. Two properties of a group's responses are then independent of magnitude: the *coherence*, the mean cosine between each individual's displacement and the leave-one-out group mean, which is one when everyone moves alike and zero when directions are random; and the group *direction* (Fig. 4a). The leave-one-out construction matters because including an individual in the mean it is compared with inflates agreement in small groups.

In 20 insulin-sensitive individuals, coherence after 4 h of insulin was 0.77: healthy people answer insulin along a shared direction (Fig. 4b). Coherence was lower in cohorts sampled earlier, 0.20 at 30 min, 0.24 at 1 h after a mixed meal and 0.19 at 2 h, and higher at 4 h in two independent cohorts (Fig. 4c). These are separate cohorts rather than a within-person time course, but the pattern is consistent with a response whose shared direction is established over hours rather than minutes; the same ordering appears gene by gene within the single cohort sampled at two times (Fig. 4d).

The genes that carry the shared 4-h direction, replicated across two healthy cohorts, are not the canonical insulin-signalling genes alone but a broader programme: immediate-early and stress transcription factors, metallothioneins and HMOX1, chaperones, the amino-acid transporter SLC7A5, and BHLHE40, the transcription factor through which the circadian clock couples to metabolism (Fig. 4e). In the same individuals, insulin repressed the clock-output genes DBP and NR1D2 and induced PER2, and did so only at 4 h (Fig. 4f). Because a second biopsy taken near the first induces immediate-early genes irrespective of the stimulus, we verified that this artefact cannot explain the clock finding: immediate-early genes respond at every sampling time, whereas clock-output genes respond only at 4 h (Fig. 4g). All groups underwent the same two-biopsy protocol, so comparisons between groups are unaffected.

## Insulin resistance preserves magnitude but reduces directional coherence and alters the clock-output response

In the same cohort, 20 insulin-resistant and 15 diabetic individuals responded to insulin with magnitudes indistinguishable from the sensitive group (Fig. 5a). What differed was the agreement between individuals: coherence fell from 0.77 to 0.35 in insulin-resistant and 0.45 in diabetic muscle, beyond the distribution obtained by permuting group labels (Fig. 5b). Because hybridisation batch was unevenly distributed across groups in this 2007 dataset^7^, we repeated the analysis in the 35 individuals whose two biopsies were processed in the same batch: the effect sizes were unchanged, and at this sample size the coherence comparisons no longer reach significance while the direction contrast does (Fig. 5c). Individually, insulin-sensitive people point in a narrow cone around the group direction, whereas insulin-resistant and diabetic people point in all directions (Fig. 5d). Low coherence was reproduced at baseline in two independent arms of a prediabetes trial profiled by RNA sequencing in a different laboratory^17^ (Fig. 5e).

The canonical insulin programme was intact in these individuals: TXNIP, KLF15 and IRS2 were repressed and PPP1R3B and HES1 induced in every participant of the prediabetic cohort. What differed was specific. Of the 55-gene healthy programme, a third was lost in insulin-resistant and diabetic muscle and essentially none was inverted (Fig. 5i). Among clock-output genes, the group-by-insulin interaction was significant for DBP, PER2 and NR1D2 in the insulin-resistant comparison and borderline for five genes in the diabetic comparison, whereas BHLHE40 responded in all three groups (Fig. 5f). The same contrast appeared in primary myotubes from donors with and without diabetes exposed to chronic insulin^14^, where DBP was repressed in cells from normoglycaemic donors more than in diabetic cells while BHLHE40 responded in both, and the amplitude of clock-output oscillations was approximately 30% lower in diabetic cells (Fig. 5g,h). Resting expression of these same genes carries no information about insulin sensitivity (Fig. 3e), so this is a property of the response and not of the state.

## Adipose tissue fails differently, and interventions restore magnitude without demonstrably restoring coherence

Adipose tissue sampled before and during a clamp showed the opposite pattern. Obese and insulin-resistant adipose tissue responded along the healthy direction, with cosines of 0.88 and 0.95, at approximately half the magnitude (Fig. 6a–c). Two years after bariatric surgery, in the same 23 women^18^, the magnitude of the response had returned to the level of non-obese controls while coherence had not (Fig. 6a). Ten weeks of nicotinamide mononucleotide, which improved clamp-measured insulin sensitivity in prediabetic women^17^, did not change coherence either (Fig. 6d). At the available sample sizes a restoration smaller than approximately 0.25 would not be detected (Fig. 6e), so these are absences of evidence rather than evidence of absence; what they establish is that coordination is not restored in step with the measures that these interventions were designed to improve.

# Discussion

Three statements summarise this work. The resting transcriptome of human muscle is, at the resolution of these cohorts, uninformative about insulin sensitivity. The healthy response to insulin is coordinated across individuals, is established over hours, and includes late modulation of clock-output genes. Insulin resistance preserves the magnitude of that response and the canonical signalling programme but reduces its directional coherence and alters the clock-output response.

The first statement helps explain an old frustration. Signatures of insulin resistance in resting muscle have rarely replicated beyond the oxidative programme; in 49 clamped individuals from two laboratories we find none at a false discovery rate below 0.1, with the oxidative markers present as nominal controls. If the disease lives in the response, resting profiles of any size will keep finding little, and the informative design is two biopsies around a perturbation. A compatible conclusion has been reached from phosphoproteomics, where the signatures of insulin resistance are highly individual^15,16^. Our result adds a geometry to that observation: individuals respond with comparable magnitude but not in agreement with one another.

The second and third statements connect literatures that have developed separately. Insulin is a timing signal for peripheral clocks^10,11^; the muscle clock governs insulin sensitivity^12^; misalignment lowers insulin sensitivity in humans^13^; and the clock of diabetic muscle cells is intrinsically disrupted^14^. We observe the corresponding step in vivo: modulation of clock output is part of the normal 4-h response of human muscle, it is altered in insulin resistance, and BHLHE40, the clock's metabolic effector, keeps responding. DBP, TEF and HLF are direct BMAL1 targets and the PAR-bZip output through which the clock times metabolism^19^; their altered response, with resting expression unchanged, points to the coupling rather than to the oscillator. A response that takes hours to establish and involves clock output is the kind of response expected to lose agreement between individuals when that coupling is disturbed, which is what we observe. We stop short of causal language: these are associations in cross-sectional cohorts with paired sampling, and the direction of causation between coupling and coherence is not established here.

The adipose result shows that organs fail differently, muscle by direction and adipose tissue by magnitude, which argues against treating insulin resistance as one molecular phenotype across tissues. That bariatric surgery restored adipose magnitude without restoring coherence, and that nicotinamide mononucleotide improved sensitivity without changing it, suggests that coordination is a property distinct from the measures these interventions target; the sample sizes allow this to be stated as a hypothesis and not more.

A recent synthesis of skeletal-muscle biology argues that muscle health is defined not by the maximal activation of any pathway but by their coordination over time, and lists the mechanisms that explain inter-individual variability in insulin sensitivity among phenotypically similar people as an open question^20^. The measure used here addresses that question directly and requires only what a clamp study already collects: two biopsies and a group of participants. If it proves robust in larger cohorts, the coherence of the response is a candidate endpoint for interventions aimed at the timing of metabolic signals, such as meal and insulin scheduling or sleep, which current endpoints do not capture.

## Limitations of the study

The central comparison rests on one clamp cohort from 2007 in which hybridisation batch was unevenly distributed across groups; restricting to same-batch pairs preserves the effect sizes but not the significance of the coherence contrasts at approximately 11 individuals per group, and the diabetic direction in particular needs replication. Two sequential biopsies induce an immediate-early response independent of the stimulus; this does not affect comparisons between groups but inflates the apparent insulin specificity of the healthy gene list. The time course in Fig. 4c combines independent cohorts rather than repeated sampling of the same individuals. The prediabetic replication is RNA sequencing from a single laboratory; the adipose cohorts are small; the myotube comparison has seven and five donors. GTEx donors are post mortem, so part of the shared individual state may reflect terminal illness, although the in vivo replicate in 19 surgical patients argues that the phenomenon is not only agonal. Restoration analyses are underpowered to detect partial recovery. Pancreatic islet could not be included in the response analyses because no in vivo perturbation with paired sampling exists for that organ in humans; an ex vivo equivalent is proposed in the Supplementary Information. All cohorts are cross-sectional at the group level, and no claim here is causal.

# Methods

Methods are written purpose first: each subsection states what question the procedure was designed to answer before describing how it was implemented. Code implementing every step, with fixed random seeds, is available as described under Code availability; each numbered script corresponds to one subsection below.

## Cohorts and preprocessing

*Purpose.* To assemble every public human dataset that could address either the resting state or the response, without selecting cohorts on their results.

Twenty-three cohorts were used: islet (GSE76895, GSE164416, GSE50244, GSE50398), skeletal muscle (GSE18732, GSE25462, GSE182120, GSE22309, GSE9105, GSE7146, GSE231509, GSE157988), adipose tissue (GSE27951, GSE64567, GSE135134, GSE20950, GSE26637, GSE59034 and the CAGE data of Rydén et al.), myotubes (GSE182117), whole blood and monocytes (GSE66306, GSE156993, GSE21321), and GTEx v11 (skeletal muscle, subcutaneous and visceral adipose tissue, pancreas, liver and whole blood, with the open sample and subject annotations). Cohorts were included if they provided human transcriptomes with either a glycaemic or insulin-sensitivity phenotype, or paired sampling around a perturbation; none was excluded after analysis.

Processed expression was used as deposited where the authors provided it. Sequencing counts were converted to log2 counts per million with an expression filter requiring the gene to exceed one count per million in at least half the samples. Probes were mapped to gene symbols with the platform annotation and averaged within symbol. Where a dataset contained technical factors that were unevenly distributed across groups, these were either regressed out or used to define a restricted analysis, as stated for each analysis.

## Coexpression networks and their descriptors

*Purpose.* To describe how each organ is organised, rather than which genes it expresses, so that organs with different transcriptomes can be compared on the same footing; and to test whether that organisation changes along the glycaemic axis within an organ.

For each organ and stage, a signed weighted network was built on the genes shared across the organs being compared, using biweight midcorrelation raised to a soft-thresholding power chosen once in the healthy state by the scale-free topology criterion and held fixed for the other stages so that differences between stages cannot arise from the thresholding. Because the number of individuals differs between stages and every network descriptor depends on sample size, all descriptors were computed on subsamples of equal size drawn without replacement, repeated 20 times, and reported as mean and 95% confidence interval across repetitions.

Three descriptors were used, chosen because each answers a different question. Global efficiency, the mean inverse shortest weighted path, measures how directly information can travel across the network. Mean effective resistance, from the pseudoinverse of the Laplacian, measures the same traffic as a diffusion problem and is sensitive to the loss of redundant routes rather than of individual edges. Spectral entropy of the normalised Laplacian measures how evenly the organisation is distributed across modes. Modules were detected by blockwise hierarchical clustering of the topological overlap matrix. Values are reported both absolutely and relative to a degree-preserving configuration null, so that changes attributable only to the degree distribution are excluded.

## The cellular sheaf: comparing the organisation of different organs

*Purpose.* To ask whether organs that express different genes are nevertheless organised alike, and whether an individual occupies a consistent position across their organs. The obstacle is that expression levels are not comparable between organs; the solution adopted here is to compare the geometry of each organ's network rather than its expression.

Each organ contributes a *stalk*: the subspace spanned by the leading eigenvectors of its normalised Laplacian, which summarises the organisation of its network independently of overall expression. Organs are compared after an orthogonal alignment, because a rotation of the spectral basis is arbitrary and should not count as a difference. Alignments are fitted by Procrustes to a neutral reference so that no organ is privileged. The energy is the mean squared Grassmann distance between aligned stalks, so that it is zero when the organs share one architecture and grows as they diverge; it does not depend on the sign or scale of the eigenvectors.

Two nulls were used because they answer different questions. The gene-correspondence null permutes the gene labels of each organ independently, preserving every within-organ property and destroying only the correspondence between organs; it tests whether the shared architecture is more than the coincidence of two well-structured networks. The donor-shuffle null permutes individuals between organs, preserving each organ's structure and each individual's profile within an organ; it tests whether an individual's position is their own. In the paired-tissue analyses, expression was residualised for age, sex, death classification, RNA integrity number, ischaemic time and sequencing batch before networks were built, and the analysis repeated with and without the technical terms to show what they change.

Individual position was summarised by canonical correlation between the principal components of two organs computed on the same donors, with the donor-shuffle null as reference. Blood-to-tissue prediction used ten blood components to predict five tissue components under five-fold cross-validation, reporting R² on held-out donors only.

## Quasi-potential landscape

*Purpose.* To test the common assumption that health and diabetes are two stable states separated by a barrier, which if true would justify searching for a tipping point.

Individuals were embedded in the principal subspace of their covariate-adjusted, stage-balanced expression. The density was estimated with a Gaussian kernel and converted to a quasi-potential by taking its negative logarithm. Basins were counted by sublevel-set persistence rather than by counting local minima, because a minimum created by sampling noise and a genuine basin are indistinguishable without a notion of how deep a basin must be to survive resampling. The persistence threshold was set to twice the standard deviation of the persistence of the leading basin across 200 bootstrap resamples, and basins were required to survive at three kernel bandwidths so that no conclusion depends on one smoothing choice. As an independent check, Gaussian mixture models were fitted and the association between the assigned component and glycaemic stage was measured by Cramér's V; a mixture that is favoured by the Bayesian information criterion but unrelated to stage is evidence against bistability of the disease even when the data are formally bimodal. The critical-transition index was computed as the mean absolute gene–gene correlation divided by the mean absolute donor–donor correlation at equal sample size.

## Geometry of the within-person response

*Purpose.* To compare responses between groups in a way that distinguishes a group that responds weakly from a group that responds strongly but inconsistently, which conventional differential expression cannot do.

For each dataset, all samples were projected into the principal subspace of the 800 most variable genes (five components), and each individual's response was taken as the displacement between their paired biopsies in that subspace. Two quantities were then computed per group. The *coherence* is the mean cosine between each individual's displacement and the mean displacement of their group computed without that individual; the leave-one-out construction removes the self-inclusion bias that inflates agreement in small groups. The group *direction* is the mean displacement, and the *magnitude* is the mean Euclidean norm of the individual displacements.

Significance was assessed by permutation throughout, because the sampling distribution of a cosine-based statistic in small groups is not analytically tractable. Comparisons between groups permute individual labels between the two groups being compared, 2,000 times, and the null distribution is shown in the figures rather than summarised. Comparisons before and after an intervention permute the two time points within each participant, 3,000 times, which preserves the pairing. Per-gene group-by-insulin interactions permute group labels 5,000 times and are corrected by the Benjamini–Hochberg procedure across the nine clock-output genes tested, which were fixed in advance of the analysis.

Where a dataset contained a technical factor confounded with group, the analysis was repeated on the subset in which the factor is balanced, and both results are reported; this applies to hybridisation batch in GSE22309. Power for the intervention analyses was computed for a two-sided test at α = 0.05 assuming the between-individual standard deviation of coherence observed across these cohorts.

## Differential expression and gene sets

*Purpose.* To provide the canonical read-out against which the response-based analysis is compared, using the analysis the field would apply by default.

Differential expression at rest used a two-sided Welch *t*-test per gene on expression residualised for the covariates of each study, with Benjamini–Hochberg correction across all genes tested. Association with insulin sensitivity used Spearman correlation between residualised expression and the clamp M-value.

Gene sets were curated rather than taken from a database, because the question is whether specific, mechanistically defined programmes respond, and a scan over thousands of overlapping terms would not answer it. The sets are listed in Supplementary Table 3. Enrichment is the mean absolute statistic of the set compared with 2,000 random sets of the same size drawn from the same data, which controls for set size and for the overall distribution of the statistic in that dataset; *P* is the proportion of random sets reaching the observed value.

## Reporting and reproducibility

Exact sample sizes, statistics and *P* values are given in the figure legends. Random seeds are fixed in the code. Analyses were run with Python 3.10 (NumPy, SciPy, pandas, scikit-learn) and R 4.3 (WGCNA, limma, GEOquery); versions are pinned in the repository.

# Data availability

All data analysed are public. GEO accessions are listed in Methods and in the data documentation of the code repository. GTEx v11 expression and open-access annotations are from the GTEx portal. The adipose clamp data of Rydén et al. are from the export accompanying that publication. No new data were generated.

# Code availability

All code is available at https://github.com/Danpc11/t2d_landscape, including the preparation of every cohort, the network, sheaf and landscape implementations, the response-geometry statistics, the classical differential-expression comparison, and the scripts that generate every figure from the resulting tables.

# References

1. DeFronzo, R. A. & Tripathy, D. Skeletal muscle insulin resistance is the primary defect in type 2 diabetes. *Diabetes Care* **32**, S157–S163 (2009).
2. Richter, E. A., Bilan, P. J. & Klip, A. A comprehensive view of muscle glucose uptake: regulation by insulin, contractile activity, and exercise. *Physiol. Rev.* **105**, 1867–1945 (2025).
3. James, D. E., Stöckli, J. & Birnbaum, M. J. The aetiology and molecular landscape of insulin resistance. *Nat. Rev. Mol. Cell Biol.* **22**, 751–771 (2021).
4. Højlund, K. et al. [Selective and heterogeneous defects downstream of Akt in insulin-resistant human muscle — to be completed with the authors' preferred reference].
5. Mootha, V. K. et al. PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. *Nat. Genet.* **34**, 267–273 (2003).
6. Patti, M. E. et al. Coordinated reduction of genes of oxidative metabolism in humans with insulin resistance and diabetes: potential role of PGC1 and NRF1. *Proc. Natl Acad. Sci. USA* **100**, 8466–8471 (2003).
7. Wu, X. et al. The effect of insulin on expression of genes and biochemical pathways in human skeletal muscle. *Endocrine* **31**, 5–17 (2007).
8. Coletta, D. K. et al. [Effect of acute physiological hyperinsulinaemia on gene expression in human skeletal muscle in vivo — GSE9105; to be completed].
9. Fuentes, E. N. et al. [Blunted induction of NR4A receptors in insulin-resistant human muscle — to be completed].
10. Tuvia, N. et al. Insulin directly regulates the circadian clock in adipose tissue. *Diabetes* **70**, 1985–1999 (2021).
11. Crosby, P. et al. Insulin/IGF-1 drives PERIOD synthesis to entrain circadian rhythms with feeding time. *Cell* **177**, 896–909 (2019).
12. Dyar, K. A. et al. Muscle insulin sensitivity and glucose metabolism are controlled by the intrinsic muscle clock. *Mol. Metab.* **3**, 29–41 (2014).
13. Wefers, J. et al. Circadian misalignment induces fatty acid metabolism gene profiles and compromises insulin sensitivity in human skeletal muscle. *Proc. Natl Acad. Sci. USA* **115**, 7789–7794 (2018).
14. Gabriel, B. M. et al. Disrupted circadian oscillations in type 2 diabetes are linked to altered rhythmic mitochondrial metabolism in skeletal muscle. *Sci. Adv.* **7**, eabi9654 (2021).
15. Bang, [initials] et al. [Personalised phosphoproteomic signatures of insulin resistance in human muscle — to be completed].
16. Levate, G. et al. Molecular signatures of skeletal muscle insulin resistance: bringing personalised diabetes treatment a step closer. *Signal Transduct. Target. Ther.* (2025).
17. Yoshino, M. et al. Nicotinamide mononucleotide increases muscle insulin sensitivity in prediabetic women. *Science* **372**, 1224–1229 (2021).
18. Rydén, M. et al. [Adipose tissue transcriptome before and after bariatric surgery with hyperinsulinaemic clamp — Cell Reports 2016; to be completed].
19. Gachon, F. et al. [PAR-bZip transcription factors DBP, TEF and HLF as clock output — to be completed].
20. Vainshtein, A. et al. The hallmarks of skeletal muscle health. *Nat. Metab.* **8**, 1843–1870 (2026).

> **Note on references.** Entries marked "to be completed" are those whose exact bibliographic details could not be verified from within the analysis environment and must be checked against the primary source before submission; the content attributed to each is accurate but the citation string is provisional. All other entries were verified.

\newpage

# Figure legends

**Fig. 1 | Metabolic organs share transcriptional architecture and individual position, without robust separation by glycaemic stage.**
**a**, A coexpression network of the kind the sheaf is built on: the 150 most variable genes of healthy skeletal muscle (GSE25462, *n* = 15), with edges drawn for the strongest 3% of the biweight-midcorrelation adjacency raised to the soft-thresholding power β = 6, node size proportional to connectivity and colour given by hierarchical clustering into four modules. The sheaf does not use the drawing: it uses the spectral subspace of this network's normalised Laplacian. **b**,**c**, Network organisation by glycaemic stage in the four discovery cohorts, computed on the 800 genes shared by all organs at equal *n* per stage (subsampling, 20 repetitions; mean ± 95% CI). Global efficiency (**b**) falls and effective resistance (**c**) rises from healthy to T2D in islet (0.067 → 0.053 and 0.69 → 0.77) and adipose tissue (0.163 → 0.156 and 0.068 → 0.071); muscle is flat and liver, with four samples per group, is uninformative. **d**, The cellular sheaf used throughout. Each organ's coexpression network contributes a stalk, the spectral subspace of its normalised Laplacian (*r* = 8 components); edges carry orthogonal alignments fitted by Procrustes to a neutral reference; the energy is the mean squared Grassmann distance between aligned stalks and is zero when the organs share one architecture. **e**, Sheaf energy by glycaemic stage in three discovery cohorts (islet GSE76895, muscle GSE18732, adipose GSE27951; 800 shared high-variance genes, equal *n* per stage, 20 repetitions). Circles, observed; squares, mean of the gene-correspondence null (*B* = 200); error bars, 95% CI of the null. z = −10.6, −7.5 and −7.0 (*P* = 0.005 each, the smallest value attainable with *B* = 200); stage-to-stage differences are not significant (*P* = 0.15–0.53). **f**, The same statistic in GTEx v11 donors with muscle, subcutaneous adipose tissue and pancreas from the same individual (*n* = 253), after residualising age, sex, Hardy classification, RNA integrity, ischaemic time and batch. Grey, null (*B* = 200); blue line, observed (E = 0.95 versus 1.91 ± 0.04 s.d.; z = −26). **g**, Position of each donor along the canonical variate shared by muscle and adipose tissue, as a rank traced across the three organs (one line per donor, coloured by rank in muscle). Left, observed; right, donors shuffled between organs. First canonical correlation 0.93, 0.91 and 0.89 against 0.23 ± 0.04 under shuffling (*P* = 0.002, 500 permutations); in 19 living surgical patients with paired adipose depots (GSE20950), 0.93 against 0.59 ± 0.13 (*P* = 0.0005). **h**, Effect of technical adjustment in GTEx: |z| of the sheaf energy (blue, left axis) and the muscle–adipose canonical correlation (orange, right axis) before and after adding RNA integrity, ischaemic time and batch (z from −31.2 to −26.1; ρ 0.929 to 0.932). **i**, Cross-validated R² of whole blood predicting the same donor's organ state (ten blood components predicting five organ components, five-fold, *n* = 224): 0.30, 0.23 and 0.23. **j**, Number of basins of the quasi-potential above the bootstrap stability threshold in each cohort at three kernel bandwidths (0.7×, 1.0× and 1.4× the Silverman rule); no cohort retains two basins across bandwidths. **k**, Covariate-adjusted embedding of 77 islet donors (GSE50244) coloured by HbA1c stratum; Cramér's V between mixture component and stratum is 0.17. **l**, Pearson *r* between each gene's residual expression and the canonical variate shared by muscle and adipose tissue, for the 28 genes with the largest effect of consistent sign in all three organs.

**Fig. 2 | Resting geometry does not separate disease stages; insulin reveals tissue-specific response defects.**
**a**, Cohorts analysed per organ, at rest and with biopsies before and during a hyperinsulinaemic–euglycaemic clamp; ranges are *n* per stage (rest) or per group (clamp). **b**, Sheaf energy between the three organs by stage, as in Fig. 1e: E = 1.75, 1.80 and 1.81, with no significant stage-to-stage difference (*P* = 0.15, 0.27 and 0.53; permutation of stage labels, *B* = 500). **c**, Critical-transition index by stage, one line per cohort, coloured by organ; no peak at the intermediate stage in any organ (*P* ≥ 0.19, permutation). **d**, Ratio of state dispersion in T2D relative to healthy with the 95% CI of 2,000 bootstrap resamples at equal *n*; every interval crosses 1. **e**, Filled squares, ratio of mean response magnitude (resistant or diabetic over sensitive); open circles, 1 + the difference in coherence between the same groups. Adipose tissue loses magnitude and keeps direction; muscle keeps magnitude and loses direction. **f**, Summary per organ; muscle is the only organ with paired insulin biopsies in all three states. **g**, Canonical read-out for comparison: genes differentially expressed at rest at FDR < 0.1 (two-sided Welch *t*-test on residualised expression, Benjamini–Hochberg over all genes). Islet GSE164416 (*n* = 14/34: 1,079 of 14,077), islet GSE50244 (*n* = 39/11: 451 of 15,428), muscle GSE25462 (*n* = 15/10: 2 of 48,496 probes), adipose METSIM (*n* = 142/68: 7,332 of 13,718). **h**, Set enrichment of the same contrasts: z of the mean |*t*| of each curated set relative to 2,000 random sets of equal size; asterisks mark *P* < 0.05. These resting contrasts do not measure the paired response to insulin, where the muscle signal is most evident (Figs. 4 and 5).

**Fig. 3 | Resting muscle shows limited gene-level associations with insulin sensitivity.**
**a**, Design (GSE182120): resting vastus lateralis biopsies from 49 individuals in two laboratories, 24 with normal glucose tolerance and 25 with type 2 diabetes, all with clamp-measured insulin sensitivity (M-value 34.0 ± 15.8 versus 15.3 ± 9.3 mg kg⁻¹ min⁻¹, mean ± s.d.); expression residualised on laboratory, age and BMI. **b**, Quantile–quantile plot of the genome-wide association with the M-value (Spearman's ρ, 21,595 genes): no gene at FDR < 0.1; 14 genes at *P* < 0.01 where ~215 are expected by chance; known markers highlighted (PPARGC1A ρ = 0.33, *P* = 0.021). **c**, Volcano plot of T2D versus normal glucose tolerance at rest; no gene at FDR < 0.1; dashed line, nominal *P* = 0.05. **d**, Mean |ρ| with the M-value for curated gene sets (size in parentheses) against 2,000 random sets of equal size (grey, central 95%); only the oxidative programme separates from the null (*P* = 0.012). **e**, Spearman's ρ with the M-value for the clock-output genes of Figs. 4 and 5; all *P* > 0.2.

**Fig. 4 | The healthy response to insulin is coordinated across individuals and includes late modulation of clock-output genes.**
**a**, Definition of the measure: displacement is the difference between the insulin and basal biopsy in the principal subspace of the 800 most variable genes; coherence is the mean cosine of each displacement to the leave-one-out group mean. **b**, Displacements of 20 insulin-sensitive individuals after 4 h of insulin (GSE22309) projected onto the group mean direction and its principal orthogonal direction; coherence 0.77, mean magnitude 16.3, net group response 13.1. **c**, Coherence against time after the stimulus in four independent healthy datasets: 0.20 at 30 min (GSE9105, *n* = 12), 0.24 at 1 h after a mixed meal (GSE231509, *n* = 7), 0.19 at 2 h (GSE7146, *n* = 6), 0.73 at 4 h (GSE9105) and 0.77 at 4 h (GSE22309, *n* = 20); the connecting line is not a within-person time course. **d**, Paired *t* per gene at 30 min versus 4 h in 12 healthy men (GSE9105); blue, the 55 genes of the replicated healthy programme (|*t*| > 4 in the insulin-sensitive group of GSE22309 and the same sign with |*t*| > 3 at 4 h in GSE9105). **e**, Gene sets in the healthy 4-h response (GSE22309, insulin-sensitive): mean |paired *t*| against 2,000 random sets of equal size (grey, central 95%). **f**, Paired *t* for clock-output genes at 30 min and 4 h (GSE9105) and at 4 h in GSE22309. **g**, Control for the repeat-biopsy artefact: mean paired *t* of immediate-early genes (0.5, 4.9 and 5.8) and of clock-output genes (−0.1, 1.8 and 3.2) at each time point.

**Fig. 5 | Insulin resistance preserves response magnitude but reduces directional coherence and alters clock-output gene responses.**
**a**, Magnitude of each individual's displacement by group (GSE22309; log scale; bars, mean ± s.d.): 16.3 (sensitive, *n* = 20), 18.5 (resistant, *n* = 20) and 13.5 (T2D, *n* = 15); no group differs from the sensitive group (*P* > 0.4, permutation). **b**, Loss of coherence relative to the sensitive group (0.42 and 0.32) against the distribution from permuting group labels (violins, 2,000 permutations); *P* = 0.007 and *P* = 0.037. **c**, The same comparison restricted to the 35 individuals whose two biopsies were processed in the same batch (filled; 11, 13 and 11 per group) alongside all pairs (open): effect sizes preserved (0.76, 0.45 and 0.46), *P* = 0.10 and *P* = 0.19; the direction contrast remains significant (cosine 0.14 against a null of 0.62, *P* = 0.023). **d**, Direction of each person's response relative to the healthy mean direction, which points north, as a rose of angles (20° bins, symmetrised about the vertical axis); group cosine between sensitive and T2D 0.20 against 0.85 under permutation (*P* = 0.001). **e**, Coherence across all muscle datasets analysed, including two arms of a prediabetes trial at baseline (GSE157988): 0.77 and 0.73 (healthy or sensitive), 0.35 (resistant), 0.42 and 0.02 (prediabetic arms) and 0.45 (T2D); *n* beside each bar. **f**, Response to insulin of each clock-output gene by group, as absolute response on the radius (dark, repressed; light, induced). Asterisks, FDR < 0.05 for the group × insulin interaction versus the sensitive group; daggers, FDR < 0.1 for T2D (permutation, 5,000; Benjamini–Hochberg over nine genes). Resistant versus sensitive: DBP *P* = 0.003, FDR = 0.013; PER2 *P* = 0.004, FDR = 0.013; NR1D2 *P* = 0.012, FDR = 0.024; BHLHE40 *P* = 0.047, FDR = 0.070. T2D versus sensitive: DBP, TEF, HLF, NR1D2 and BHLHE40 at *P* = 0.022–0.045, all FDR = 0.054. **g**, Effect of chronic insulin on selected clock-output genes and TXNIP in primary myotubes from donors with normal glucose tolerance (*n* = 7) and type 2 diabetes (*n* = 5) (GSE182117); DBP *t* = −3.3 versus −1.7. **h**, Amplitude of the 24-h oscillation in the same myotubes (cosinor fit per donor, control arm; mean ± s.d.); ~30% lower in T2D cells (DBP 1.03 ± 0.31 versus 0.71 ± 0.18; *P* ≈ 0.15, two-sided Mann–Whitney). **i**, Fate of the 55-gene healthy programme: kept (same sign, |*t*| > 2), lost (|*t*| < 2) or inverted; 38% lost in resistant and 42% in T2D muscle, essentially none inverted.

**Fig. 6 | Adipose response magnitude is reduced; recovery of coordination after intervention remains inconclusive.**
**a**, Magnitude (grey, left axis) and coherence (blue, right axis) of the adipose response to a clamp in 23 non-obese individuals, 23 obese individuals and the same 23 women two years after bariatric surgery (Rydén et al.): magnitude 7.4, 5.3 and 6.8; coherence 0.37, 0.19 and 0.14. The loss of coherence with obesity does not reach significance (*P* = 0.11, permutation) and surgery does not restore it (*P* = 0.71, paired permutation). **b**, The same quantities in an independent adipose clamp cohort (GSE26637; 5 sensitive and 5 resistant): magnitude 28.4 versus 15.6, coherence 0.83 versus 0.63; the resistant group keeps the direction of the sensitive group (cosine 0.95). **c**, Magnitude ratio in the two adipose cohorts and in muscle: 0.71 and 0.55 versus 1.14. **d**, Coherence before and after 10 weeks of nicotinamide mononucleotide or placebo in prediabetic women (GSE157988; 11 and 12 participants); neither arm changes (*P* = 0.77 and *P* = 0.82, paired permutation, 3,000). **e**, Change in coherence detectable with 80% power as a function of group size (two-sided, α = 0.05, between-individual s.d. 0.3); dashed lines mark the two intervention cohorts. **f**, Summary: healthy organs respond along one shared direction; insulin-resistant muscle responds in scattered directions and insulin-resistant adipose tissue in the same direction with half the magnitude.

\newpage

# Supplementary Information

## Supplementary Note 1 | Why the sheaf, and what it adds over simpler comparisons

A simpler way to ask whether two organs are organised alike is to correlate their gene–gene correlation matrices directly. That comparison fails for two reasons. First, it is dominated by genes with high expression in both tissues, which are largely housekeeping genes, so a positive answer would be uninformative. Second, it has no natural null: permuting genes breaks both the within-organ and the between-organ structure, so a significant result cannot be attributed to the correspondence between organs. The spectral formulation solves both problems. It compares subspaces rather than entries, which removes the dependence on individual gene pairs; and it admits a null that destroys only the correspondence between organs, leaving each organ's own structure intact. The orthogonal alignment is required because the eigenbasis of a Laplacian is defined only up to rotation within degenerate subspaces, so an unaligned comparison would report differences that are purely conventional.

The cost of this formulation is interpretability at the level of single genes, which is why the gene-level counterpart is reported separately (Fig. 1l) by projecting genes onto the shared individual variate.

## Supplementary Note 2 | The missing organ: a proposal for islet

No in vivo perturbation with paired sampling exists for human pancreatic islet, which is why islet appears in the resting analyses but not in the response analyses. The closest equivalent is an ex vivo perturbation with a paired control from the same donor preparation. Two public datasets have that design: islets from 26 non-diabetic donors exposed to palmitate, high glucose or their combination with a paired control and a four-day washout (GSE159984), and islets from 13 donors exposed to palmitate or control medium for 48 h (GSE53949). The first is preferable because the paired design is per preparation, several stimuli are available, and the washout provides a recovery arm comparable to the bariatric-surgery and nicotinamide-mononucleotide analyses in Fig. 6. Applying the coherence statistic to these datasets would test whether the loss of coordination is specific to muscle or is a general property of perturbed metabolic tissue; we regard this as the natural next analysis rather than a claim of the present work.

## Supplementary Table 1 | Cohorts analysed

Accession, organ, platform, design, number of individuals per group, phenotype available, and the analysis in which each cohort is used. [Generated from `docs/REPLICATION.md`.]

## Supplementary Table 2 | Claims tested and their outcome

Every claim examined during the analysis, with the test applied and the outcome: retained, weakened or withdrawn. Four claims were withdrawn after re-testing and do not appear in the manuscript: a critical window at the prediabetic stage based on a Fisher–Rao peak, which did not survive a null that preserves the sampling density; a contraction of the occupied state space with disease, which reversed in two islet cohorts; a blood-based read-out of disease label, which showed no signal in three cohorts; and immune decoupling between organs, which did not survive correction. [Generated from `docs/AUDIT.md`.]

## Supplementary Table 3 | Curated gene sets

Membership of the ten gene sets used in Figs. 2h, 3d and 4e, with the rationale for each set.

## Supplementary Table 4 | Response geometry per dataset and group

Magnitude, coherence (leave-one-out and naive), net group response and *n* for every dataset and group analysed. [Generated from `results/response/response_geometry_by_group.tsv`.]

## Supplementary Table 5 | Per-gene clock-output interactions

Response per group, interaction *P* and FDR for the nine clock-output genes in each comparison. [Generated from `results/response/GSE22309_clock_interaction.tsv`.]

## Supplementary Figure 1 | Audit of withdrawn claims

Observed Fisher–Rao peak against the density-preserving null; state-space dispersion by stage in every cohort; longitudinal dispersion with and without batch adjustment; depth of blood samples relative to the healthy reference in three cohorts.

\newpage

# Figures

![](figs/Fig1_systemic_state_visual.png)

**Fig. 1**

\newpage

![](figs/Fig2_tissue_search_visual.png)

**Fig. 2**

\newpage

![](figs/Fig3_resting_blind_visual.png)

**Fig. 3**

\newpage

![](figs/Fig4_healthy_coordinated_visual.png)

**Fig. 4**

\newpage

![](figs/Fig5_IR_fragments_visual.png)

**Fig. 5**

\newpage

![](figs/Fig6_adipose_interventions_visual.png)

**Fig. 6**
