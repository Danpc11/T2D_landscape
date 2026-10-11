# Glucocorticoid tone is a prominent component of shared transcriptional variation between human metabolic tissues

## Abstract

**Aims/hypothesis** Transcriptional studies of insulin resistance examine one tissue at a time. We asked what is shared between the transcriptomes of different tissues of the same person, whether that shared component has the structure of a co-expressed programme, and whether it is metabolically informative.

**Methods** We analysed paired subcutaneous adipose and skeletal muscle transcriptomes from 232 non-diabetic adults (AAGMEx), 28 women (GSE41168) and 647 GTEx donors, extended to 20 GTEx tissues. Gene sets were specified from published literature before analysis and were never derived from the data under test. Set-level statistics were compared against random gene sets matched for size and, where stated, for expression level and variability. Cohorts were combined by random-effects meta-analysis [13] of Fisher-transformed correlations with bootstrap standard errors. Co-expression networks [14] were built on cross-tissue concordance and, separately, within each tissue using identical soft-thresholding and the same 232 participants. The frozen signature was tested against a within-person, placebo-controlled hydrocortisone challenge (n=10). Cross-tissue prediction was cross-validated with all learning confined to training folds and disjoint genes in predictor and outcome.

**Results** Conventional differential expression varied by three orders of magnitude between cohorts, and of the genes detected in the three cohorts with appreciable signal only six were shared by all three; in the one cohort with clamp-measured sensitivity in every participant, 739 associated genes fell to none on adjustment for age or BMI alone. Shared variation between a person's tissues separated into two components of different nature. A genetic component — genes silenced by common structural variants — gave a pooled ρ of 0.807 (95% CI 0.774–0.835, I²=0%), essentially identical in two independent living cohorts. Physiological components were weaker: clock output 0.272 (0.217–0.325, I²=0%), insulin-regulated core 0.197 (0.140–0.253, I²=0%) and a canonical glucocorticoid signature 0.219 (0.088–0.343, I²=67%), led by ZBTB16 and FKBP5. Against size-matched null sets all three were significant (P ≤ 0.007); against null sets additionally matched on expression level and variability, which have twice the background concordance, the genetic component (P = 0.002) and clock output (P = 0.014) remained significant whereas the glucocorticoid signature (P = 0.100) and insulin core (P = 0.116) did not. The core oscillator and the oxidative programme were formally null. The glucocorticoid component did not form a co-expression module, and its genes sat at the 0.22 connectivity percentile within their own tissue networks, whereas the oxidative programme behaved as a genuine module (0.60–0.69). The component was present in 17 of 20 GTEx tissues but absent in the adrenal gland (ρ with the axis 0.08) and decoupled in liver, which has the highest HSD11B1 expression. A frozen score rose under hydrocortisone (Δ=+0.58, P=0.0002), as did the insulin-regulated core (Δ=+0.48, P=0.0002). In GTEx, circumstance of death shifted the axis by more than one SD.

**Conclusions/interpretation** Shared transcriptional variation between human metabolic tissues comprises a genetic component of technical rather than physiological interest and a glucocorticoid-responsive programme whose genes are coordinated by a circulating signal rather than by a shared transcriptional programme. The programme is readable from a subcutaneous fat biopsy but not from blood. Post-mortem cohorts inflate it, and estimates from them should not be transferred to living people without adjustment.

**Keywords** Adipose tissue · FKBP5 · Glucocorticoids · GTEx · Insulin resistance · Multi-tissue transcriptomics · Skeletal muscle

## Research in context

**What is already known about this subject?**
- Glucocorticoid excess causes insulin resistance, and FOXO3 integrates insulin and glucocorticoid signalling onto a shared catabolic programme in skeletal muscle.
- Transcriptional dysregulation associated with gluco-metabolic traits has been mapped separately in human adipose tissue and muscle.
- Inhibitors of 11β-HSD1 reduce intracellular cortisol and have been tested in type 2 diabetes.

**What is the key question?**
- What is shared between the transcriptomes of different tissues of the same person, and is that shared component metabolically informative?

**What are the new findings?**
- Shared variation has two components of different nature: a near-perfectly concordant genetic one (ρ ≈ 0.80, reproducible to two decimal places across cohorts) and a moderate physiological one comprising clock output and glucocorticoid signalling.
- The glucocorticoid component is absent in the adrenal gland and decoupled in liver, the two organs with autonomous cortisol handling.
- The circumstance of death shifts the axis by more than one standard deviation in GTEx.

**How might this impact on clinical practice in the foreseeable future?**
- A subcutaneous fat biopsy, but not blood, provides a readable tissue-level measure of this programme, which could complement metabolic endpoints when evaluating pharmacodynamic engagement of drugs acting on glucocorticoid signalling.

---

## Introduction

Insulin resistance involves several tissues, yet their transcriptional alterations are usually studied one at a time. Examining what is shared between organs places individual differences in a systemic context. Resting profiles, however, do not show how a tissue responds to a physiological demand, so connecting the study of transcriptional organisation between tissues with the analysis of individual responses to stimulation is necessary.

That connection matters because the disease is defined by a response. Diagnosing insulin resistance requires infusing the hormone and measuring glucose disposal [1, 2], whereas the molecular biology that describes it comes almost entirely from fasting biopsies. After two decades of resting profiles, the only reproducible finding in muscle is reduced expression of the PGC-1α-associated oxidative programme [3, 4].

Three reasons make skeletal muscle and subcutaneous adipose tissue the pair to examine. Together they account for most insulin-stimulated glucose disposal and most regulated lipid storage. They are the only two human tissues for which paired biopsies under a defined perturbation exist in more than one cohort. And they are the two for which transcriptome and clamp-derived insulin sensitivity have been measured in the same individuals.

We asked a question that the single-tissue literature does not pose: among the thousands of genes expressed in both tissues, which ones place a given person consistently relative to others in both of their organs at once? The answer, we show, has two components of entirely different nature.

---

## Methods

### Cohorts

[Table 1 lists all cohorts.] The principal analyses use three sources of paired muscle and subcutaneous adipose transcriptomes from the same individuals. AAGMEx [6] (GSE95674 adipose, GSE95675 muscle) comprises 256 and 255 samples from non-diabetic African American adults aged 18–60; 232 individuals have both tissues, biopsied at the same visit after an overnight fast on Illumina HumanHT-12 v4 arrays. GSE41168 comprises 28 non-obese women with normal glucose tolerance with both tissues sampled before intervention. GTEx v11 [5] provides 647 donors with both tissues and, for the tissue survey, 20 tissues with at least 40 donors each.

The pharmacological test uses GSE148640: ten patients on hydrocortisone replacement, each sampled in adipose tissue under hydrocortisone and under placebo (within-person, placebo-controlled).

### Gene sets

All gene sets were specified from published literature before being tested and were never derived from the data under test. The glucocorticoid set comprises FKBP5, TSC22D3, ZBTB16, KLF15, DDIT4, SESN1, GLUL, SGK1 and CEBPD. The insulin-regulated core comprises TXNIP, PDK4, HK2, IRS2, ANGPTL4 and RRAD; it shares no gene with the glucocorticoid set. Clock output, core oscillator, oxidative phosphorylation and the genetic reference set (GSTM1, GSTT1, ERAP2, and in GTEx also HLA-DQA1 and HLA-DRB4) are specified in ESM Table 1.

### Cross-tissue concordance

For each gene expressed in both tissues, concordance is the Spearman correlation across individuals between its expression in muscle and in adipose tissue. The background is the mean of this quantity over 2,500–3,000 randomly drawn genes. Set-level significance is the proportion of 3,000–5,000 random sets of the same size whose mean concordance equals or exceeds the observed value. Where stated, random sets were additionally matched to each target gene on mean expression and standard deviation, by sampling from the 80 nearest genes in that two-dimensional space.

### Axis score and tissue survey

Within each tissue, expression was residualised on RNA integrity number and ischaemic time, then standardised per gene; the axis score is the mean of the standardised glucocorticoid genes. Internal coherence is the mean pairwise Spearman correlation among the set's genes within a tissue. For GTEx analyses of demographic associations, age and sex were not residualised, being the variables of interest; death circumstance (Hardy scale) was examined both as a covariate and as a stratifier.

### Cross-tissue prediction

Five-fold cross-validation with ridge regression. All learning — covariate adjustment, standardisation and model fitting — was confined to the training folds, and both tissues of a donor were always assigned to the same fold. The predictor comprised the nine glucocorticoid genes measured in adipose tissue; the outcome comprised eight different genes measured in muscle, with no gene appearing on both sides.

### Statistical analysis and reproducibility

Benjamini–Hochberg correction [12] is applied as the reverse cumulative minimum over p-values sorted ascending. Enrichment against MSigDB [11] v2026.1.Hs used the same random-set procedure rather than a competitive test against the whole genome, because the analysed gene list is already selected for cross-tissue concordance. Code and intermediate results are available at [repository].

---

## Results

### The conventional analysis, and its limit

We first applied standard differential expression across the cohorts (Fig. 1c). The yield varies by three orders of magnitude between studies, from 6,930 genes at FDR < 0.1 in METSIM adipose tissue to none in two muscle cohorts, and the overlap between studies is slight: of the genes detected in the three cohorts with appreciable signal, only six are shared by all three (Fig. 4b).

In the one cohort with clamp-measured insulin sensitivity in every participant, the unadjusted association yields 739 genes at FDR < 0.1, with 85 of 211 mitochondrial genes among them, reproducing the published enrichment [3, 4]. The association disappears entirely on adjustment for age or for BMI, either one alone (Fig. 1d). BMI itself has no associated gene, so the signal arises from variance shared between adiposity, age and sensitivity rather than from a signature of adiposity. With partial correlation and both sides residualised, no gene of 21,595 survives at FDR < 0.1.

### The pipeline recovers published findings

Before asking a different question we verified that the analysis recovers what is already known (Fig. 2a). Insulin induction of HK2 is blunted in insulin-resistant participants of GSE13070 (+656 vs +165, P = 0.012), as reported for that cohort [7]. ANGPTL4 is concordant between adipose tissue and muscle in AAGMEx (ρ = 0.33, P = 0.002), as reported for that cohort [18]. The oxidative programme is associated with insulin sensitivity before adjustment. Pathway enrichment of the differential expression recovers the expected themes per cohort (Fig. 2b).

### Shared variation between a person's tissues has two components

Among 34,675 genes measured in both tissues of 232 individuals, cross-tissue sharing is rare: most genes show none (background ρ = 0.037), 1.4% exceed ρ = 0.5 and 0.2% exceed ρ = 0.7 (Fig. 1e, Fig. 3c).

Against that background, six pre-specified sets behave very differently (Fig. 1f, Fig. 3a). In a random-effects meta-analysis across the two living cohorts, the strongest component by far is genetic — genes silenced by common structural variants — at pooled ρ = 0.807 (95% CI 0.774–0.835) with I² = 0%: 0.808 in AAGMEx and 0.800 in GSE41168, two cohorts differing in platform, population and size. Among physiological components, clock output is pooled ρ = 0.272 (0.217–0.325, I² = 0%), the insulin-regulated core 0.197 (0.140–0.253, I² = 0%) and the canonical glucocorticoid signature 0.219 (0.088–0.343), the last with substantial heterogeneity (I² = 67%; Fig. 3b). The core oscillator is formally null (0.035, −0.175 to 0.241) as is the oxidative programme (0.039, −0.041 to 0.118).

Two null models were used. Against random sets matched only for size, all four components with non-zero pooled estimates are significant (P ≤ 0.007). Against random sets additionally matched on mean expression and standard deviation, only the genetic component (P = 0.002) and clock output (P = 0.014) remain significant; the glucocorticoid signature gives P = 0.100 and the insulin core P = 0.116. The matched null is the more conservative of the two for a measurement reason rather than a biological one: its background concordance is 0.091 against 0.044 for the size-matched null, because highly expressed and highly variable genes are better measured and therefore correlate better between tissues. It is also modestly circular — it draws known glucocorticoid targets at 6.5 times the chance rate against 1.4 times for the size-matched null — although none of the twenty genes most closely matched to FKBP5 is a known target. We report both values throughout and note that a nine-gene set has little power against a background twice as high.

The glucocorticoid component is led by ZBTB16 (ρ = 0.441) and FKBP5 (0.424) (Fig. 1 legend data; Fig. 2c–d), and in an unbiased ranking of all genes by cross-tissue concordance FKBP5 is the single highest-ranking annotated gene.

### The shared component is not a co-expression module

If the glucocorticoid component were a transcriptional programme, its genes should co-vary with one another. Building a network on cross-tissue concordance (869 genes with ρ > 0.30) yields eight modules, and only the genetically driven ones emerge (Fig. 3d,e). Module 8 contains GSTM1 and GSTM2, adjacent at 1p13, together with GSTT1, GSTT2, GSTT2B, DDT and DDTL, all at 22q11; module 3 contains ERAP2 and MICA from the extended MHC together with HPR. The network therefore groups genes by chromosomal locus and structural-variant class, an external validation that it detects real structure. The genetic set falls 2 of 3 into a 15-gene module (P = 8.8 × 10⁻⁴), whereas glucocorticoid and clock genes do not form a module (2 of 9 and 2 of 7 enter the network; P = 0.86 against chance).

Within-tissue networks built on the same 232 subjects with identical soft-thresholding in both tissues show the same thing from another angle (Fig. 3f–h). Glucocorticoid genes sit at the 0.22 connectivity percentile in adipose tissue and 0.22 in muscle; clock output at 0.10 in muscle, the insulin core at 0.11, the genetic set at 0.04–0.08. All are significantly less connected than chance. The only set behaving as a genuine co-expressed module is the oxidative programme (percentile 0.60 and 0.69). Connectivity itself correlates 0.75 between tissues, so the architecture is preserved and this is not an artefact of a poorly built network.

Genes that place an individual consistently across their tissues are therefore the ones that do *not* belong to co-expressed programmes within a tissue — consistent with coordination by a circulating signal rather than a shared transcriptional programme.

### Where the component exists, and what it represents

Across 20 GTEx tissues the glucocorticoid set shows internal coherence between 0.16 and 0.58 (Fig. 2f). Two tissues fall outside the pattern and they are the two predicted to: the adrenal gland, which produces cortisol, has the lowest coherence (0.158) and essentially no correlation with the systemic axis (ρ = 0.078 against 0.59–0.95 elsewhere); the liver, with the highest HSD11B1 expression (8.06 vs 1.3–4.4), is likewise decoupled (0.190). Coherence correlates with glucocorticoid receptor expression across tissues (ρ = +0.34).

Against expression- and variability-matched random sets, the axis is positively associated with ubiquitin-mediated degradation, muscle atrophy genes (ρ = 0.447), autophagy/FOXO (0.416), metallothioneins (0.441) and lipolysis (0.211), and negatively with extracellular matrix, respiratory electron transport and muscle contraction (Fig. 4a,c). Inflammation (0.125, P = 0.27) and the unfolded protein response are not associated. The oxidative set averages −0.023 but is not homogeneous: PPARGC1A (+0.216) and ESRRA (+0.118) are positive while respiratory chain components are negative.

Inferred transcription-factor activity ranks FOXO3 highest (ρ = 0.66 with the glucocorticoid score and 0.67 with the insulin core), followed by the androgen receptor (0.62) and CEBPB (0.55), with BCL6 negative (Fig. 2e). The FOXO3 regulon contains TXNIP (ρ = 0.72 with the axis) and FBXO32 (0.45), linking the insulin-regulated core and the atrophy programme to this factor through externally curated regulatory annotation (ESM Fig. 1). The androgen-receptor signal should be read with caution: its top-ranked targets are ZBTB16 and FKBP5, the canonical glucocorticoid genes, and the two steroid receptors bind near-identical response elements, so this most likely reflects the glucocorticoid receptor. Motif-defined target sets point the same way but with weaker support after correction over 853 sets (FOXO3_01 z = +2.61, P = 0.0067, FDR = 0.056; CEBP_C FDR = 0.095); MEF2 is negative and does survive (z = −4.85, FDR = 0.015) (Fig. 4d).

### The signature responds to a controlled glucocorticoid challenge

Applied without refitting to a within-person placebo-controlled hydrocortisone challenge in ten living patients, the frozen score rises (Δ = +0.582, P = 0.0002), as does the insulin-regulated core (Δ = +0.482, P = 0.0002) (Fig. 5a). ZBTB16 is the most induced gene of the experiment (Δ = +2.05, P = 4 × 10⁻⁸), followed by PER1, TXNIP, FKBP5 and PDK4 (Fig. 5b). Three features qualify the reading: ARNTL rises (+0.78, P = 0.007), so cortisol does move the core oscillator even though it is not part of the between-person axis; the atrophy genes do not move at this dose and interval; and clock output as a set does not move although PER1 responds strongly.

### Circumstance of death inflates the axis in post-mortem tissue

In GTEx the axis score follows a monotonic gradient with the Hardy scale, from −0.718 in violent rapid death to +0.459 in slow death, exceeding one standard deviation (Fig. 5c). Ventilator cases, half the cohort and frequent recipients of clinical corticosteroids, score +0.141. The glucocorticoid concordance in GTEx (0.478) is correspondingly higher than in either living cohort. Age shows no association (ρ = 0.016; 0.004 adjusted); women score higher than men (+0.124 vs −0.060, P = 0.001).

### A fat biopsy reports the axis; blood does not

With disjoint genes in predictor and outcome and all learning confined to training folds, adipose tissue predicts the muscle axis at cross-validated R² = 0.601; adjusting additionally for circumstance of death gives 0.505, against 0.237 for a model using age, sex and Hardy class alone (Fig. 5d). Blood predicts muscle at R² = 0.137, and the axis is less than half as coherent in whole blood as in solid peripheral tissues (Fig. 5e).

## Discussion

Shared transcriptional variation between a person's metabolic tissues is not a single phenomenon. It separates into a genetic component of near-perfect concordance and physiological components of moderate concordance, and conflating them would misstate both.

The genetic component is biologically unsurprising: a deletion polymorphism determines expression identically in every tissue. Its value is methodological. Its reproducibility across two cohorts differing in platform, population and size — 0.808 and 0.800, pooled I²=0% — makes it a calibrator. A multi-tissue dataset in which GSTT1 and ERAP2 fail to show high cross-tissue concordance has a technical problem, because that signal must be present. The same conclusion emerges without supervision: the only modules recovered from the concordance network are genetic, and they group genes by chromosomal locus and structural-variant class — GSTM1 with GSTM2 at 1p13, GSTT1, GSTT2, GSTT2B, DDT and DDTL at 22q11, ERAP2 with MICA in the extended MHC.

Among physiological components, clock output was the most consistent across cohorts (I² = 0%) and the only one besides the genetic component to survive the expression-matched null; the glucocorticoid signature was the most heterogeneous (I² = 67%) and did not survive that control. We therefore describe the glucocorticoid programme as a prominent component rather than the dominant one, and place clock output alongside it rather than beneath it; we have not quantified what fraction of total cross-tissue covariation each explains under a common procedure, and the clock result carries a design caveat discussed below.

The identity of the glucocorticoid component rests on four independent lines. The canonical signature genes are its strongest contributors, and in an unbiased ranking of all genes by cross-tissue concordance FKBP5 is the highest-ranking annotated gene. The programme responds to a controlled hydrocortisone challenge in living people, with the signature frozen beforehand. It is absent in the adrenal gland and decoupled in liver, the two organs that do not simply receive circulating cortisol — an internal validation we did not seek. And its functional associations are those of the glucocorticoid-FOXO3 axis: ubiquitin-mediated degradation, atrophy genes, autophagy and lipolysis, with extracellular matrix and respiratory electron transport opposed.

That the component is *not* a co-expression module, and that its genes are among the least connected nodes of their own tissue networks, is the observation we consider most informative. A set of genes can share the property of tracking an individual across tissues without co-varying with one another, and that is what a circulating hormonal signal predicts: each gene responds to the same external input rather than to its neighbours. The oxidative programme, which does behave as a module, provides the contrast. This also explains why a conventional WGCNA analysis finds nothing here: we were looking for a module where, by this logic, none should exist.

The circuit itself is established. FOXO3 mediates crosstalk between IRS-1/PI3K/Akt and glucocorticoid signalling and coordinates ubiquitin and atrogin-1 expression [15]; maximal glucocorticoid induction of human PDK4 requires FOXO3 and FOXO1 binding at its promoter [16]; KLF15 induces atrogin-1, MuRF1, FOXO1 and FOXO3 [17]; and ANGPTL4 is a direct glucocorticoid receptor target [19] reported concordantly associated with gluco-metabolic phenotypes in both adipose and muscle of this same cohort [6, 18]. Our contribution is not the circuit but its position: that this programme, rather than inflammation or mitochondrial dysfunction, is a prominent component of between-person covariation across organs. Inferred transcription-factor activity ranks FOXO3 highest and its curated regulon contains both TXNIP and FBXO32, but inferred activity does not establish causality, and enrichment of binding annotations over a list selected by expression is external evidence rather than independent experimental validation. The androgen-receptor signal in the same analysis is most likely the glucocorticoid receptor: its top-ranked targets are ZBTB16 and FKBP5, and the two steroid receptors bind near-identical response elements.

The Hardy gradient is a caution with implications beyond this study. A shift exceeding one standard deviation between violent and slow death means any post-mortem estimate of this axis — its amplitude, its distribution across tissues, its cross-tissue predictability — carries an agonal contribution. The axis is not an artefact: it is present in living people in two independent cohorts, and responds to a controlled challenge. But its magnitude in GTEx is not physiological, and the gap between GTEx (0.478) and the living cohorts (0.168, 0.303) quantifies that inflation. Studies using GTEx to estimate hormone-responsive programmes should stratify or adjust for death circumstance, and should note that recorded circumstance is not the same as documented corticosteroid exposure.

The clinical reading we consider defensible is narrow. Inhibitors of 11β-HSD1 reduce intracellular cortisol regeneration [22] and have been evaluated in type 2 diabetes; INCB-13739 improved HbA1c, glucose and HOMA-IR [20], which is not equivalent to clamp-measured sensitivity, and AZD4017 did not meet its glucose-disposal endpoint [21]. Our data cannot explain those outcomes and do not demonstrate hidden clinical benefit. What they show is that the programme those drugs act on is measurable at tissue level, is shared across organs, responds to a controlled glucocorticoid challenge, and is accessible from a subcutaneous fat biopsy (cross-validated R²=0.50 after adjustment for death circumstance, with disjoint genes) but not from blood (R²=0.14). A tissue transcriptional readout could therefore complement metabolic endpoints when assessing pharmacodynamic engagement and between-patient heterogeneity. The poor performance of blood acquires a mechanistic reading from published work reporting reduced basal GILZ and FKBP5 with blunted dexamethasone responses in mononuclear cells [23] of people with type 2 diabetes, consistent with cellular glucocorticoid resistance that differs between leucocytes and solid tissue.

### Limitations

Clock output was the most consistent physiological component across cohorts, but in every cohort both biopsies were taken at the same visit, so the two samples share time of day. None of the available cohorts records collection time, so we cannot separate shared circadian phase from shared sampling time, and we do not interpret the clock result mechanistically.

Inferred transcription-factor activity does not establish FOXO3 as the causal regulator, and enrichment of binding annotations over a list selected by expression is external evidence rather than independent experimental validation.

The association of the axis with type 2 diabetes requires replication outside the cohorts examined here; in 71 living individuals with measured HOMA-IR the glucocorticoid signature did not outperform random gene sets.

Finally, the insulin sensitivity phenotypes of AAGMEx are not deposited with the expression data, so we could not relate axis position to clamp- or FSIGT-derived sensitivity in the largest living cohort.


---

## References

1. DeFronzo RA, Tripathy D (2009) Skeletal muscle insulin resistance is the primary defect in type 2 diabetes. Diabetes Care 32(Suppl 2):S157–S163
2. DeFronzo RA, Tobin JD, Andres R (1979) Glucose clamp technique: a method for quantifying insulin secretion and resistance. Am J Physiol 237:E214–E223
3. Mootha VK, Lindgren CM, Eriksson KF et al (2003) PGC-1α-responsive genes involved in oxidative phosphorylation are coordinately downregulated in human diabetes. Nat Genet 34:267–273
4. Patti ME, Butte AJ, Crunkhorn S et al (2003) Coordinated reduction of genes of oxidative metabolism in humans with insulin resistance and diabetes. Proc Natl Acad Sci USA 100:8466–8471
5. GTEx Consortium (2020) The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369:1318–1330
6. Sharma NK, Chuang Key CC, Civelek M et al (2016) Tissue-specific and genetic regulation of insulin sensitivity-associated transcripts in African Americans. J Clin Endocrinol Metab 101:1455–1468
7. Sears DD, Hsiao G, Hsiao A et al (2009) Mechanisms of human insulin resistance and thiazolidinedione-mediated insulin sensitization. Proc Natl Acad Sci USA 106:18745–18750
8. Rubenstein AB, Smith GR, Raue U et al (2020) Single-cell transcriptional profiles in human skeletal muscle. Sci Rep 10:229
9. De Micheli AJ, Spector JA, Elemento O, Cosgrove BD (2020) A reference single-cell transcriptomic atlas of human skeletal muscle tissue reveals bifurcated muscle stem cell populations. Skelet Muscle 10:19
10. Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J (2019) Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Res 29:1363–1375
11. Liberzon A, Birger C, Thorvaldsdóttir H, Ghandi M, Mesirov JP, Tamayo P (2015) The Molecular Signatures Database hallmark gene set collection. Cell Syst 1:417–425
12. Benjamini Y, Hochberg Y (1995) Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc B 57:289–300
13. DerSimonian R, Laird N (1986) Meta-analysis in clinical trials. Control Clin Trials 7:177–188
14. Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9:559
15. Zheng B, Ohkawa S, Li H, Roberts-Wilson TK, Price SR (2010) FOXO3a mediates signaling crosstalk that coordinates ubiquitin and atrogin-1/MAFbx expression during glucocorticoid-induced skeletal muscle atrophy. FASEB J 24:2660–2669
16. Connaughton S, Chowdhury F, Attia RR et al (2010) Regulation of pyruvate dehydrogenase kinase isoform 4 (PDK4) gene expression by glucocorticoids and insulin. Mol Cell Endocrinol 315:159–167
17. Shimizu N, Yoshikawa N, Ito N et al (2011) Crosstalk between glucocorticoid receptor and nutritional sensor mTOR in skeletal muscle. Cell Metab 13:170–182
18. Das SK, Sharma NK, Chu AY et al (2018) Adipose tissue and skeletal muscle transcript signatures associated with gluco-metabolic phenotypes. [Complete from PMC5821540 — Transcriptional regulatory mechanisms in adipose and muscle tissue associate with composite gluco-metabolic phenotypes]
19. Chen TC, Benjamin DI, Kuo T et al (2017) The glucocorticoid-Angptl4-ceramide axis induces insulin resistance through PP2A and PKCζ. Sci Signal 10:eaai7905
20. Rosenstock J, Banarer S, Fonseca VA et al (2010) The 11-β-hydroxysteroid dehydrogenase type 1 inhibitor INCB13739 improves hyperglycemia in patients with type 2 diabetes mellitus inadequately controlled by metformin monotherapy. Diabetes Care 33:1516–1522
21. Othonos N, Marjot T, Woods C et al (2023) 11β-HSD1 inhibition in men mitigates prednisolone-induced adverse effects in a proof-of-concept randomised double-blind placebo-controlled trial. Nat Commun 14:1025
22. Morgan SA, McCabe EL, Gathercole LL et al (2014) 11β-HSD1 is the major regulator of the tissue-specific effects of circulating glucocorticoid excess. Proc Natl Acad Sci USA 111:E2482–E2491
23. [Blood glucocorticoid-resistance reference — to be completed: reduced basal GILZ/FKBP5 and blunted dexamethasone response in mononuclear cells in type 2 diabetes]

*Note: references 18 and 23 require completion from the primary sources before submission; all others were verified against the publisher record.*

---

## Figure legends

**Fig. 1 Shared transcriptional variation between a person's muscle and adipose tissue separates into a genetic and a physiological component.** (**a**) Cross-tissue concordance (Spearman ρ across individuals, for each gene, between its expression in muscle and in subcutaneous adipose tissue), averaged over six gene sets specified from published literature before analysis. Bars are grouped by cohort: GTEx (post mortem, n = 647), AAGMEx (n = 232) and GSE41168 (n = 28). Dotted lines give the background, the mean concordance of randomly drawn genes in that cohort. (**b**) The genetic component in AAGMEx, gene by gene: genes whose expression is governed by common structural variants. (**c**) The glucocorticoid signature in AAGMEx, gene by gene; grey bars fall below background.

**Fig. 2 The glucocorticoid component across 20 GTEx tissues, and its absence where cortisol is handled autonomously.** (**a**) Internal coherence of the glucocorticoid gene set within each tissue (mean pairwise Spearman correlation among its genes). Red: adrenal gland and liver, with their correlation to the systemic axis annotated. Grey: whole blood. (**b**) Coherence against HSD11B1 expression; liver regenerates intracellular cortisol and is decoupled. (**c**) Coherence against glucocorticoid receptor (NR3C1) expression.

**Fig. 3 The signature responds to a controlled glucocorticoid challenge, and is inflated by circumstance of death in post-mortem tissue.** (**a**) Change in the frozen score under hydrocortisone versus placebo, within person, in adipose tissue of ten patients (GSE148640). Asterisks: significance against random gene sets of the same size (\*p<0.05, \*\*p<0.01, \*\*\*p<0.001). (**b**) Individual genes; note that the core oscillator (ARNTL) responds although it is not part of the between-person axis, and that atrophy genes do not respond to an acute dose. (**c**) Glucocorticoid axis score by Hardy death-circumstance category in GTEx.

**Fig. 4 What programme the axis represents, and how far it can be read from one tissue to another.** (**a**) The five most positively and five most negatively associated MSigDB Hallmark sets. (**b**) Curated pathway sets tested against random sets matched for size, mean expression and variability. (**c**) Cross-validated R² for predicting the axis in one tissue from another, with disjoint genes in predictor and outcome and all learning confined to training folds.

---

## Tables

**Table 1** Cohorts analysed: accession, tissue, design, n, platform and role in the analysis. [To be completed from ESM Table 1.]

---

## Declarations

**Data availability** All datasets are publicly available under the accessions listed in Table 1. Analysis code and intermediate results: [repository, DOI].

**Funding** [To be completed.]

**Authors' relationships and activities** The authors declare no relationships or activities that might bias, or be perceived to bias, their work.

**Contribution statement** [To be completed.]
