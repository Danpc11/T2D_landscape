# Electronic Supplementary Material

## Contents

**ESM Table 1** Cohorts analysed: accession, cohort name, tissue, design, number of individuals, platform and role in the study.

**ESM Table 2** Definition of every gene set used, with its members and provenance. All sets were specified from published literature before being tested and were never derived from the data under test.

**ESM Table 3** Cross-tissue concordance of all 34,675 genes measured in both tissues (AAGMEx, n = 232). This is the primary quantity from which Figs 1, 2 and 3 derive.

**ESM Table 4** Per-cohort set-level concordance with bootstrap 95% confidence intervals and standard errors.

**ESM Table 5** Random-effects meta-analysis: pooled estimate, 95% CI, τ², I² and Cochran's Q per gene set.

**ESM Table 6** MSigDB Hallmark enrichment against the axis score, all 50 sets.

**ESM Table 7** Reactome enrichment, all evaluable pathways.

**ESM Table 8** Motif-defined transcription-factor target enrichment, all evaluable sets.

**ESM Table 9** DoRothEA inferred transcription-factor activity against the glucocorticoid and insulin-core scores, all 193 regulons of confidence A–C.

**ESM Table 10** Module assignment of the 869 genes entering the cross-tissue concordance network.

**ESM Table 11** Within-tissue connectivity of all 4,004 network genes in adipose tissue and muscle.

**ESM Table 12** Axis coherence, concordance with the systemic axis, and NR3C1 and HSD11B1 expression across 20 GTEx tissues.

**ESM Fig. 1** Circularity controls and sensitivity analyses. (a) The association between the glucocorticoid score and the insulin-regulated core is not driven by single genes: removing FKBP5, GLUL or ZBTB16 in turn changes ρ from 0.765 to 0.759, 0.738 and 0.762. (b) Eight atrophy and autophagy genes, none of which contributes to any score, give mean ρ = 0.521 with the axis against −0.012 for size-matched null sets and −0.015 for null sets additionally matched on mean expression and standard deviation (P = 0.002 against both). (c) Cell-composition adjustment using three independent single-cell references. (d) Axis coherence and concordance with the systemic axis across 20 GTEx tissues. (e) Full within-tissue connectivity distributions, with glucocorticoid genes marked.

## Additional methods

### Bootstrap confidence intervals

For each cohort and gene set, individuals were resampled with replacement 400 times; the reported interval is the 2.5–97.5 percentile range of the resampled set-level mean concordance, and the standard error its standard deviation. Resampling is over individuals, not genes, because individuals are the unit of replication.

### Random-effects meta-analysis

Set-level concordances were Fisher-transformed, with standard errors propagated by the delta method. Between-study variance τ² was estimated by DerSimonian–Laird; I² is (Q − k + 1)/Q, truncated at zero. Pooled estimates were back-transformed for reporting.

### Co-expression networks

Two networks were built. The cross-tissue concordance network used the 869 genes with cross-tissue ρ > 0.30, clustering on the correlation of their per-individual profiles averaged over the two tissues (average linkage, eight clusters). The within-tissue networks used the 4,004 most variable genes common to both tissues, with the same soft-thresholding exponent (β = 6) and the same 232 participants in each tissue, so that connectivity is comparable between them. Connectivity was expressed as a within-network percentile and compared against 3,000 random gene sets of the same size.

An earlier network analysis across 24 cohort-by-condition combinations was examined and not used: mean degree correlated with sample size at ρ = −0.92 and with the soft-thresholding exponent at ρ = −0.79, and both differed between studies, so the metrics were not comparable across cohorts.

### Analyses evaluated and not used

GSE17992 (dexamethasone for four days in muscle and adipose tissue) was evaluated and excluded for processing reasons rather than for its result: platform identifiers in the deposited annotation do not match those in the raw files, sentinel values of ±10⁵ are interleaved with the data, and dye-swap orientation is not reconstructible from the metadata. In the deposited matrix FKBP5 does not change after four days of dexamethasone, which is implausible and indicates the problem.
