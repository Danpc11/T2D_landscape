# Supplementary Information

## Supplementary Note 1 | Why the sheaf, and what it adds over simpler comparisons

Correlating the gene–gene correlation matrices of two organs directly fails for two reasons: it is dominated by genes highly expressed in both, largely housekeeping genes, and it admits no null that isolates the between-organ correspondence, since permuting genes breaks the within-organ structure too. The spectral formulation compares subspaces rather than entries and admits a null that destroys only the correspondence between organs. Orthogonal alignment is required because the eigenbasis of a Laplacian is defined only up to rotation within degenerate subspaces. The cost is interpretability at the level of single genes, which is why the gene-level counterpart is reported separately (Fig. 1l).

## Supplementary Note 2 | Relationship to directional coherence in single-cell screens

The statistic used here has an antecedent in single-cell CRISPR screens, where the mean cosine between each cell's displacement and the mean response direction has been proposed as a measure of perturbation stability, and where it was found to track effect magnitude closely (Spearman ρ = 0.84–0.98 across six datasets). Our results differ in three ways that are worth stating explicitly.

First, between people the two quantities are largely dissociated (ρ = 0.33 across 29 groups, and inverted within the key cohort, where the group with the largest displacement has nearly the lowest alignment). This is not a contradiction: isogenic cells in one culture differ in far fewer respects than people differ.

Second, at the sample sizes of biopsy studies the naive statistic is badly biased, returning 0.41 at n = 4 and 0.19 at n = 20 in the absence of any shared direction, so the leave-one-out form is necessary. With thousands of cells the bias is negligible, which is why it has not been an issue in that setting.

Third, a group-level statistic discards the per-individual observations. Replacing group coherence by per-person alignment, on identical data, changes P from 0.10 to 0.013. We recommend the per-person form for human cohorts.

## Supplementary Note 3 | What concentration looks like across perturbations

Eleven contrasts analysed with the same pipeline place the insulin result in context.

| Perturbation | Time scale | Tissue | n | Coherence |
|---|---|---|---|---|
| Insulin, clamp, 4 h (sensitive) | hours | muscle | 20 | 0.77 |
| Acute exercise, before training | hours | muscle | 25 | 0.80 (P = 0.0005) |
| Acute exercise, after training | hours | muscle | 23 | 0.70 (P = 0.0005) |
| Acute exercise, 3 h recovery | hours | muscle | 17–20 | 0.45–0.53 |
| Acute exercise | hours | adipose | 12–14 | 0.31–0.38 |
| Training, weeks | weeks | muscle | 24 | 0.31 |
| Training, weeks | weeks | adipose | 13 | −0.09 (P = 0.72) |
| Diurnal contrast, aligned | 12 h | muscle | 12 | 0.10 |
| Diurnal contrast, misaligned | 12 h | muscle | 11 | −0.07 |
| Cold acclimation, 10 days | days | adipose | 7 | −0.01 (P = 0.54) |
| Palmitate ex vivo, 4 days | days | islet | 5 | 0.42 |

Concentration near 0.8 is what a healthy acute response looks like in muscle, which places the 0.77 observed under insulin in context and shows it is not an artefact of the 2007 platform, since the exercise value comes from an independent cohort on another platform. The acute-versus-chronic contrast is demonstrated within the same participants, tissue and platform. Circadian misalignment does not provide a mechanistic link: the contrast that design measures is evening versus morning twelve hours apart, which is not a concentrated response even in the aligned condition, so the comparison is uninformative rather than negative.

## Supplementary Note 4 | The missing organ: islet

No in vivo perturbation with paired sampling exists for human pancreatic islet. The closest equivalent is ex vivo perturbation with a paired control per donor preparation (Extended Data Fig. 2). Islet responses to palmitate, high glucose and their combination show intermediate concentration (0.38 to 0.90) at three to five preparations per contrast, which neither supports nor excludes a difference, and the stimulus differs in kind from an acute hormonal signal.

## Supplementary Note 5 | Analyses attempted and not supported

Recorded so that the negative results are available. A splicing axis suggested by candidate inspection (induction of SRSF1 and SRSF5 in healthy muscle, repression of CLK1 in insulin-resistant muscle) does not survive a set-level test (P = 0.39 and 0.077 against size-matched random sets). Correlating each gene's response with a person's alignment returns half the transcriptome and is circular, since alignment is computed from the same matrix. A conjunction of thresholds identifying genes that "gain" a response is not a test and those genes are unremarkable under the proper interaction test. Transcription-factor activity inference gains power over single genes (39 regulons at FDR < 0.05 versus 2 genes) but none survives restriction to same-batch pairs, and regulons share targets so the count is not 39 independent findings.

## Supplementary Tables

1. Cohorts analysed: accession, study, organ, design, perturbation, groups, participants, platform and the analysis each is used in.
2. Claims tested and their outcome, including four withdrawn after re-testing.
3. Curated gene sets and the rationale for each.
4. Response geometry per dataset and group: magnitude, coherence, alignment, n.
5. Per-gene clock-output interactions with P and FDR.
