# Supplementary Information

## Supplementary Note 1 | Comparing the organisation of different organs

Correlating the gene–gene correlation matrices of two organs directly is dominated by genes highly expressed in both, and admits no null that isolates the between-organ correspondence. Representing each organ by the leading spectral subspace of its coexpression network and comparing after orthogonal alignment solves both problems: it compares subspaces rather than entries, and admits a null that permutes gene labels within each organ, destroying only the correspondence between them. Orthogonal alignment is required because the eigenbasis of a Laplacian is defined only up to rotation within degenerate subspaces.

Applied to eleven GTEx tissues, the shared individual position is general rather than metabolic: all 30 available tissue pairs show canonical correlations far above the shuffled null, including muscle–ileum (0.87) and liver–adrenal (0.89). Metabolic pairs show it more strongly (excess over null 0.62 versus 0.47; P = 0.001), and replacing pancreas with stomach raises the sheaf energy from 0.95 to 1.21. We therefore treat this as context for the response analyses rather than as a metabolic finding.

## Supplementary Note 2 | Response geometry in human cohorts: practical notes

The statistic used here has an antecedent in single-cell CRISPR screens, where the mean cosine between each cell's displacement and the mean response direction measures perturbation stability and was found to track effect magnitude closely (ρ = 0.84–0.98). Three differences matter when the unit is a person rather than a cell.

Between people the two quantities are largely dissociated (ρ = 0.33 across 29 groups, inverted within the key cohort). Isogenic cells in one culture differ in far fewer respects than people do.

At the sample sizes of biopsy studies the naive statistic is badly biased, returning 0.41 at n = 4 in the absence of any shared direction; with thousands of cells this bias is negligible, which is why it has not been an issue in that setting. The leave-one-out form is necessary here.

A group-level statistic discards the per-individual observations. Replacing group coherence by per-person alignment, on identical data, changes P from 0.10 to 0.013.

## Supplementary Note 3 | Concentration across perturbations

| Perturbation | Time scale | Tissue | n | Coherence |
|---|---|---|---|---|
| Insulin, clamp, 4 h (sensitive) | hours | muscle | 20 | 0.77 |
| Acute exercise, before training | hours | muscle | 25 | 0.80 |
| Acute exercise, after training | hours | muscle | 23 | 0.70 |
| Acute exercise, 3 h recovery | hours | muscle | 17–20 | 0.45–0.53 |
| Acute exercise | hours | adipose | 12–14 | 0.31–0.38 |
| Training, weeks | weeks | muscle | 24 | 0.31 |
| Training, weeks | weeks | adipose | 13 | −0.09 |
| Diurnal contrast | 12 h | muscle | 11–12 | −0.07 to 0.10 |
| Cold acclimation, 10 days | days | adipose | 7 | −0.01 |
| Palmitate ex vivo, 4 days | days | islet | 5 | 0.42 |
| Insulin, 100 nM, 0.5–2 h | hours | myotubes | 12–24 | −0.25 to 0.27 |

Concentration near 0.8 is what a healthy acute response looks like in intact muscle; the acute-versus-chronic contrast is demonstrated within the same participants, tissue and platform; and the isolated myocyte reaches at most a third of the tissue value.

## Supplementary Note 4 | The islet

No in vivo perturbation with paired sampling exists for human pancreatic islet. The closest equivalent is ex vivo perturbation with a paired control per donor preparation (Extended Data Fig. 2), which shows intermediate concentration at three to five preparations per contrast and uses a stimulus that differs in kind from an acute hormonal signal.

## Supplementary Tables

1. Cohorts analysed: accession, study, organ, design, perturbation, groups, participants, platform, and the analysis each is used in.
2. Claims tested and their outcome.
3. Curated gene sets and the rationale for each.
4. Response geometry per dataset and group.
5. Per-gene clock-output interactions with P and FDR.
