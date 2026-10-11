# Submission kit — Diabetologia

**Title** Glucocorticoid tone is a prominent component of shared transcriptional variation between human metabolic tissues

## Contents

    01_manuscript/     MANUSCRIPT_Diabetologia.docx and .md (5,081 words)
    02_figures/        Fig1–Fig5, vector PDF and 600-dpi PNG, 180 mm wide
    03_tables/         Table 1, cohorts (TSV and Markdown)
    04_supplementary/  ESM_methods + 12 ESM tables + ESM Fig. 1
    05_code/           repository release v2.0 and result tables

## Ready

- Structured abstract (Aims/hypothesis, Methods, Results, Conclusions/interpretation)
- Research in context box
- 23 numbered references, 20 in-text citations
- Five main figures, vector PDF with embedded fonts
- Twelve supplementary tables and one supplementary figure
- Code: 25 analysis scripts, orchestrator, smoke test 15/15, no absolute paths

## Required before submission

1. **Author list, affiliations, corresponding author.** Not yet in the manuscript.
2. **Funding and contribution statements.** Placeholders present.
3. **Archive DOI for the code.** Connect the repository to Zenodo *before* creating the GitHub
   release, then cite the DOI under Data availability.
4. **Two incomplete references.** [18] Das et al., the AAGMEx transcript signature paper
   (PMC5821540). [23] the blood glucocorticoid-resistance reference (reduced basal GILZ and FKBP5
   with blunted dexamethasone response in mononuclear cells in type 2 diabetes).
5. **Two schematics.** Placeholders reserved in Fig. 1a (study design) and Fig. 5f (model).
6. **Figure font sizes.** Several panels use 4.4–5.0 pt labels; Diabetologia expects roughly 6 pt
   after reduction. Fig. 2, with eight panels including the chord diagram, is the one at risk.
   Either move its bottom row to supplementary, or reduce the chord diagram to two or three
   transcription factors so the labels can be enlarged.

## Points a reviewer will probe, and where they are addressed

- **Circularity.** Gene sets come from published literature, never from the data under test
  (Methods; ESM Table 2). Removing FKBP5 changes ρ from 0.765 to 0.759 (ESM Fig. 1a). Genes
  outside any score reproduce the association against expression-matched nulls (ESM Fig. 1b).
- **Two null models.** Both are reported throughout. The matched null has twice the background
  concordance for measurement reasons and is modestly circular (6.5× enrichment for known
  glucocorticoid targets against 1.4× for the size-matched null). Neither value is presented alone.
- **Post-mortem confounding.** The Hardy gradient exceeds one SD (Fig. 5c). Three estimates are
  reported separately: GTEx and two living cohorts. Cross-tissue prediction is given with and
  without adjustment for circumstance of death.
- **Information leakage.** All learning inside training folds, both tissues of a donor in the same
  fold, disjoint genes in predictor and outcome (Methods; Fig. 5d).
- **Clock result.** Both biopsies share time of day and no cohort records collection time; stated
  in Limitations and not interpreted mechanistically.
- **Causality of FOXO3.** Explicitly not claimed. Inferred activity and motif enrichment are
  described as external evidence, not experimental validation.

## What did not survive, and is reported as such

The original hypothesis — that insulin resistance degrades the coordination of the muscle
transcriptional response to insulin — depends entirely on one cohort. Removing GSE22309 moves b₁
from −0.58 to −0.215 and the interaction from P = 0.0023 to 0.15; GSE13070 (68 pairs, 5-h clamp)
does not replicate, and GSE24215 (resistance induced within person by bed rest) runs the other way.
It is retained as exploratory, not as a finding.
