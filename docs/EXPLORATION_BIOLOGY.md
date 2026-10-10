# Biological exploration: what the data support, and what they do not

Record of an exploratory search for novel biology and druggable targets in the response data.
Written so that the attempts that failed are documented alongside the one that survived, because the
failures determine what the manuscript may claim.

## 1. No single-gene discovery survives genome-wide correction

A genome-wide gene × insulin interaction test (Welch on the per-person response, Benjamini–Hochberg
over 8,632 genes, GSE22309) gives:

| Comparison | FDR < 0.05 | FDR < 0.1 | *P* < 0.001 (8 expected) |
|---|---|---|---|
| Insulin-sensitive vs insulin-resistant | 2 (TBX1, CBARP) | 299 | 95 |
| Insulin-sensitive vs T2D | 1 (MARCKS) | 1 | 57 |

There is a clear excess of small *P* values — about twelve times what chance predicts — but no
individual gene is robust, and the two or three that pass are not interpretable as a mechanism. This is
the expected signature of a distributed effect and is consistent with the central claim of the paper:
the difference lives in the coordination of many genes, not in any one of them. **No gene-level target
can be nominated from these data.**

Reassuringly, the top-ranked genes are not batch artefacts: 19 of the top 20 for the resistant
comparison, and 20 of 20 for the diabetic comparison, remain nominally significant in the same-batch
subset.

## 2. The clock-gene result is a pre-specified hypothesis, not a discovery

DBP, PER2 and NR1D2 reach FDR < 0.05 only because the correction is over the nine clock-output genes
chosen in advance from the circadian literature. Genome-wide they do not stand out. This is a
legitimate targeted test and must be reported as such: it is a confirmation of a prior hypothesis in a
new setting, not a screen result. Any text implying that the clock genes emerged from the data is wrong.

## 3. Attempts that failed, and why

**Splicing axis.** Candidate inspection suggested that insulin induces SRSF1 and SRSF5 in healthy
muscle while CLK1 is repressed in insulin-resistant muscle, which would be attractive because CLK1
phosphorylates SR proteins and has inhibitors in clinical development. It does not hold up. A set-level
test of 19 splicing and SR-kinase genes gives *P* = 0.39 (resistant) and *P* = 0.077 (diabetic) against
size-matched random sets. The individual *P* values (SRSF5 0.017, CLK1 0.018, RSRP1 0.0004) are what
one expects from 8,632 genes without correction. **Discarded.**

**Correlating gene responses with individual alignment.** Correlating each gene's response with how
closely that person's displacement aligns with the healthy direction gives 4,323 of 8,632 genes at
FDR < 0.05. This is circular: the alignment is computed from the same expression matrix, so the
correlation is arithmetic, not biological. **Discarded**, and the approach should not be used.

**Four genes that "gain" a response.** A conjunction filter (|t| < 2 in sensitive, |t| > 4 in
resistant, |t| > 3 in diabetic) returns CTAGE5, RSRP1, CLK1 and NR1D2. A conjunction of thresholds is
not a test; under the proper interaction test these genes are unremarkable. **Discarded as a finding**,
retained only as a description.

## 4. What does survive, and what it implies therapeutically

Two descriptive observations at the level of the replicated 55-gene healthy programme, which are
group-level comparisons and not circular:

- **What is lost** in both resistant and diabetic muscle includes PRNP, DNAJA1, BZW1, IER2, HOMER1,
  GOSR1, PIK3IP1, COL6A1 and LTBP4: chaperones, immediate-early genes, a Golgi SNARE involved in
  vesicle trafficking, a negative regulator of PI3K, and matrix genes whose repression by insulin fails.
- **What persists** in all three groups: MYOD1, VEGFA, MT1X, BHLHE40, RRAD, DNAJB5, ETS2, KLF10, JUNB.

These lists are hypothesis-generating only; none of their members survives genome-wide correction
individually.

**Therapeutic implication.** The data do not nominate a drug target, and the manuscript should not
imply one. What they support is a statement about *strategy rather than molecule*: if what fails is the
coordination of a response that takes hours to assemble and involves clock output, then interventions
that act on the timing of metabolic signals — meal timing, the timing of insulin administration, sleep
and circadian alignment — are the class worth testing, and the coherence of the response is an endpoint
that current trials do not measure. Two interventions analysed here, bariatric surgery and nicotinamide
mononucleotide, improved conventional measures without detectably restoring coordination, which is
consistent with that reading but underpowered to establish it.

Any claim stronger than this requires the replication cohort described in the Discussion: paired
insulin-stimulated biopsies in a modern sequenced cohort, which exist but are not deposited.
