# Gene-level exploration and candidate targets

Every gene list below is derived from one of the four replicated findings (R1–R4). Evidence level is stated for each. None of this is causal; these are the genes that *carry* the systemic phenomena we measured, and therefore the first places to look.

## A. Genes that carry the shared individual state (R1, GTEx, 253 paired donors)

Genes whose residual expression (after age, sex, death type, RIN, ischaemic time) tracks the canonical variate shared by muscle and adipose tissue, with the same sign in muscle, adipose and pancreas (|r| ≥ 0.42 in all three):

| Module | Genes | Known biology |
|---|---|---|
| Glucocorticoid / stress response | **FKBP5**, MT1X, MT2A, SERPINA3, **XBP1** | FKBP5 is the glucocorticoid-receptor co-chaperone; MT/SERPINA3 are glucocorticoid- and IL-6-induced; XBP1 is the ER-stress transcription factor |
| Inflammatory signalling | **IL1RL1** (ST2, IL-33 receptor), **STAT3**, IFNGR1, F2RL3 | IL-33/ST2 governs adipose regulatory T cells; STAT3 is the IL-6 effector |
| Insulin-signalling brake | **PTPN1** (PTP1B) | Dephosphorylates the insulin receptor; a long-pursued T2D target |
| Ribosome biogenesis / growth | DDX21, NCL, WDR43, RRP1, RRP12, PPAN, DIMT1, RCL1, EIF4A1 | mTORC1 output |
| Lipid droplet | PLIN2 | ectopic lipid storage |

Reading: the state that a person carries into every organ is, to first approximation, a **systemic stress–inflammation–glucocorticoid axis coupled to growth signalling**. Caveat: GTEx donors are post mortem, so part of this axis is terminal illness; it is nonetheless the same axis that clinical physiology implicates in insulin resistance, and the one to test first in living paired-tissue cohorts (FUSION).

## B. Genes whose network position becomes tissue-specific in T2D (R1 decay, discovery sheaf)

Highest per-gene incoherence in T2D across islet, muscle and adipose: **C1QA, C1QB, C3AR1, ITGB2, LSP1, IL1R1, HLA-DOA, S100A4, HPGDS, STK17B**. Complement and innate immunity reorganise differently in each organ. Not FDR-significant at n = 33–118; replicated in direction only.

## C. The insulin-response programme that survives in prediabetes, and the one that does not (R4)

**C1. Still coordinated in prediabetic muscle** (GSE157988, 23 women, clamp, |t| > 5 across subjects):
- Induced: PPP1R3B, PIK3R1, HES1, SPRY4, HMOX1, G0S2, MYOD1, IL6R, RRAD, NFIL3, MT1E, DNAJB1
- Repressed: **TXNIP, KLF15, IRS2, SLC27A1**, and the entire clock output **PER1, PER2, PER3, NR1D1, NR1D2, DBP, TEF, HLF, CIART**

This is the canonical insulin transcriptional programme (TXNIP/KLF15 repression, PPP1R3B/HES1 induction). It is intact. The loss of coordination we measure at the whole-transcriptome level is therefore **not** in the canonical programme but in the non-canonical, individual part of the response. That has a direct consequence: targeting canonical insulin signalling would not be expected to restore coordination, which is what the NMN and bariatric results showed.

**C2. Healthy meal response lost in T2D muscle** (GSE231509, 1 h after a mixed meal, |t| > 4 in healthy, |t| < 1 in T2D): **PFKFB3** (glycolytic flux regulator), **CEBPD**, DNAJB4, ANKRD37, PHC2, LRIG1, GNA13, JMJD1C, RBM20, DMD, BMPR1B. Preserved in T2D: ANGPTL4, DEPP1, CCN1/CCN2, DBP, CIART, DUSP16. The genes that stop responding are regulators of glycolytic flux, chaperoning and chromatin rather than core insulin signalling; n = 7 per group, hypothesis level.

**C3. The clock–insulin coupling.** Insulin acutely represses the muscle clock output genes in every subject, prediabetic included (C1), and the diurnal rhythm of adipose gene expression is reduced in T2D (GSE104674, the authors' own finding). The coordinated response we measure is built over hours (0.29 at 30 min → 0.77 at 4 h); a time-structured, clock-coupled programme is what would be expected to fragment across individuals when the clock is desynchronised.

## D. Islet genes already changed in the prediabetic window (R3, GSE50244 vs GSE164416)

Few and weakly replicated at this n: **SST** (somatostatin, δ-cell) down and **CAPN2** up in prediabetes in both cohorts; F3, CD9, NEDD9, EPHB1 up in one. The prediabetic window is visible as a distribution-level phenomenon (Fisher peak) before it is visible gene by gene.

## E. Candidate targets, ranked by how much of our evidence converges on them

| Target | Our evidence | Existing therapeutic status | What our framework adds |
|---|---|---|---|
| **IL-1 axis (IL1R1, IL1RL1/ST2, IL-33)** | B (tissue-specific reorganisation) + A (systemic axis) | anakinra improved glycaemia in a small trial; canakinumab (CANTOS) lowered inflammation without preventing diabetes | Reorganisation is tissue-specific: a systemic blocker acts on one axis in organs that have decoupled; predicts that re-coupling, not CRP, is the endpoint to measure |
| **Glucocorticoid signalling (FKBP5)** | A, strongest systemic module | FKBP51 antagonists (SAFit2) improve insulin signalling and adiposity in mice; no human trial | The systemic individual state is largely a glucocorticoid-stress axis; FKBP5 is its most druggable node |
| **PTP1B (PTPN1)** | A | multiple inhibitors failed on pharmacology, not biology; trodusquemine in trials | PTPN1 tracks the whole-organism state, not one tissue; supports systemic rather than organ-targeted dosing |
| **Complement C3a/C3aR** | B | C3aR antagonists preclinical; C3 associates with insulin resistance | Which organ decouples first for the complement module is measurable with the sheaf and would order the intervention |
| **TXNIP** | C1 (most coordinated insulin-repressed gene, intact in prediabetes) | verapamil represses TXNIP; trials in type 1 and early type 2 diabetes | Intact response means TXNIP is a marker of preserved canonical signalling, not of the coordination defect |
| **Clock (REV-ERB / NR1D1-2, PER)** | C1 + C3 | REV-ERB agonists preclinical; chronotherapy | The coordination defect is time-structured; timing of insulin/meal exposure is a modifiable variable with no drug needed |
| **PFKFB3** | C2 | PFKFB3 inhibitors exist (oncology) | Loss of meal-induced PFKFB3 response in T2D muscle; glycolytic flux gating as the earliest lost response |
| **ER stress (XBP1)** | A | TUDCA/4-PBA improved insulin sensitivity in small trials | Part of the systemic axis |

## F. What this means for therapy, in one paragraph

The canonical insulin programme (TXNIP, KLF15, IRS2, PPP1R3B) still responds in prediabetic muscle; what is lost is the coordination of the rest of the response between individuals, and neither surgery nor NMN restored it while both restored magnitude or sensitivity. The genes that carry the systemic individual state are a glucocorticoid–inflammation–growth axis (FKBP5, IL1RL1, STAT3, PTPN1, XBP1), and the genes that decouple between organs in diabetes are complement and innate immunity. Together these point away from "more insulin signalling in one organ" and towards interventions on the systemic stress–immune axis, timed with the clock, with re-coordination of the response across individuals and organs as the endpoint. Each line of this is a hypothesis with the dataset that would test it named in `REPLICATION.md`; FUSION (living paired tissues with glycaemic staging) is the one that tests A and B together.
