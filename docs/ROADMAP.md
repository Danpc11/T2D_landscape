# Roadmap: replication, the remission test, and the path to a clinical tool

## Part 1 — Replication cohorts (test F1, F2, F3, F5 on data that did not generate them)

Priority order. "Control" = continuous glycaemic variable available per subject.
All are bulk transcriptomes of the same tissues as the discovery cohorts.

| Tissue | Accession | Platform | n (stages) | Control | Why |
|---|---|---|---|---|---|
| **Islet** | **GSE164416** | RNA-seq (GPL16791) | 133 living donors: ND 18, IGT 41, T3cD 35, T2D 39 | HbA1c (in the paper's supplementary tables; not deposited in GEO) | Largest IGT group available; living donors metabolically phenotyped; decisive test of P2/F5 in islet |
| **Islet** | **GSE50244** | RNA-seq | 89 organ donors | HbA1c in GEO characteristics | Independent second replicate; continuous control for Fisher |
| Islet | GSE38642, GSE25724 | Affymetrix | 63 (9 T2D); 13 (6 T2D) | none | Small; use only for F2 sign |
| **Muscle** | GSE25462 | HG-U133 Plus 2 | 50: 25 ND (15 with family history), 25 T2D | fasting glucose, insulin | Independent of GSE18732; same platform family |
| Muscle | GSE22309 | HG-U133A | 110 (basal + insulin clamp): ND, IR, T2D | clamp M-value | Three stages; insulin-stimulated arm is a bonus (does the response contract?) |
| **Adipose SAT** | GSE64567 | Illumina | 70 (NGT/IGT/T2D) | HbA1c, glucose | Direct replicate of GSE27951 design |
| Adipose | GSE20950 | HG-U133 Plus 2 | 39 (SAT + omental), insulin sensitive vs resistant | HOMA-IR | Two depots per subject: within-subject coherence between SAT and VAT |
| Liver | GSE23343 | HG-U133 Plus 2 | 17 (10 T2D) | none | Only to replace the obesity-stratified GSE15653; still underpowered |

**Within-donor multi-tissue (the real test of F1).** GTEx v8/v10 has skeletal
muscle, subcutaneous and visceral adipose, liver and whole pancreas from the
*same* donors (hundreds), with a type 2 diabetes field (MHT2D) in the
protected phenotypes (dbGaP phs000424). This is the only resource where the
sheaf can be built on paired tissues rather than on independent cohorts, and
where F1 becomes a within-person statement. Requires a dbGaP application
(≈ 6–10 weeks); expression is open, the diabetes label is not. Start the
application now; it is the single biggest upgrade available.

**What replication must reproduce, in order of importance.**
1. F2 sign and magnitude: dispersion and entropy decrease healthy → T2D on the
   principal subspace (pre-specify r = 3, 800 HV genes, same covariates).
2. F5: no stable second basin; Cramér's V between GMM component and stage
   < 0.2; I_c monotone (no IGT peak).
3. F3: I_c(T2D) > I_c(healthy).
4. Fisher peak inside the prediabetic range when a continuous control exists.
5. F1 can only be replicated with paired tissues (GTEx) or by rebuilding the
   sheaf on independent replicate cohorts (islet GSE164416 + muscle GSE25462 +
   adipose GSE64567) and recovering z ≪ 0.

**Pipeline changes needed.** `01` hard-codes the four accessions and their
condition regexes. Add a `cohorts.tsv` (accession, tissue, platform, stage
column, regex per stage, control column, covariate columns) read by `01`,
`sheaf_coherence.py` (`STATE_ORDERS`) and `landscape.py`. RNA-seq cohorts
(GSE164416, GSE50244) need a count → log-CPM/vst path before the common
pipeline (voom or DESeq2 vst, then the same gene filtering).

## Part 2 — The remission test (prediction: re-expansion and re-coupling)

Paired biopsies before and after bariatric surgery, same subjects.

| Accession | Tissue | Design | n | Why |
|---|---|---|---|---|
| **GSE59034** | SAT | before / 2 years after RYGB (+ 5 years in the long-term series) | 16 paired (+ follow-up) | Longest follow-up; the "does the state re-expand and stay expanded" question |
| **GSE66921 / GSE66306** | SAT and peripheral monocytes, RNA-seq | before / 3 months after | 22 women (SAT), 19 (monocytes) | Early time point; monocytes give a *blood* readout of the same subjects (bridge to Part 3) |
| GSE29409 / GSE29410 | SAT + omental | before / short term after | 5 + 3 | Small; two depots per subject |
| GSE72158 | SAT | before / after | ~20 | Third replicate |
| DiaBar cohort (Berlin; PMC12309261) | SAT | baseline with 12-month remission outcome | ~100 | Remission vs non-remission at baseline: does baseline dispersion/depth predict who remits? (contact authors; not all in GEO) |

**Analysis (pre-specified).** For each paired dataset: (i) dispersion and
entropy of the before and after groups in a common embedding, paired
bootstrap of the difference; prediction: after > before. (ii) I_c before vs
after; prediction: after < before. (iii) Within-subject displacement vector
before → after; its projection on the healthy → T2D axis learned in the
discovery adipose cohort (GSE27951), batch-corrected; prediction: moves
towards healthy. (iv) In the 2-/5-year series: is the re-expansion sustained,
and does it track weight regain? (v) In cohorts with remission outcome:
baseline depth (−ln p̂ under the healthy reference) as a predictor of
remission, against DiaRem/DiaBar.

A positive result here is what turns the theory from descriptive into
dynamic: the same subject moving back up the funnel.

## Part 3 — From framework to clinical tool

### 3.1 What the tool would measure

Two per-patient quantities, computed from one tissue sample against a fixed
healthy reference distribution built from the replication cohorts:

- **Depth** d(x) = −ln p̂_healthy(x): how far the patient's transcriptional
  state is from the healthy region (Mahalanobis-type, on the organised
  subspace). Interpretable as "position in the funnel".
- **Fragmentation** f(x): the patient's residual (off-principal-axes)
  distance to the healthy manifold, the per-patient analogue of I_c.

And, when two or more tissues from the same patient are available (surgery,
research settings): **coherence** c(x) from the sheaf, the per-patient
analogue of F1.

What it is *not*: a diagnostic of T2D (glycaemia does that for free). Its
plausible clinical use is earlier and different: (a) staging within
prediabetes, where the transcriptome changes fastest and glycaemia is
uninformative (Fisher peak at HbA1c ≈ 6.3 %, glucose 5.5–5.8); (b) predicting
remission or relapse under an intervention; (c) a surrogate endpoint for
trials: does a drug re-expand and re-couple, or only move glycaemia.

### 3.2 The obstacle, and the way through it

Tissue biopsies are not a screening tool. The clinical version must work on
**blood**. Nothing in our data says the contraction is visible in blood; it
has to be tested. Resources, in order:

1. GSE66306 (monocytes before/after surgery, paired with SAT in GSE66921):
   does blood show the same re-expansion as adipose in the same people?
2. Whole-blood T2D cohorts with prediabetes and controls (GSE15932, GSE21321,
   GSE156993, GSE9006): is depth monotone with HbA1c in blood?
3. Population cohorts with baseline blood transcriptomes and incident T2D at
   follow-up (KORA F4/FF4, Framingham Offspring RNA-seq, MESA): the
   prospective test — does baseline depth predict conversion beyond HbA1c,
   BMI and age? This is the study that would justify clinical development,
   and it requires data-access applications (KORA, dbGaP for FHS/MESA).

If blood does not carry the signal, the tool stays a research/surgical
instrument (adipose biopsy at surgery, where DiaBar already shows SAT
expression predicts remission) and the framework remains a trial endpoint.

### 3.3 Development stages and what each needs

| Stage | Deliverable | Evidence required | Time |
|---|---|---|---|
| 0 (now) | Framework + discovery results | done | — |
| 1 | Replication in independent cohorts (Part 1) | F2, F5, F3 reproduced in ≥ 2 tissues | 2–3 months |
| 2 | Remission test (Part 2) | re-expansion after surgery; depth predicts remission | 2–3 months, parallel |
| 3 | Reference distribution + depth/fragmentation score, versioned, with uncertainty | calibration across platforms (RNA-seq vs arrays); batch robustness | 2 months after 1 |
| 4 | Blood feasibility | monotone depth in blood cohorts | 3 months |
| 5 | Prospective validation | incident-T2D prediction beyond clinical variables, external cohort, pre-registered | 12–18 months incl. access |
| 6 | Clinical-grade software | locked model, documented intended use, software-as-a-medical-device quality system (IEC 62304, ISO 13485) if diagnostic claims are made; otherwise research-use-only | after 5 |

Stages 1–3 are what we can do with public data and the current code.
Stage 5 is where a clinician and a biostatistician have to be on the team
and where ethics approvals and data-access agreements take the time, not
the computation.

### 3.4 The packages (Bioconductor / Python) in this plan

They belong at stage 3: a `reference` object (healthy distribution per
tissue and platform, with version and provenance), `depth()`,
`fragmentation()`, `coherence()` with bootstrap uncertainty, and the
landscape/sheaf functions already written. Building them before stage 1 is
premature; building them after stage 3 without the reference object is
pointless.

## Immediate actions

1. Download GSE164416 counts + supplementary clinical table; GSE50244;
   GSE25462; GSE64567; GSE59034; GSE66921/GSE66306.
2. Add `cohorts.tsv` support to `01` and the Python modules; add an RNA-seq
   ingestion path (vst).
3. Run F2/F5/F3 on each replication cohort with the pre-specified settings.
4. Run the paired before/after analysis on GSE59034 and GSE66921.
5. Submit the GTEx dbGaP application (paired multi-tissue, MHT2D).
6. In parallel, write the paper on discovery + replication + remission test;
   the clinical tool is the Discussion's last section, not a claim.
