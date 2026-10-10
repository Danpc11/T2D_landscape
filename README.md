# t2d_landscape

Code and documents for the paper *A coordinated transcriptional response to insulin is a tissue-level property of human skeletal muscle and is lost in insulin resistance* — a reanalysis of 26 public human studies (32 accessions, 1,717 participants) across islet, skeletal muscle, adipose tissue and blood, asking whether the information about insulin sensitivity is in how muscle *is* or in how it *responds*.

## What the paper shows

1. In resting muscle from 49 clamp-phenotyped individuals, no gene is associated with insulin sensitivity at FDR < 0.1. The 739 genes obtained without covariate adjustment, reproducing the published analysis of that cohort, disappear when expression is adjusted for **either age or BMI** alone.
2. Healthy muscle responds to insulin along a shared direction established over ~4 h, which includes late modulation of clock-output genes (DBP, PER2, NR1D2; pre-specified set of nine).
3. Insulin-resistant and diabetic muscle show no detectable reduction in magnitude but a loss of directional concentration: per-person alignment falls by 0.49 (*P* = 0.006 within hybridisation batch) and the von Mises–Fisher concentration κ falls from 8.2 to 3.6.
4. The loss is specific to insulin. In a hierarchical model over 140 participants and four cohorts, impairment reduces κ by 44% under insulin and not at all under acute exercise (interaction *P* = 0.005).
5. **The coordinated response belongs to the tissue, not the myocyte.** Fibre type does not account for it (adjusting leaves the effect at 0.40, *P* = 0.002) but mononuclear composition does (0.46 → 0.08, *P* = 0.49; random covariates leave it at 0.44). And in primary myotubes from 24 donors given 100 nM insulin, coherence never exceeds 0.27 against 0.77 in intact muscle, and does not differ between healthy and diabetic donors at any time point.
6. Tissues of one person share transcriptional architecture and position, but this is general rather than metabolic: all 30 GTEx tissue pairs show it, metabolic pairs more strongly. No cohort shows two stable states and there is no monotonic change with glycaemic stage.

## Repository

```
README.md               this file
THEORY.md               framework and formal definitions behind the measures (sheaf coherence, landscape, response geometry)
CHANGES.md              changelog
docs/
  MANUSCRIPT.md             manuscript draft (Nature Metabolism, Analysis)
  DRAFT_CellMetab_format.md earlier draft in Cell Press format (Highlights, eTOC, STAR Methods)
  FIGURE_PLAN.md            what each figure claims and which panels prove it
  FIGURE_LEGENDS.md         legend for every figure, panel by panel
  AUDIT.md                  every claim re-tested: survives, weakened or withdrawn
  REPLICATION.md            every cohort analysed and its verdict
  TARGETS.md                gene-level exploration and candidate targets
  REVIEW_RESPONSE.md        status of every point raised in the editorial review
  ROADMAP.md, THEORY_RESULTS.md, DRAFT_STORY.md
run_all.sh              discovery pipeline (R networks + python/sheaf_coherence.py + python/landscape.py)
run_replication.sh      every other analysis in the paper, in order (python/analyses/01–18)
requirements.txt        pinned Python versions used to produce the published numbers
                        (scanpy, igraph and leidenalg are needed only by step 19a)
renv_packages.txt       R packages used by run_all.sh
tests/smoke_test.py     checks the core statistics behave as claimed; needs no data
R_scripts_*/            coexpression networks, metrics, bootstrap, drivers, limma
python/sheaf_coherence.py, python/landscape.py
python/lib/             geo.py (GEO readers, platform annotation), response.py (response geometry, LOO coherence, permutation tests)
python/figures/         make_figures.py (Fig 1–6 from results/; T2D_PANEL_TITLES=1 draws
                        panel titles for internal review, off by default)
python/analyses/        01 cohorts · 02 landscape replication · 03a-c GTEx · 04 insulin response ·
                        05 resting muscle vs M (partial correlation) · 06 myotubes · 07 supplementary ·
                        08 audit · 09 classical differential expression · 10 network panel ·
                        11 exercise specificity · 12 confounding and power · 13 cohort inventory ·
                        14 sheaf across organ sets (specificity control) · 15 discovery sheaf by stage ·
                        16 per-person alignment (primary analysis) · 17 von Mises-Fisher model ·
                        18 hierarchical vMF across cohorts · 19a cell-type signature from scRNA-seq ·
                        19 tissue versus cell (composition and myotubes)
python/simulation/      synthetic validations and power
data/README.md          every input file and where to download it
```

## Reproduce

```bash
# 1. put the public inputs in data/raw/ as listed in data/README.md
# 2. discovery (needs R with WGCNA, limma, GEOquery; ~hours)
bash run_all.sh --workers=32
# 3. everything else (python 3.10+; pip install -r requirements.txt; ~45 min)
bash run_replication.sh          # add --quick for a fast, lower-resolution pass
python tests/smoke_test.py       # sanity check of the core statistics, needs no data
```

Each step writes to its own `results/<subfolder>`: `replication/`, `gtex/`, `discovery/`, `response/`,
`resting/`, `myotubes/`, `supplementary/`, `audit/`, `de/`, `networks/`, `power/`, `inventory/`.
Do not set `T2D_OUT` globally; the orchestrator unsets it so that each step uses its own default.
Then build the figures from those tables only:

```bash
python python/figures/make_figures.py                       # figures/Fig1..Fig6, 250 mm wide (drafting)
T2D_FIG_WIDTH_MM=180 python python/figures/make_figures.py  # Nature double-column width (submission)
```

`data/raw/`, `data/export/`, `results/` and `figures/` are not versioned.

Every script reads its inputs and writes its outputs through environment variables, so the whole
pipeline can run against directories elsewhere without editing code:
`T2D_RAW` (downloaded inputs, default `data/raw`), `T2D_EXPORT` (cohorts prepared by step 01,
default `data/export`), `T2D_RES` / `T2D_OUT` (results, default `results/...`) and `T2D_FIG`
(figures, default `figures`).

## Status

Manuscript in `docs/MANUSCRIPT.md` (Article format), supplementary notes in `docs/SUPPLEMENTARY.md`,
figure legends in `docs/FIGURE_LEGENDS.md`.

What the central claim rests on, stated plainly. The insulin contrast comes from one cohort, GSE22309
(2007, n = 20/20/15), in which hybridisation batch is unevenly distributed; the primary analysis is
therefore restricted to the 35 individuals whose two biopsies were processed in the same run, where
the effect is undiminished and the batch-only control is null. The specificity to insulin is estimated
jointly over 140 participants from four cohorts. Components of the result replicate in GSE9105,
GSE157988, GSE182117 and GSE182120. A modern sequenced clamp cohort with paired biopsies would settle
the main contrast and does not currently exist in public repositories: the deposited samples of the
largest recent clamp study are fasted only.

`docs/AUDIT.md` lists what survives, what is weakened and what was withdrawn; `docs/REVIEW_RESPONSE.md`
tracks the points raised in editorial review; `docs/EXPLORATION_BIOLOGY.md` records the analyses that
were tried and did not hold, so they are not repeated.

## Data and licence

All inputs are public (GEO, GTEx portal, Rydén et al. 2016 export). Code: MIT. Documents: CC BY 4.0.
