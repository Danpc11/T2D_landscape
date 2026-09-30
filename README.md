# T2D Landscape

**Attractor landscape of the multi-tissue transcriptome in type 2 diabetes**

Healthy metabolic homeostasis and type 2 diabetes (T2D) are treated as two attractors of the multi-organ transcriptional system, with impaired glucose tolerance (IGT) as the transition state between them. The pipeline asks whether this bistability — predicted by physiological models of T2D (Topp et al. 2000; Ha, Satin & Sherman 2016) — leaves a measurable footprint in the tissue transcriptomes of four metabolically coupled organs, and whether the transition is a single-tissue event or a loss of coordination between tissues. Full framework and pre-specified predictions P0–P5: `THEORY.md`.

## Central objects

1. **Quasi-potential landscape** U(x) = −ln P(x) per tissue on a low-dimensional embedding of the transcriptome, conditioned on covariates (`python/landscape.py`): persistent basins (Morse theory / sublevel-set persistence), barrier asymmetry, critical-transition index, Fisher–Rao information along HbA1c/glucose, and the non-equilibrium flux J from a Schrödinger bridge between healthy and T2D distributions (Hodge decomposition of the inferred drift).
2. **Cross-tissue coherence** via a cellular sheaf over the tissue graph (`python/sheaf_coherence.py`): stalks are aligned spectral embeddings of each tissue's coexpression network; the sheaf energy measures whether network reorganization is systemic or tissue-specific.
3. **Network geometry of each attractor** (`R_scripts_*`): global efficiency, communicability, effective resistance, strength entropy and modularity of weighted bicor networks, reported at equal sample size and relative to a strength-preserving configuration null.

### Network metrics

| Metric               | Symbol | Definition                                  | Reading                                   |
|----------------------|--------|---------------------------------------------|-------------------------------------------|
| Global efficiency    |  EG    | Mean over pairs of 1/d_ij                   | Short-path accessibility                  |
| Communicability      |  Ḡ     | Mean off-diagonal entry of exp(D⁻¹/²WD⁻¹/²) | All-path propagation                      |
| Effective resistance |  R̄     | 2·tr(L⁺)/(p−1)                              | Redundancy / robustness to edge removal   |
| Strength entropy     |  Hb    | −Σ pᵢ log pᵢ, pᵢ = kᵢ/Σkⱼ                   | Evenness of connectivity                  |

The **Network Organization Index** NOI = z(EG) + z(Ḡ) − z(R̄) + z(Hb) is computed on equal-n metrics; `NOI_rel` uses the metrics relative to the configuration null (structure, not density). These are descriptive of the attractors' geometry (prediction S in `THEORY.md`); the inferential predictions rest on the landscape and sheaf modules.

### Two-Branch Architecture

The pipeline runs two parallel branches that address complementary scientific questions:

**Shared branch (scripts 02–05):** Operates on a fixed subnetwork of up to 800 high-variance genes common to all four datasets. Enables direct cross-tissue comparison of network metrics and gene driver scores on identical nodes.

**Full branch (scripts 02b–04b):** Operates on each dataset's complete transcriptome (p = 1,851–21,755 genes). Captures tissue-specific network architecture at full resolution, at the cost of losing cross-dataset node-level comparability.

---

## Datasets: The Metabolic Quartet

| GEO accession | Tissue            | Technology               | n   | States                              |
|---------------|-------------------|--------------------------|-----|-------------------------------------|
| GSE76895      | Pancreatic islets | Affymetrix GPL570        | 83  |(32 ND, 15 IGT, 36 T2D)              |
| GSE18732      | Skeletal muscle   | Affymetrix (Ensembl CDF) | 118 |(47 ND, 26 IGT, 45 T2D)              |
| GSE15653.     | Liver             | Affymetrix GPL570        | 18  |(5 Lean, 4 Obese-noT2D, 9 Obese-T2D) |
| GSE27951      | Adipose tissue    | Affymetrix GPL570        | 33  |(NGT, IGT, T2D)                      |

T3cD samples in GSE76895 are excluded. NGT in GSE27951 is treated as the healthy reference equivalent to ND.

---

## Pipeline Structure

```
01_download_qc_preprocess.R
        │
        ├─ data/processed/{acc}_processed.rds       (shared branch input)
        ├─ data/processed/{acc}_processed_full.rds  (full branch input)
        └─ data/processed/high_variance_genes.rds   (shared subnetwork nodes)
        │
        ├── SHARED BRANCH ──────────────────────────────────────────────────────┐
        │                                                                       │
        ▼                                                                       │
02_networks_modularity_metrics.R                                    │
        │  Fixed subnetwork (≤800 genes)                                        │
        │  bicor adjacency + WGCNA soft threshold                               │
        │  EG, Ḡ, R̄, Hb, NOI (global z-score)                                   │
        │  Modules: hclust + cutree, maximize n_mod × Q                         │
        └─ results/networks/, results/metrics/                                  │
                │                                                               │
                ▼                                                               │
03_bootstrap_nulls_reference.R                                                    │
        │  Bootstrap n=100: IC for all 4 metrics                                │
        │  Permutation tests n=1000 (EG+Hb): p_min=0.001                        │
        │  Configuration reference: W* ∝ (kᵢkⱼ)^(1/α)                               │
        └─ results/bootstrap/, results/nulls/, results/reference/                 │
                │                                                               │
                ▼                                                               │
04_gene_drivers_and_enrichment.R                                                │
        │  Node metrics: strength, eigencentrality, participation               │
        │  PTI, IRI, TRI (exact per-gene recomputation)                       │
        │  KO-support; GO:BP + KEGG enrichment (top-50)                         │
        └─ results/gene_drivers/, results/enrichment/                           │
                │                                                               │
                ├── FULL BRANCH ──────────────────────────────────────────────┐ │
                │                                                             │ │
                ▼                                                             │ │
        02b_full_networks_biology_optimized.R                                 │ │
                │  Complete transcriptome (p = 1k–22k genes)                  │ │
                │  bicor ONE pass + scale-free β selection                    │ │
                │  EG: top-0.1% edges (proportional sparsification)           │ │
                │  Ḡ, R̄: top-3000 hub genes subnetwork                        │ │
                │  Modules: blockwiseModules (WGCNA standard)                 │ │
                └─ results/full_networks/, results/full_metrics/              │ │
                        │                                                     │ │
                        ▼                                                     │ │
                03b_full_bootstrap_and_nulls_optimized.R                      │ │
                        │  Bootstrap n=100 (top-2000 genes): 4 metrics        │ │
                        │  Permutations n=1000 (EG+Hb): p_min=0.001           │ │
                        │  Optimum: Ḡ, R̄ on pre-aligned top-3000 subnet       │ │
                        └─ results/full_bootstrap/, results/full_nulls/,      │ │
                           results/full_optimum/                              │ │
                                │                                             │ │
                                ▼                                             │ │
                        04b_full_gene_drivers_and_enrichment.R                │ │
                                │  PTI/IRI/TRI/KO on top-3000 genes           │ │
                                │  mclapply COW for large W matrices          │ │
                                │  Top-1% edge sparsification for metrics     │ │
                                └─ results/full_gene_drivers/,                │ │
                                   results/full_enrichment/                   │ │
                                                                              │ │
        ◄─────────────────────────────────────────────────────────────────────┘ │
        │                                                                       │
        ◄───────────────────────────────────────────────────────────────────────┘
        │
        ▼
04c_expression_differential.R
        │  limma (correct for RMA-normalized microarrays)
        │  Contrasts: all states vs healthy + IGT→T2D transition
        │  combined_score = rank_pct(PTI+IRI+TRI) × |logFC|
        └─ results/differential_expression/
                │
                ▼
05_cross_dataset_summary.R
        │  Consolidates both branches
        │  Trend tables: EG_change, R̄_change, Ḡ_change, Hb_change
        │  Cross-tissue recurrent driver genes
        └─ results/summary/
```

---

## Mathematical Invariants

These invariants are enforced identically across all scripts and both branches:

**Effective resistance:**
```
R̄ = 2·tr(L⁺)/(p−1)
```
Derivation: Kirchhoff index Kf = Σ_{i<j} R_ij = p·tr(L⁺); mean over the p(p−1)/2 pairs gives 2·tr(L⁺)/(p−1). (An earlier "correction" to `/p` was wrong; numerically the difference is (p−1)/p.)

**Global efficiency:** EG = mean over pairs of 1/d_ij (a spurious ×2 factor was removed).

**Cross-tissue coherence (cellular sheaf):** see `python/sheaf_coherence.py` and `CHANGES.md`. Stalks are aligned spectral embeddings of each tissue's network; restriction maps are O(r) Procrustes rotations; the sheaf energy per stage measures how tissue-specific the network reorganization is.

**Soft-thresholding power (β):** Selected once **per dataset** on the healthy reference state as the smallest β ∈ {1,…,20} with scale-free fit R² ≥ 0.80 (WGCNA convention); fallback β = 6. The same β is used for every state, bootstrap and permutation of that dataset, so that metrics compare biology rather than β.

**Configuration-model reference (formerly "configuration reference"):**
```
W*ᵢⱼ ∝ (kᵢ · kⱼ)^(1/α),  subject to Σᵢ<ⱼ (Wᵢⱼ)^α = C
```
where kᵢ are nodal strengths of the healthy reference network and C is its wiring budget. Up to the exponent this is the weighted configuration (Chung–Lu) model: the network with the same strength sequence and no structure. Deviation from it therefore measures the amount of structure (modularity), and the pipeline also reports every metric relative to this null (`*_rel`) to separate structure from density.

**NOI (global z-score):** Always computed across all datasets and states jointly, never within a single dataset. This ensures cross-dataset comparability of sign and magnitude.

---

## Key Methodological Decisions

**bicor vs. Pearson:** Biweight midcorrelation (bicor) is robust to outlier samples common in human tissue microarray data. `maxPOutliers = 0.1` is set consistently across all scripts.

**blockwiseModules in 02b, hclust in 02:** For p ≤ 800 (shared branch), hclust on `1−W` is exact and fast. For p > 3,000 (full branch), TOM-based blockwiseModules is the published WGCNA standard (Langfelder & Horvath 2008). Edge-sparsified fast_greedy was rejected because sparsification can fragment genuine modules.

**EG + Hb for permutation tests, not all 4 metrics:** Gbar (expm, O(p³)) and Rbar (eigen, O(p³)) called 1,000 × 2 × 2 × 4 = 16,000 times would require hours. EG and Hb together detect the same group-level differences under label permutation H0. Gbar and Rbar are reported from observed networks (02b), not the null distribution.

**Rank-percentile for combined_score:** Min-max scaling introduces an edge artifact where the gene with the lowest topological score always gets `combined_score = 0` regardless of |logFC|. Rank normalization to [0,1] is monotonic, has no edge artifacts, and maps naturally to "topological percentile × DE magnitude."

**TRI/KO:** Computed by exact recomputation per gene (EG+Hb in the shared branch, p ≤ 800; R̄+Hb in the full branch, p ≤ 1500). The earlier "rank-1 Sherman–Morrison update" was invalid (ΔL is not rank 1) and was removed. These scores are topological sensitivities of the correlation network to a node, not causal rescue/knock-out evidence.

**Gbar and Rbar on top-3000 hub genes (02b):** Both are computed with a full eigendecomposition on the top-3,000-strength subnetwork (seconds). Spectral truncation was removed: the spectrum of D⁻¹ᐟ²WD⁻¹ᐟ² lies in [−1, 1], so exp(λ) does not decay and truncating to k eigenvectors gives 30–50 % error.

---

## Computational Requirements

| Script |   RAM  | CPU (40 workers) | Bottleneck                      |
|--------|--------|------------------|---------------------------------|
| 01     | ~4 GB  | ~20 min          | GEO download                    |
| 02     | ~2 GB  | ~5 min           | bicor (p≤800)                   |
| 02b    | ~32 GB | ~30–60 min       | bicor (p=14k), blockwiseModules |
| 03     | ~4 GB  | ~10 min          | 100 bootstrap × 4 metrics       |
| 03b    | ~16 GB | ~20 min          | 1000 permutations (EG+Hb)       |
| 04     | ~8 GB  | ~5 min           | TRI rank-1 (p≤800)              |
| 04b    | ~32 GB | ~15 min          | TRI rank-1 (p≤1500)             |
| 04c    | ~4 GB  | ~5 min           | limma + enrichGO                |
| 05     | ~2 GB  | <1 min           | File consolidation              |

### Required R packages

```r
# CRAN
install.packages(c(
  "WGCNA", "data.table", "dplyr", "tidyr", "igraph",
  "expm", "Matrix", "parallel", "foreach", "doParallel",
  "stringr", "tibble"
))

# Recommended (graceful fallback if absent)
install.packages(c("RSpectra", "doRNG"))

# Bioconductor
BiocManager::install(c(
  "limma", "GEOquery", "clusterProfiler",
  "org.Hs.eg.db", "AnnotationDbi", "enrichplot",
  "hthgu133a.db"
))
```

---

## Execution Order

```bash
# 1. Download, QC, normalize (required first)
Rscript 01_download_qc_preprocess.R

# 2a. Shared branch (sequential)
Rscript 02_networks_modularity_metrics.R
Rscript 03_bootstrap_nulls_reference.R
Rscript 04_gene_drivers_and_enrichment.R

# 2b. Full branch (HPC; can run in parallel with 2a)
Rscript 02b_full_networks_biology_optimized.R
Rscript 03b_full_bootstrap_and_nulls_optimized.R
Rscript 04b_full_gene_drivers_and_enrichment.R

# 3. Differential expression + integration (requires both branches)
Rscript 04c_expression_differential.R

# 4. Cross-dataset summary (requires all previous steps)
Rscript 05_cross_dataset_summary.R
```

Individual datasets can be targeted: `Rscript 02_networks.R GSE76895 GSE18732`

---

## Output Files

```
results/
├── qc/                         # Sample counts, phenotype tables, gene lists
├── metrics/                    # Shared: EG, Ḡ, R̄, Hb, NOI per state
├── full_metrics/               # Full: same metrics at full resolution
├── networks/                   # Shared: adjacency matrices W, β
├── full_networks/              # Full: adjacency matrices W, β
├── modules/                    # Shared: module membership per gene
├── bootstrap/                  # Shared: bootstrap metric distributions
├── nulls/                      # Shared: permutation p-values (EG, Hb, NOI)
├── optimum/                    # Shared: W*, ΔE, ΔR, ΔG, ΔH
├── full_bootstrap/             # Full branch equivalents
├── full_nulls/
├── full_optimum/
├── gene_drivers/               # PTI, IRI, TRI, KO, class, GO/KEGG
├── full_gene_drivers/
├── enrichment/
├── full_enrichment/
├── differential_expression/    # limma results + combined_score integration
└── summary/                    # Cross-dataset consolidated tables
```

### Key output files

| File | Content |
|------|---------|
| `results/metrics/all_network_metrics.tsv`           | Main network metrics (shared branch) |
| `results/full_metrics/all_full_network_metrics.tsv` | Full-resolution metrics                  |
| `results/summary/cross_dataset_trends.tsv`              | EG_change, R̄_change per tissue/branch    |
| `results/summary/all_deviation_from_optimum.tsv`        | ΔE, ΔR, ΔG, ΔH relative to W*            |
| `results/gene_drivers/{acc}_gene_driver_scores.tsv`     | Per-gene PTI/IRI/TRI/class               |
| `results/differential_expression/{acc}_full_integrated_drivers_DE.tsv` | combined_score            |
| `results/summary/all_integrated_drivers_DE.tsv`         | Cross-tissue integrated ranking          |

---

## Gene Driver Score Definitions

| Score | Full name              | Definition                  | High value means |
|-------|------------------------|-----------------------------|------------------|
| PTI | Phase Transition Index   | Node topology change ND→IGT | Gene rewires most at the critical transition |
| IRI | Irreversibility Index    | Node topology change ND→T2D | Gene is most altered in established disease  |
| TRI | Topological Rescue Index | ΔEG + ΔHb if gene's edges restored to healthy | Restoring this gene improves network most.  |
| KO  | KO-support               | ΔΔEG if gene is knocked out | Gene sustains the diseased network topology  |
| combined_score | Integrated driver | rank_pct(PTI+IRI+TRI) × \|logFC\| | Topologically important AND differentially expressed |

Gene classes (mutually exclusive, evaluated in priority order):

| Class                    | Criteria            | Biological interpretation                   |
|--------------------------|---------------------|---------------------------------------------|
| `transition_rescue`      | high PTI + high TRI | Early driver with rescue potential          |
| `irreversibility_rescue` | high IRI + high TRI | Late driver with rescue potential           |
| `irreversibility_lock`   | high IRI + high KO  | Sustains the diseased attractor             |
| `transition_driver`      | high PTI            | Reorganization driver without rescue effect |
| `other`                  | none of the above   | Background topology                         |

---

## Reproducibility

- All parallel scripts: `set.seed(1234)` with `RNGkind("L'Ecuyer-CMRG")`
- doRNG (if installed) provides per-iteration reproducibility in foreach loops
- mclapply uses L'Ecuyer-CMRG streams automatically via `set.seed()` before forking
- β values are saved in network `.rds` files and reused identically in all downstream scripts
- All scripts accept dataset targets as command-line arguments for partial re-runs

---


## Architecture

The repository is split by what each language does best:

| Layer | Language | Why | Files |
|---|---|---|---|
| Download, QC, annotation, coexpression networks, WGCNA modules, bootstrap/permutation nulls, limma DE, enrichment | **R** | GEOquery / WGCNA / limma / clusterProfiler have no Python equivalent of the same maturity | `R_scripts_*/` |
| Cross-tissue coherence (cellular sheaf), quasi-potential landscape (persistence, score, Fisher, Schrödinger bridge, Hodge), power simulations | **Python** | numpy/scipy linear algebra, TDA and OT tooling | `python/` |
| Orchestration | bash | one entry point, `--quick` smoke test | `run_all.sh` |

Interface between layers: `01` writes `export_sheaf/` (expression TSV per tissue restricted to common high-variance genes + pheno); the Python modules read only that. Theory: `THEORY.md`. Change log: `CHANGES.md`.

```bash
bash run_all.sh --quick             # ~10 min, one tissue, reduced iterations
bash run_all.sh --workers=32        # full run, ~4 h
```

## Sample-size and design caveats

Metrics are also reported at equal n within each dataset (`*_sub`), since smaller groups yield noisier correlations and therefore denser networks. GSE15653 (n = 5/4/9) is excluded from the intermediate-state question. Power analysis for the sheaf coherence analysis with the real group sizes is in `python/simulation/power_curve.py`.

## Author 

Daniel Pérez Calixto

Instituto Nacional de Médicina Genómica

Contact info: dperez@inmegen.gob.mx
