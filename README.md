# T2D landscape

Code and documents for the paper *Insulin resistance is a loss of coordination, not of signal* — a reanalysis of 23 public human transcriptomic cohorts (islet, skeletal muscle, adipose tissue, blood; GTEx paired tissues; hyperinsulinaemic-clamp biopsies) asking whether the information about insulin sensitivity is in how muscle *is* or in how it *responds*.

## What the paper shows

1. The resting muscle transcriptome of 49 clamped individuals carries no gene-level information about insulin sensitivity (0 genes at FDR < 0.1 vs clamp M-value).
2. Healthy muscle responds to insulin along one direction shared by individuals, assembled over ~4 h, which includes a reset of the peripheral clock output (DBP, PER2, NR1D2).
3. Insulin-resistant and diabetic muscle respond with the same magnitude but without a shared direction, and the clock output no longer responds.
4. Adipose tissue fails by magnitude, not direction; bariatric surgery restores magnitude, not coordination.
5. The organs of one person share one transcriptional state (cellular-sheaf coherence) that drifts continuously with disease, without a second basin.

## Repository

```
README.md               this file
THEORY.md               framework and formal definitions behind the measures (sheaf coherence, landscape, response geometry)
CHANGES.md              changelog
docs/                   manuscript draft, audit of every claim, replication table, gene-level exploration, roadmap
run_all.sh              discovery pipeline (R networks + python/sheaf_coherence.py + python/landscape.py)
run_replication.sh      every other analysis in the paper, in order (python/analyses/01–08)
R_scripts_*/            coexpression networks, metrics, bootstrap, drivers, limma
python/sheaf_coherence.py, python/landscape.py
python/lib/             geo.py (GEO readers, platform annotation), response.py (response geometry, LOO coherence, permutation tests)
python/analyses/        01 cohorts · 02 landscape replication · 03 GTEx · 04 insulin response · 05 resting muscle vs M · 06 myotubes · 07 supplementary · 08 audit
python/simulation/      synthetic validations and power
data/README.md          every input file and where to download it
```

## Reproduce

```bash
# 1. put the public inputs in data/raw/ as listed in data/README.md
# 2. discovery (needs R with WGCNA, limma, GEOquery; ~hours)
bash run_all.sh --workers=32
# 3. everything else (python 3.10+, numpy, scipy, pandas, openpyxl; ~30 min)
bash run_replication.sh
```

Outputs land in `results/` (`response/`, `resting/`, `myotubes/`, `gtex/`, `replication/`, `supplementary/`, `audit/`). Figures and manuscript numbers are generated from these TSVs only.

## Status

Manuscript in preparation. The central claim (loss of coordination and of the insulin→clock coupling in insulin resistance) rests on GSE22309 (2007, n = 20/20/15, batch-uneven) with replication of its components in GSE9105, GSE157988, GSE182117 and GSE182120; see `docs/AUDIT.md` for what survives, what is weakened and what was withdrawn, and `docs/REPLICATION.md` for every cohort.

## Data and licence

All inputs are public (GEO, GTEx portal, Rydén et al. 2016 export). Code: MIT. Documents: CC BY 4.0.
