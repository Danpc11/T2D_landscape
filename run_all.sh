#!/usr/bin/env bash
# Orquestador. Uso: bash run_all.sh [--workers=N] [--quick]
# R  : descarga/QC, redes de coexpresion (WGCNA), metricas, bootstrap, drivers, limma.
# Py : coherencia cross-tejido (haz), paisaje cuasi-potencial (U, J), simulaciones.
set -euo pipefail
WORKERS=32; QUICK=0
for a in "$@"; do case $a in --workers=*) WORKERS="${a#*=}";; --quick) QUICK=1;; esac; done
mkdir -p logs
run() { echo "=== $*"; "$@" > "logs/$(basename "$2" .R).log" 2>&1 || { echo "FALLO: $*  (ver logs/)"; exit 1; }; }
if [ $QUICK -eq 1 ]; then
  R_FLAGS_02="--p_sub=300 --n_sub_reps=3"; R_FLAGS_03="--workers=$WORKERS --n_boot=5 --n_perm=20"
  PY_SHEAF="--n_genes 300 --reps 5 --B 50 --no_loto --no_sens"; PY_LAND="--n_genes 300 --reps 5 --B 40"; DS="GSE27951"
else
  R_FLAGS_02=""; R_FLAGS_03="--workers=$WORKERS"
  PY_SHEAF="--n_genes 800 --reps 30 --B 500"; PY_LAND="--n_genes 800 --reps 20 --B 200"; DS=""
fi
[ -d data/processed ] || run Rscript R_scripts_preprocessing/01_download_qc_preprocess.R
run Rscript R_scripts_sub_networks/02_networks_modularity_metrics.R $R_FLAGS_02 $DS
run Rscript R_scripts_sub_networks/03_bootstrap_nulls_reference.R $R_FLAGS_03 $DS
run Rscript R_scripts_sub_networks/04_gene_drivers_and_enrichment.R --workers=$WORKERS $DS
if [ $QUICK -eq 0 ]; then
  run Rscript R_scripts_full_networks/02b_full_networks_biology.R --workers=$WORKERS
  run Rscript R_scripts_full_networks/03b_full_bootstrap_and_nulls.R --workers=$WORKERS
  run Rscript R_scripts_full_networks/04b_full_gene_drivers_and_enrichment.R --workers=$WORKERS
  run Rscript R_scripts_full_networks/04c_expression_differential.R
fi
run Rscript R_scripts_sub_networks/05_cross_dataset_summary.R
echo "=== python: coherencia (haz)"
python python/sheaf_coherence.py --export_dir export_sheaf --out results/sheaf --beta 6 $PY_SHEAF > logs/sheaf.log 2>&1
echo "=== python: paisaje (U, J)"
python python/landscape.py --export_dir export_sheaf --out results/landscape $PY_LAND > logs/landscape.log 2>&1
if [ $QUICK -eq 0 ]; then
  for b in 4 8; do python python/sheaf_coherence.py --export_dir export_sheaf --out results/sheaf_beta$b --beta $b --n_genes 800 --reps 20 --B 200 --no_loto --no_sens > logs/sheaf_b$b.log 2>&1; done
  python python/landscape.py --export_dir export_sheaf --out results/landscape_nocovar --no_covar $PY_LAND > logs/landscape_nocovar.log 2>&1
fi
echo "Listo. Resumen: results/summary/, results/sheaf/main_*.tsv, results/landscape/landscape_summary.tsv"
