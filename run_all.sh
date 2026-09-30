#!/usr/bin/env bash
# Orchestrator. Usage: bash run_all.sh [--workers=N] [--quick]
# R  : download/QC, coexpression networks (WGCNA), metrics, bootstrap, drivers, limma.
# Py : cross-tissue coherence (sheaf), quasi-potential landscape (U, J).
set -uo pipefail
WORKERS=32; QUICK=0; FROM=""
for a in "$@"; do case $a in --workers=*) WORKERS="${a#*=}";; --quick) QUICK=1;; --from=*) FROM="${a#*=}";; esac; done
# --from=PASO reanuda desde ese paso (02, 03, 04, 02b, 03b, 04b, 04c, 05, sheaf, landscape)
STEPS="02 03 04 02b 03b 04b 04c 05 sheaf landscape sheaf_b4 sheaf_b8 landscape_nocovar landscape_nobalance"
SKIP=1; [ -z "$FROM" ] && SKIP=0
# Seguridad: no permitir dos corridas simultaneas (escriben en los mismos results/ y logs/)
LOCK=logs/run_all.lock; mkdir -p logs
if [ -e "$LOCK" ] && kill -0 "$(cat "$LOCK")" 2>/dev/null; then echo "Ya hay una corrida de run_all.sh en marcha (PID $(cat "$LOCK")). Abortando."; exit 1; fi
echo $$ > "$LOCK"; trap 'rm -f "$LOCK"' EXIT
mkdir -p logs
run() {  # run <logname> <cmd...>
  local name="$1"; local log="logs/$1.log"; shift
  if [ $SKIP -eq 1 ]; then [ "$name" = "$FROM" ] && SKIP=0 || { echo "--- omitido (--from=$FROM): $name"; return 0; }; fi
  echo "=== $*"
  if ! "$@" > "$log" 2>&1; then echo "FAILED: $*"; echo "--- last lines of $log:"; tail -15 "$log"; exit 1; fi
}
if [ $QUICK -eq 1 ]; then
  R02="--p_sub=300 --n_sub_reps=3"; R03="--workers=$WORKERS --n_boot=5 --n_perm=20"
  PY_SHEAF="--n_genes 300 --reps 5 --B 50 --no_loto --no_sens"; PY_LAND="--n_genes 300 --reps 5 --B 40"; DS="GSE27951"
else
  R02=""; R03="--workers=$WORKERS"
  PY_SHEAF="--n_genes 800 --reps 30 --B 500"; PY_LAND="--n_genes 800 --reps 20 --B 200"; DS=""
fi

# --- 01: download + QC (only if processed data are missing) ---
[ -d data/processed ] || [ -n "$FROM" ] || run 01 Rscript R_scripts_preprocessing/01_download_qc_preprocess.R

# --- export_sheaf/: regenerate from data/processed if missing (no download) ---
if [ ! -f export_sheaf/high_variance_genes_ordered.tsv ]; then
  echo "=== regenerating export_sheaf/ from data/processed/"
  Rscript -e '
    suppressPackageStartupMessages(library(data.table)); dir.create("export_sheaf", showWarnings = FALSE)
    hv <- readRDS("data/processed/high_variance_genes.rds")
    for (f in list.files("data/processed", pattern = "_processed\\.rds$", full.names = TRUE)) {
      acc <- sub("_processed\\.rds$", "", basename(f)); o <- readRDS(f)
      e <- o$expr[intersect(hv, rownames(o$expr)), , drop = FALSE]
      fwrite(data.table(gene = rownames(e), e), file.path("export_sheaf", paste0(acc, "_expr.tsv")), sep = "\t")
      fwrite(o$pheno, file.path("export_sheaf", paste0(acc, "_pheno.tsv")), sep = "\t")
    }
    fwrite(data.frame(gene = hv), "export_sheaf/high_variance_genes_ordered.tsv", sep = "\t")
    cat("export_sheaf/ ready\n")' > logs/export_sheaf.log 2>&1 || { echo "FAILED: export_sheaf regeneration"; tail -15 logs/export_sheaf.log; exit 1; }
fi

run 02  Rscript R_scripts_sub_networks/02_networks_modularity_metrics.R $R02 $DS
run 03  Rscript R_scripts_sub_networks/03_bootstrap_nulls_reference.R $R03 $DS
run 04  Rscript R_scripts_sub_networks/04_gene_drivers_and_enrichment.R --workers=$WORKERS $DS
if [ $QUICK -eq 0 ]; then
  run 02b Rscript R_scripts_full_networks/02b_full_networks_biology.R --workers=$WORKERS
  run 03b Rscript R_scripts_full_networks/03b_full_bootstrap_and_nulls.R --workers=$WORKERS
  run 04b Rscript R_scripts_full_networks/04b_full_gene_drivers_and_enrichment.R --workers=$WORKERS
  run 04c Rscript R_scripts_full_networks/04c_expression_differential.R
fi
run 05  Rscript R_scripts_sub_networks/05_cross_dataset_summary.R
run sheaf     python python/sheaf_coherence.py --export_dir export_sheaf --out results/sheaf --beta 6 $PY_SHEAF
run landscape python python/landscape.py --export_dir export_sheaf --out results/landscape $PY_LAND
if [ $QUICK -eq 0 ]; then
  for b in 4 8; do run sheaf_b$b python python/sheaf_coherence.py --export_dir export_sheaf --out results/sheaf_beta$b --beta $b --n_genes 800 --reps 20 --B 200 --no_loto --no_sens; done
  run landscape_nocovar   python python/landscape.py --export_dir export_sheaf --out results/landscape_nocovar --no_covar $PY_LAND
  run landscape_nobalance python python/landscape.py --export_dir export_sheaf --out results/landscape_nobalance --no_balance $PY_LAND
fi
echo "Done. Summaries: results/summary/, results/sheaf/main_*.tsv, results/landscape/landscape_summary.tsv"
