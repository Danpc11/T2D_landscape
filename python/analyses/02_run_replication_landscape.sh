#!/usr/bin/env bash
# Corre landscape.py sobre cada cohorte de replicacion exportada por 01. Uso: bash 02_run_replication_landscape.sh [--quick]
set -euo pipefail; cd "$(dirname "$0")/../.."
Q=""; [ "${1:-}" = "--quick" ] && Q="--n_genes 300 --reps 5 --B 40"
for acc in GSE164416 GSE50244 GSE50398 GSE25462 METSIM; do
  python python/analyses/../landscape.py --export_dir data/export/$acc --out results/replication/$acc --tissues $acc --n_genes 800 --reps 15 --B 100 $Q || echo "FALLO $acc"
done
