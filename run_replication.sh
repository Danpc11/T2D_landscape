#!/usr/bin/env bash
# Reproduce every analysis beyond the discovery pipeline (which is run_all.sh). Outputs under results/.
set -uo pipefail; cd "$(dirname "$0")"; export T2D_RAW=${T2D_RAW:-data/raw}; mkdir -p logs results
step() { echo "=== $1"; shift; "$@" > "logs/$(basename "$1" .py).log" 2>&1 || { echo "FAILED (see logs/)"; }; }
step "01 prepare replication cohorts"   python python/analyses/01_prepare_replication_cohorts.py
step "02 landscape on replication"      bash   python/analyses/02_run_replication_landscape.sh ${1:-}
step "03a GTEx prepare"                 python python/analyses/03a_gtex_prepare.py
step "03b GTEx paired coherence"        python python/analyses/03b_gtex_paired.py
step "03c GTEx systemic axis genes"     python python/analyses/03c_gtex_axis_genes.py
step "04 insulin response geometry"     python python/analyses/04_insulin_response.py
step "05 resting muscle vs M"           python python/analyses/05_resting_muscle_vs_M.py
step "06 myotube clock"                 python python/analyses/06_myotube_clock.py
step "07 supplementary"                 python python/analyses/07_remission_and_blood.py
step "08 audit: Fisher null"            python python/analyses/08_audit_fisher_null.py
echo "done: results/response, results/resting, results/myotubes, results/gtex, results/replication, results/supplementary, results/audit"
