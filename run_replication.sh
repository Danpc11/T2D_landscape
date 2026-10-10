#!/usr/bin/env bash
# Reproduce every analysis beyond the discovery pipeline (which is run_all.sh). Outputs under results/.
set -uo pipefail; cd "$(dirname "$0")"; export T2D_RAW=${T2D_RAW:-data/raw} T2D_EXPORT=${T2D_EXPORT:-data/export}
unset T2D_OUT   # cada paso escribe en su propio results/<subcarpeta>; fijarlo globalmente los mezclaria
mkdir -p logs results
step() { local name="$1"; shift; local script=""; for a in "$@"; do case "$a" in *.py|*.sh) script="$a";; esac; done
  local log="logs/$(basename "${script:-$name}" | sed 's/\.[a-z]*$//').log"; echo "=== $name"; "$@" > "$log" 2>&1 || echo "    FAILED (see $log)"; }
step "01 prepare replication cohorts"   python python/analyses/01_prepare_replication_cohorts.py
step "02 landscape on replication"      bash   python/analyses/02_run_replication_landscape.sh ${1:-}
step "03a GTEx prepare"                 python python/analyses/03a_gtex_prepare.py
step "03b GTEx paired coherence"        python python/analyses/03b_gtex_paired.py
step "03c GTEx systemic axis genes"     python python/analyses/03c_gtex_axis_genes.py
step "14 sheaf across organ sets"       python python/analyses/14_sheaf_organ_sets.py
step "15 discovery sheaf by stage"      python python/analyses/15_discovery_sheaf.py
step "04 insulin response geometry"     python python/analyses/04_insulin_response.py
step "05 resting muscle vs M"           python python/analyses/05_resting_muscle_vs_M.py
step "06 myotube clock"                 python python/analyses/06_myotube_clock.py
step "07 supplementary"                 python python/analyses/07_remission_and_blood.py
step "08 audit: Fisher null"            python python/analyses/08_audit_fisher_null.py
step "09 classical DE + enrichment"     python python/analyses/09_classical_de.py
step "10 network panel inputs"          python python/analyses/10_network_panel.py
step "11 exercise specificity"          python python/analyses/11_exercise_specificity.py
step "12 power and confounding"         python python/analyses/12_power_and_confounding.py
step "13 cohort inventory"              python python/analyses/13_cohort_inventory.py
step "16 alignment (primary)"           python python/analyses/16_alignment_primary.py
step "17 von Mises-Fisher model"        python python/analyses/17_vmf_model.py
step "18 hierarchical vMF (pooled)"     python python/analyses/18_hierarchical_vmf.py
step "figures"                          python python/figures/make_figures.py
echo "done: results/{response,resting,myotubes,gtex,replication,supplementary,audit,de} and figures/"
