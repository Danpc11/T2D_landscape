#!/usr/bin/env python3
"""Analisis primario: alineamiento por persona.

La coherencia de grupo resume n individuos en un solo numero, de modo que una comparacion entre
dos grupos dispone de dos observaciones. El alineamiento por persona conserva un valor por
participante y compara n frente a n, usando el mismo objeto geometrico. Incluye el control
decisivo: pseudo-grupos definidos solo por lote de hibridacion.
Salida: results/response/alignment_per_person.tsv, alignment_tests.tsv
"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo, response as R
OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
rows, tests = [], []

def analyse(name, D, groups, reference, batches=None, subsets=("all",)):
    groups = np.asarray(groups)
    al = R.alignment_per_person(D, groups == reference)
    mag = np.linalg.norm(D, axis=1)
    for a, g, m in zip(al, groups, mag): rows.append(dict(dataset=name, group=g, alignment=a, magnitude=m))
    others = [g for g in pd.unique(groups) if g != reference]
    for sub in subsets:
        mask = np.ones(len(D), bool) if sub == "all" else np.asarray(batches, bool)
        for o in list(others) + (["+".join(map(str, others))] if len(others) > 1 else []):
            sel = mask & np.isin(groups, [reference] + o.split("+"))
            if (groups[sel] == reference).sum() < 4 or (groups[sel] != reference).sum() < 4: continue
            lb = np.where(groups[sel] == reference, 0, 1)
            obs, p, nm = R.perm_test_alignment(D[sel], lb, rng)
            tests.append(dict(dataset=name, subset=sub, contrast=f"{reference} vs {o}",
                              n_ref=int((lb == 0).sum()), n_other=int((lb == 1).sum()),
                              diff_alignment=obs, p=p, null_mean=nm))
    return al, mag

# ---- musculo bajo insulina (GSE22309) ----
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); e = geo.probes_to_genes(e, geo.read_gpl_annot("GPL91.annot.gz"))
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"})
p["subj"] = np.arange(len(p)) // 2; p["run"] = p.title.str.extract(r"(Run\d+)")[0].fillna("none")
rec = []
for s in p.subj.unique():
    a = p[(p.subj == s) & (p.agent == "untreated")]; b = p[(p.subj == s) & (p.agent == "insulin")]
    if len(a) and len(b): rec.append((a.gsm.iloc[0], b.gsm.iloc[0], a.grp.iloc[0], a.run.iloc[0] == b.run.iloc[0], a.run.iloc[0]))
Z = R.embed(e); D = np.array([Z.loc[b].values - Z.loc[a].values for a, b, _, _, _ in rec])
grp = np.array([r[2] for r in rec]); same = np.array([r[3] for r in rec]); run = np.array([r[4] for r in rec])
al, mag = analyse("GSE22309_muscle_4h", D, grp, "IS", batches=same, subsets=("all", "sameRun"))
# control: pseudo-grupos por lote, ignorando el estado clinico
rs = pd.Series(run); big = [r for r, c in rs.value_counts().items() if c >= 8][:2]
if len(big) == 2:
    m = np.isin(run, big); lb = np.where(run[m] == big[0], 0, 1)
    obs, pv, nm = R.perm_test_alignment(D[m], lb, rng)
    tests.append(dict(dataset="GSE22309_muscle_4h", subset="control", contrast=f"pseudo-groups {big[0]} vs {big[1]}",
                      n_ref=int((lb == 0).sum()), n_other=int((lb == 1).sum()), diff_alignment=obs, p=pv, null_mean=nm))
# disociacion de la magnitud
r_, p_ = stats.spearmanr(mag, al)
tests.append(dict(dataset="GSE22309_muscle_4h", subset="all", contrast="rho(magnitude, alignment)",
                  n_ref=len(D), n_other=0, diff_alignment=r_, p=p_, null_mean=np.nan))
pd.DataFrame(rows).to_csv(f"{OUT}/alignment_per_person.tsv", sep="\t", index=False)
t = pd.DataFrame(tests); t.to_csv(f"{OUT}/alignment_tests.tsv", sep="\t", index=False)
print(t.round(4).to_string(index=False))
