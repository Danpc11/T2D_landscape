#!/usr/bin/env python3
"""Modelo von Mises-Fisher de las direcciones de respuesta.

Los analisis anteriores describen la coordinacion con estadisticos (coherencia, alineamiento).
Aqui se ajusta el modelo generativo que les corresponde: la direccion de respuesta de cada persona
es una realizacion de una distribucion von Mises-Fisher sobre la esfera, con direccion media mu y
concentracion kappa. kappa cuantifica la coordinacion con un parametro interpretable, comparable
entre cohortes de distinto tamano, y permite una prueba de razon de verosimilitudes en lugar de una
permutacion. Salida: results/response/vmf_fits.tsv, vmf_tests.tsv
"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import optimize, special, stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo, response as R
OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
def log_norm_const(k, d):
    """log C_d(kappa) usando ive para evitar desbordamiento con kappa grande"""
    return (d / 2 - 1) * np.log(k) - (d / 2) * np.log(2 * np.pi) - np.log(special.ive(d / 2 - 1, k)) - k
def fit_vmf(U):
    n, d = U.shape; m = U.sum(0); Rb = np.linalg.norm(m) / n
    k = optimize.minimize_scalar(lambda x: -n * (log_norm_const(x, d) + x * Rb), bounds=(1e-3, 500), method="bounded").x
    return dict(mu=m / np.linalg.norm(m), kappa=float(k), Rbar=float(Rb), n=n, loglik=float(n * (log_norm_const(k, d) + k * Rb)))
def lrt_equal_kappa(UA, UB):
    """H0: ambos grupos comparten kappa (direcciones libres). H1: kappa distinta."""
    d = UA.shape[1]; fa, fb = fit_vmf(UA), fit_vmf(UB)
    def joint(k): return -sum(len(S) * (log_norm_const(k, d) + k * np.linalg.norm(S.sum(0)) / len(S)) for S in (UA, UB))
    k0 = optimize.minimize_scalar(joint, bounds=(1e-3, 500), method="bounded").x
    stat = 2 * ((fa["loglik"] + fb["loglik"]) - (-joint(k0)))
    return dict(LR=float(stat), p=float(stats.chi2.sf(stat, 1)), kappa_pooled=float(k0))
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); e = geo.probes_to_genes(e, geo.read_gpl_annot("GPL91.annot.gz"))
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"})
p["subj"] = np.arange(len(p)) // 2; p["run"] = p.title.str.extract(r"(Run\d+)")[0].fillna("none")
rec = []
for s in p.subj.unique():
    a = p[(p.subj == s) & (p.agent == "untreated")]; b = p[(p.subj == s) & (p.agent == "insulin")]
    if len(a) and len(b): rec.append((a.gsm.iloc[0], b.gsm.iloc[0], a.grp.iloc[0], a.run.iloc[0] == b.run.iloc[0]))
Z = R.embed(e); D = np.array([Z.loc[b].values - Z.loc[a].values for a, b, _, _ in rec])
grp = np.array([r[2] for r in rec]); same = np.array([r[3] for r in rec])
U = D / np.linalg.norm(D, axis=1, keepdims=True)
fits, tests = [], []
for sub, mask in [("all", np.ones(len(U), bool)), ("sameRun", same)]:
    for g in ["IS", "IR", "T2D"]:
        S = U[mask & (grp == g)]
        if len(S) < 5: continue
        f = fit_vmf(S); bs = [fit_vmf(S[rng.integers(0, len(S), len(S))])["kappa"] for _ in range(400)]
        fits.append(dict(subset=sub, group=g, n=f["n"], kappa=f["kappa"], kappa_lo=float(np.percentile(bs, 2.5)),
                         kappa_hi=float(np.percentile(bs, 97.5)), Rbar=f["Rbar"], dim=U.shape[1]))
    for g in ["IR", "T2D"]:
        A = U[mask & (grp == "IS")]; B = U[mask & (grp == g)]
        if len(A) < 5 or len(B) < 5: continue
        t = lrt_equal_kappa(A, B); t.update(subset=sub, contrast=f"IS vs {g}", n_IS=len(A), n_other=len(B)); tests.append(t)
pd.DataFrame(fits).to_csv(f"{OUT}/vmf_fits.tsv", sep="\t", index=False)
pd.DataFrame(tests).to_csv(f"{OUT}/vmf_tests.tsv", sep="\t", index=False)
print(pd.DataFrame(fits).round(3).to_string(index=False)); print(); print(pd.DataFrame(tests).round(4).to_string(index=False))
