#!/usr/bin/env python3
"""Cuanta coherencia puede fabricar el ruido, y cuanta real se detectaria.
Responde a tres preguntas sobre la medida, usando los propios datos (GSE22309):
 (a) sesgo: que coherencia da el estadistico cuando no hay direccion compartida (desplazamientos
     con direcciones isotropas, conservando la longitud; cambiar solo el signo conservaria los ejes
     originales y daria un nulo demasiado benevolo) y como depende del tamano de grupo;
 (b) confusion: cuanta coherencia aparece entre individuos que solo comparten lote de hibridacion,
     y cuanto se desplaza la estimacion al restringir a pares del mismo lote;
 (c) potencia: que diferencia de coherencia se detecta con 80% de potencia a cada n.
Salida: results/power/power_and_confounding.tsv
"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import geo, response as R
OUT = os.environ.get("T2D_OUT", "results/power"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz")
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"}); p["subj"] = np.arange(len(p)) // 2
p["run"] = p.title.str.extract(r"(Run\d+)")[0].fillna("none")
Z = R.embed(e)
def pairs(sub): return [(p[(p.subj == s) & (p.agent == "untreated")].gsm.iloc[0], p[(p.subj == s) & (p.agent == "insulin")].gsm.iloc[0]) for s in sub]
D = {g: R.displacements(Z, pairs(p[p.grp == g].subj.unique())) for g in ["IS", "IR", "T2D"]}
rows = []
# (a) sesgo del estadistico frente al tamano de grupo, sin direccion compartida
for n in range(4, 21):
    vals_loo, vals_naive = [], []
    for _ in range(200):
        idx = rng.choice(len(D["IS"]), n, replace=True)
        Dr = R.isotropic_null(D["IS"][idx], rng)       # misma longitud, direcciones isotropas
        vals_loo.append(R.coherence_loo(Dr)); vals_naive.append(R.coherence_naive(Dr))
    rows.append(dict(analysis="null coherence vs n", n=n, loo_mean=np.mean(vals_loo), loo_p95=np.percentile(vals_loo, 95),
                     naive_mean=np.mean(vals_naive), naive_p95=np.percentile(vals_naive, 95)))
# (b) confusion por lote: coherencia entre individuos que solo comparten lote
same = p.groupby("subj").run.nunique().eq(1); keep = set(same[same].index)
for g in ["IS", "IR", "T2D"]:
    sub_all = p[p.grp == g].subj.unique(); sub_same = [s for s in sub_all if s in keep]
    rows.append(dict(analysis="batch restriction", group=g, n=len(sub_all), coherence=R.coherence_loo(R.displacements(Z, pairs(sub_all))),
                     n_same=len(sub_same), coherence_same=R.coherence_loo(R.displacements(Z, pairs(sub_same))) if len(sub_same) > 3 else np.nan))
# pseudo-grupos definidos solo por lote, ignorando el estado clinico
runs = p.groupby("subj").run.first()
for r_ in runs.value_counts().index[:4]:
    sub = runs[runs == r_].index
    if len(sub) < 5: continue
    rows.append(dict(analysis="pseudo-group by batch", group=str(r_), n=len(sub), coherence=R.coherence_loo(R.displacements(Z, pairs(sub)))))
# (c) potencia: diferencia detectable con 80% de potencia por permutacion
for n in [8, 11, 15, 20, 30, 40]:
    for delta in [0.15, 0.25, 0.35, 0.45]:
        hits = 0
        for _ in range(int(os.environ.get('T2D_POWER_REPS', 120))):
            A = D["IS"][rng.choice(len(D["IS"]), n, replace=True)]
            B = D["IS"][rng.choice(len(D["IS"]), n, replace=True)]
            w = rng.random(n) < delta            # fraccion de individuos con direccion aleatorizada
            B = np.where(w[:, None], R.isotropic_null(B, rng), B)
            _, pv = R.perm_test_coherence(A, B, rng, B=200)
            hits += pv < 0.05
        reps = int(os.environ.get('T2D_POWER_REPS', 120)); rows.append(dict(analysis="power", n=n, delta_fraction_randomised=delta, power=hits / reps))
    pd.DataFrame(rows).to_csv(f"{OUT}/power_and_confounding.tsv", sep="\t", index=False)
# (d) equivalencia de magnitudes (TOST sobre log-magnitud, margen 0.5 en log2)
from scipy import stats as _st
lm = {g: np.log2(np.linalg.norm(D[g], axis=1)) for g in D}
for g in ["IR", "T2D"]:
    d1 = lm["IS"]; d2 = lm[g]; diff = d2.mean() - d1.mean()
    se = np.sqrt(d1.var(ddof=1) / len(d1) + d2.var(ddof=1) / len(d2)); dfree = len(d1) + len(d2) - 2
    margin = 0.5   # un factor 1.41 en magnitud
    p_lo = _st.t.sf((diff + margin) / se, dfree); p_hi = _st.t.cdf((diff - margin) / se, dfree)
    rows.append(dict(analysis="magnitude equivalence (TOST)", group=g, diff_log2=diff, se=se, margin_log2=margin,
                     p_tost=max(p_lo, p_hi), equivalent=max(p_lo, p_hi) < 0.05))
df = pd.DataFrame(rows); df.to_csv(f"{OUT}/power_and_confounding.tsv", sep="\t", index=False)
for a in df.analysis.unique(): print("==", a); print(df[df.analysis == a].dropna(axis=1, how="all").round(3).to_string(index=False))
