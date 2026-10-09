#!/usr/bin/env python3
"""Analisis complementarios (resultados retirados o en perspectiva, reportados por transparencia):
(a) GSE66306 monocitos antes/3 meses tras cirugia (dispersion, Ic, test pareado);
(b) GSE59034 SAT antes/2 anos/nunca obesas (dispersion; confundido con lote longitudinal);
(c) sangre: profundidad frente a referencia sana en GSE156993 y GSE21321;
(d) GSE129843 restriccion horaria: coherencia del desplazamiento diurno R vs U (pareado).
Salida: results/supplementary/*.txt"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..")); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import geo, response as R, landscape as L
OUT = os.environ.get("T2D_OUT", "results/supplementary"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0); log = open(f"{OUT}/supplementary_results.txt", "w")
def P(*a): print(*a); print(*a, file=log)
def disp(A): return float(np.trace(np.cov(A.T)))
# (a) monocitos
try:
    e = pd.read_csv(geo.path("GSE66306_PM_processed_counts.txt.gz"), sep="\t").drop(columns=["Ensembl Gene ID"]).set_index("Gene Name").groupby(level=0).mean()
    X = np.log2(e + 1); X = X[(X > 3).mean(axis=1) > 0.5]; Z = R.embed(X, r=3); tp = np.array([c.split("_")[1] for c in X.columns]); subj = np.array([c.split("_")[0] for c in X.columns])
    Y = X.to_numpy().T; Yz = (Y - Y.mean(0)) / (Y.std(0) + 1e-9)
    P(f"GSE66306 monocitos: dispersion T0={disp(Z.values[tp=='T0']):.1f} T3={disp(Z.values[tp=='T3']):.1f}; Ic T0={L.critical_index(Yz[tp=='T0']):.2f} T3={L.critical_index(Yz[tp=='T3']):.2f}")
except FileNotFoundError as ex: P("GSE66306 omitido:", ex)
# (b) GSE59034
try:
    e = pd.read_csv(f"{os.environ.get('T2D_EXPORT','data/export')}/GSE59034/GSE59034_expr.tsv", sep="\t", index_col=0); p = pd.read_csv(f"{os.environ.get('T2D_EXPORT','data/export')}/GSE59034/GSE59034_pheno.tsv", sep="\t")
    Z = R.embed(e, r=3); cond = p.set_index(".sample_id").condition.loc[Z.index].values
    P("GSE59034 SAT dispersion:", {g: round(disp(Z.values[cond == g]), 0) for g in ["never_obese", "obese_before", "obese_after"]}, "(antes/despues confundido con lote; no se interpreta)")
except FileNotFoundError as ex: P("GSE59034 omitido (correr 01 primero):", ex)
# (c) sangre
def depth(name, e, g, ctrl):
    Z = R.embed(e, r=3).values; g = np.array(g); c = g == ctrl; h = L.silverman(Z[c])
    dep = np.array([-np.log(L.kde(Z[c & (np.arange(len(Z)) != i)], Z[i:i + 1], h)[0][0]) for i in range(len(Z))])
    P(f"{name}: profundidad media (ref={ctrl}):", {k: (round(dep[g == k].mean(), 2), int((g == k).sum())) for k in pd.unique(g)})
try:
    p, e = geo.read_series_matrix("GSE156993_series_matrix.txt.gz"); g = p.title.str.extract(r"^(\w+)")[0]; depth("GSE156993 PBMC", e[p.gsm], g, "H")
    p, e = geo.read_series_matrix("GSE21321-GPL6883_series_matrix.txt.gz"); g = p["patient type"].str.replace(r"\s*\d+$", "", regex=True); depth("GSE21321 sangre", e[p.gsm], g, g.unique()[0])
except FileNotFoundError as ex: P("sangre omitida:", ex)
# (d) TRF
try:
    e = pd.read_csv(geo.path("GSE129843_RESTRICT.txt.gz"), sep="\t").set_index("Symbol").drop(columns=["genes", "entrez"]).apply(pd.to_numeric, errors="coerce").groupby(level=0).mean(); e = e[(e > 1).mean(axis=1) > 0.5]
    Z = R.embed(e, r=5); cols = pd.Series(e.columns); m = pd.DataFrame({"col": cols, "subj": cols.str.split(".").str[0], "cond": cols.str.split(".").str[1], "tp": cols.str.split(".").str[2].str[1:].astype(int)})
    for k in sorted(m.tp.unique())[1:]:
        res = {}
        for c in ["R", "U"]:
            prs = [(m[(m.subj == s) & (m.cond == c) & (m.tp == 1)].col.iloc[0], m[(m.subj == s) & (m.cond == c) & (m.tp == k)].col.iloc[0]) for s in m.subj.unique() if ((m.subj == s) & (m.cond == c) & (m.tp == 1)).any() and ((m.subj == s) & (m.cond == c) & (m.tp == k)).any()]
            res[c] = R.coherence_loo(R.displacements(Z, prs))
        P(f"GSE129843 T1->T{k}: coherencia R={res['R']:.2f} U={res['U']:.2f}")
except FileNotFoundError as ex: P("GSE129843 omitido:", ex)
log.close()
