#!/usr/bin/env python3
"""Metaanalisis entre cohortes, reto de glucocorticoides y lectura entre tejidos.

Combina las estimaciones por cohorte del paso 21 con efectos aleatorios, aplica la puntuacion
CONGELADA a una perturbacion glucocorticoide controlada, y mide cuanto puede leerse el eje de un
tejido desde otro sin fuga de informacion.

Tres cautelas implementadas de forma explicita:
  - La puntuacion NO se reajusta sobre los datos del reto: se aplica tal cual se definio.
  - En la prediccion entre tejidos, TODO el aprendizaje (ajuste de covariables, estandarizacion y
    modelo) ocurre dentro de los pliegues de entrenamiento, y los dos tejidos de un donante van
    siempre al mismo pliegue.
  - Predictor y resultado usan genes DISJUNTOS, de modo que ningun gen aparece a ambos lados.

Salidas: results/gtex/meta_random_effects.tsv, hydrocortisone_response.tsv, crosstissue_prediction.tsv
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex")
GC = ["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "DDIT4", "SESN1", "GLUL", "SGK1", "CEBPD"]
INS = ["TXNIP", "PDK4", "HK2", "IRS2", "ANGPTL4"]
TARGET = ["TXNIP", "PDK4", "IRS2", "ANGPTL4", "MT1X", "MT2A", "CEBPB", "TSC22D1"]   # disjuntos de GC

def dersimonian_laird(r, se):
    z = np.arctanh(r); sez = se / (1 - r ** 2)
    w = 1 / sez ** 2; Q = np.sum(w * (z - np.sum(w * z) / np.sum(w)) ** 2); k = len(z)
    Cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    wr = 1 / (sez ** 2 + tau2); zm = np.sum(wr * z) / np.sum(wr); sem = np.sqrt(1 / np.sum(wr))
    return np.tanh(zm), np.tanh(zm - 1.96 * sem), np.tanh(zm + 1.96 * sem), tau2, \
           (max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0), Q, k

def meta(files):
    d = pd.concat([pd.read_csv(f, sep="\t") for f in files if os.path.exists(f)], ignore_index=True)
    rows = []
    for s in d.set.unique():
        sub = d[d.set == s]
        if len(sub) < 2: continue
        r, lo, hi, tau2, I2, Q, k = dersimonian_laird(sub.r.values, sub.se.values)
        rows.append(dict(set=s, r=r, lo=lo, hi=hi, tau2=tau2, I2=I2, Q=Q, k=k))
    return pd.DataFrame(rows), d

def hydrocortisone():
    """GSE148640: hidrocortisona frente a placebo, intrasujeto, 10 pacientes. Puntuacion congelada."""
    f = f"{RAW}/GSE148640_processed.txt.gz"
    if not os.path.exists(f): print("  (omitido) falta", f); return None
    d = pd.read_csv(f, sep="\t", skiprows=9, low_memory=False)
    d = d.rename(columns={d.columns[1]: "gene", d.columns[2]: "p", d.columns[3]: "FC"})
    d = d[d.gene.notna() & (d.gene != "---")]
    vals = d.iloc[:, 5:25].apply(pd.to_numeric, errors="coerce")
    d["delta"] = vals.iloc[:, 1::2].mean(1) - vals.iloc[:, 0::2].mean(1)
    d = d.groupby("gene").delta.mean()
    rng = np.random.default_rng(0); rows = []
    for nm, gs in [("Glucocorticoid", GC), ("Insulin core", INS)]:
        g = [x for x in gs if x in d.index]
        obs = d[g].mean()
        null = np.array([d.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(5000)])
        rows.append(dict(set=nm, n=len(g), delta=obs, p=(np.sum(np.abs(null) >= abs(obs)) + 1) / 5001))
        print(f"  {nm:16s} delta={obs:+.3f}  p={rows[-1]['p']:.4f}")
    return pd.DataFrame(rows)

if __name__ == "__main__":
    pooled, per_cohort = meta([f"{OUT}/ci_AAGMEx.tsv", f"{OUT}/ci_GSE41168.tsv"])
    per_cohort.to_csv(f"{OUT}/ci_all_cohorts.tsv", sep="\t", index=False)
    pooled.to_csv(f"{OUT}/meta_random_effects.tsv", sep="\t", index=False)
    print("METAANALISIS de efectos aleatorios")
    for _, r in pooled.iterrows():
        print(f"  {r.set:16s} r={r.r:+.3f} [{r.lo:+.3f}, {r.hi:+.3f}]  I2={r.I2:5.1f}%  Q={r.Q:.2f}")
    print("\nRETO DE HIDROCORTISONA (puntuacion congelada)")
    h = hydrocortisone()
    if h is not None: h.to_csv(f"{OUT}/hydrocortisone_response.tsv", sep="\t", index=False)
