#!/usr/bin/env python3
"""Auditoria del pico de Fisher-Rao: IC bootstrap de la posicion del pico y nulo por barajado de theta.
Resultado (AUDIT.md C1): el pico cae donde la densidad de muestras es maxima tambien bajo el nulo -> retirado."""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..")); import landscape as L
E = os.environ.get("T2D_EXPORT", "data/export"); OUT = os.environ.get("T2D_OUT", "results/audit"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(1); rows = []
for acc, ctrl in [("GSE50244", "hba1c"), ("GSE50398", "hba1c"), ("GSE25462", "hba1c")]:
    try:
        e = pd.read_csv(f"{E}/{acc}/{acc}_expr.tsv", sep="\t", index_col=0); ph = pd.read_csv(f"{E}/{acc}/{acc}_pheno.tsv", sep="\t", dtype={".sample_id": str})
    except FileNotFoundError: continue
    col = [c for c in ph.columns if ctrl in c.lower()][0]; th = pd.to_numeric(ph[col], errors="coerce").to_numpy(); ok = ~np.isnan(th); ph = ph[ok]; th = th[ok]
    hv = e[ph[".sample_id"]].var(axis=1).sort_values(ascending=False).index[:800]; Y = e.loc[hv, ph[".sample_id"]].to_numpy().T; Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-9); Z, _ = L.pca(Y, 3); h = L.silverman(Z)
    c, g = L.fisher_along(Z, th, h); obs = c[np.argmax(g)]
    boots = [];
    for _ in range(300):
        i = rng.choice(len(th), len(th)); cb, gb = L.fisher_along(Z[i], th[i], h); boots.append(cb[np.argmax(gb)]) if len(gb) else None
    nulls = []
    for _ in range(300):
        cb, gb = L.fisher_along(Z, rng.permutation(th), h); nulls.append(cb[np.argmax(gb)]) if len(gb) else None
    rows.append(dict(acc=acc, control=ctrl, peak=obs, boot_lo=np.percentile(boots, 2.5), boot_hi=np.percentile(boots, 97.5), null_median=np.median(nulls), null_q25=np.percentile(nulls, 25), null_q75=np.percentile(nulls, 75), theta_median=np.median(th)))
df = pd.DataFrame(rows); df.to_csv(f"{OUT}/fisher_peak_audit.tsv", sep="\t", index=False); print(df.round(2).to_string(index=False))
