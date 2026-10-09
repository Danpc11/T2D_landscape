#!/usr/bin/env python3
"""Convergencia in vitro (GSE182117, Gabriel 2021): efecto de glucosa+insulina cronica sobre la salida del reloj
en miotubos NGT vs T2D y amplitud circadiana (cosinor 24 h). Salida: results/myotubes/*.tsv"""
import sys, os, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import geo
OUT = os.environ.get("T2D_OUT", "results/myotubes"); os.makedirs(OUT, exist_ok=True)
m = geo.ensembl_map_from_gtex(); c = pd.read_csv(geo.path("GSE182117_counts.tsv.gz"), sep="\t", index_col=0)
c.index = [m.get(i.split(".")[0], i) for i in c.index]; c = c[~c.index.str.startswith("ENSG")].groupby(level=0).sum(); X = geo.log_cpm(c)
cols = pd.Series(X.columns); meta = pd.DataFrame({"col": cols, "dis": cols.str.extract(r"^(NGT|T2D)")[0], "trt": cols.str.extract(r"_(CTL|HGI)_")[0], "t": cols.str.extract(r"_(\d+)h_")[0].astype(int), "ind": cols.str.extract(r"_(i\d+)$")[0]})
genes = ["DBP", "TEF", "HLF", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "BHLHE40", "ARNTL", "CLOCK", "CRY1", "CRY2", "NFIL3", "TXNIP", "KLF15", "PPP1R3B", "HES1"]
rows = []
for gene in genes:
    if gene not in X.index: continue
    rec = dict(gene=gene)
    for dis in ["NGT", "T2D"]:
        d = []
        for ind in meta[meta.dis == dis].ind.unique():
            a = meta[(meta.dis == dis) & (meta.ind == ind) & (meta.trt == "CTL")]; b = meta[(meta.dis == dis) & (meta.ind == ind) & (meta.trt == "HGI")]; tt = sorted(set(a.t) & set(b.t))
            if len(tt) >= 3: d.append(np.mean([X.loc[gene, b[b.t == t_].col.iloc[0]] - X.loc[gene, a[a.t == t_].col.iloc[0]] for t_ in tt]))
        d = np.array(d); rec[f"dHGI_{dis}"] = d.mean(); rec[f"t_{dis}"] = d.mean() / (d.std(ddof=1) / np.sqrt(len(d)) + 1e-9)
        A = []
        for ind in meta[meta.dis == dis].ind.unique():
            s = meta[(meta.dis == dis) & (meta.ind == ind) & (meta.trt == "CTL")].sort_values("t")
            if len(s) < 6: continue
            y = X.loc[gene, s.col].values; t = s.t.values; Dm = np.column_stack([np.ones(len(t)), np.cos(2 * np.pi * t / 24), np.sin(2 * np.pi * t / 24)]); bb, *_ = np.linalg.lstsq(Dm, y, rcond=None); A.append(np.hypot(bb[1], bb[2]))
        rec[f"amp_{dis}"] = np.mean(A); rec[f"amp_{dis}_sd"] = np.std(A)
    rows.append(rec)
df = pd.DataFrame(rows).set_index("gene"); df.to_csv(f"{OUT}/GSE182117_clock_HGI_and_amplitude.tsv", sep="\t"); print(df.round(2).to_string())
