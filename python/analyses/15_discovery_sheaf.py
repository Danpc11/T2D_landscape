#!/usr/bin/env python3
"""Paneles 1b,1c,1e: organizacion de red por estadio dentro de cada organo, y energia del haz entre
los tres organos de descubrimiento (cohortes distintas, personas distintas), a n igual por estadio.
Recalcula en Python lo que antes producia el pipeline de R, para que toda cifra del articulo salga
de una sola ejecucion. Salida: results/discovery/*.tsv
"""
import os, sys, gzip, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); import sheaf_coherence as SC
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo
U = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/discovery"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0); NG = int(os.environ.get("T2D_SHEAF_GENES", 800)); NPERM = int(os.environ.get("T2D_SHEAF_PERM", 100)); BETA, R = 6, 8
def stage_of(series, mapping):
    return pd.Series([mapping.get(str(v).strip().lower()) for v in series])
def load(acc, fname, stage_col, mapping, log2=False):
    p, e = geo.read_series_matrix(fname)
    st = stage_of(p[stage_col], mapping); keep = st.notna().values
    e = e[p.gsm[keep]]; st = st[keep].values
    if log2: e = np.log2(e.clip(lower=1))
    e = e[(e > np.percentile(e.values, 30)).mean(axis=1) > 0.5]
    return e, st
cfg = [("islet", "GSE76895_series_matrix.txt.gz", "diabetes status", {"nd": "healthy", "igt": "intermediate", "t2d": "T2D"}, False),
       ("muscle", "GSE18732_series_matrix.txt.gz", "glycemiagroup", {"1": "healthy", "2": "intermediate", "3": "T2D"}, False),
       ("adipose", "GSE27951-GPL570_series_matrix.txt.gz", "clinical status", {"ngt": "healthy", "igt": "intermediate", "dm": "T2D"}, True)]
data = {}
for name, fn, col, mp, lg in cfg:
    try:
        p, _ = geo.read_series_matrix(fn)
        hit = [c for c in p.columns if col in c.lower()] or [c for c in p.columns if c not in ("gsm", "title")]
        e, st = load(name, fn, hit[0], mp, lg)
        if pd.Series(st).nunique() >= 2: data[name] = (e, st); print(f"{name}: {e.shape}, estadios {pd.Series(st).value_counts().to_dict()}")
        else: print(f"{name}: no se pudo asignar estadio con la columna '{hit[0]}'")
    except Exception as ex: print(f"{name}: {type(ex).__name__} {ex}")
if len(data) < 2: sys.exit("no hay suficientes organos con estadio asignado")
genes = sorted(set.intersection(*[set(e.index) for e, _ in data.values()]))
v = np.mean([pd.DataFrame(e.loc[genes]).var(axis=1).to_numpy() / pd.DataFrame(e.loc[genes]).var(axis=1).mean() for e, _ in data.values()], 0)
hv = [genes[i] for i in np.argsort(-v)[:NG]]
rows_m, rows_s = [], []
for stage in ["healthy", "intermediate", "T2D"]:
    Xs, ns = {}, []
    for name, (e, st) in data.items():
        cols = e.columns[np.array(st) == stage]
        if len(cols) < 8: Xs = {}; break
        Xs[name] = e.loc[hv, cols].to_numpy(); ns.append(len(cols))
    if not Xs: print(f"{stage}: insuficientes muestras"); continue
    nmin = min(ns); reps = 10
    Es, mets = [], {k: [] for k in Xs}
    for rep in range(reps):
        sub = {k: x[:, rng.choice(x.shape[1], nmin, replace=False)] for k, x in Xs.items()}
        emb = []
        for k, x in sub.items():
            emb.append(SC.spectral_embedding(SC.adjacency(x, BETA), R))
        ref = emb[0]
        for _ in range(3): ref = np.mean([SC.procrustes_align(e_, ref) for e_ in emb], 0)
        Es.append(SC.sheaf_energy([SC.procrustes_align(e_, ref) for e_ in emb]))
    null = []
    sub = {k: x[:, rng.choice(x.shape[1], nmin, replace=False)] for k, x in Xs.items()}
    emb = [SC.spectral_embedding(SC.adjacency(x, BETA), R) for x in sub.values()]
    for _ in range(NPERM):
        E2 = [emb[0]] + [e_[rng.permutation(e_.shape[0])] for e_ in emb[1:]]; rf = E2[0]
        for _ in range(3): rf = np.mean([SC.procrustes_align(e_, rf) for e_ in E2], 0)
        null.append(SC.sheaf_energy([SC.procrustes_align(e_, rf) for e_ in E2]))
    z = (np.mean(Es) - np.mean(null)) / np.std(null)
    rows_s.append(dict(stage=stage, organs=",".join(Xs), n_per_stage=nmin, observed=float(np.mean(Es)), observed_sd=float(np.std(Es)),
                       null_mean=float(np.mean(null)), null_sd=float(np.std(null)), z=z, n_perm=NPERM, reps=reps))
    print(f"  {stage:13s} n={nmin:3d}  E={np.mean(Es):.3f}±{np.std(Es):.3f}  null={np.mean(null):.3f}±{np.std(null):.3f}  z={z:.1f}")
pd.DataFrame(rows_s).to_csv(f"{OUT}/discovery_sheaf_by_stage.tsv", sep="\t", index=False)
