#!/usr/bin/env python3
"""Prepara los insumos del panel de redes de la Fig 1: una red de coexpresion de ejemplo
(musculo sano, 150 genes de mayor varianza, bicor^beta) con su disposicion por fuerzas y sus
modulos, mas la copia de las metricas de red por estadio a n igual.
Salida: results/networks/<acc>_<cond>_network.pkl y all_metrics_combined.tsv"""
import os, sys, pickle, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
E = os.environ.get("T2D_EXPORT", "data/export"); OUT = os.environ.get("T2D_OUT", "results/networks"); os.makedirs(OUT, exist_ok=True)
ACC, COND, BETA, NGENES = os.environ.get("T2D_NET_ACC", "GSE25462"), os.environ.get("T2D_NET_COND", "ND"), 6, 150
e = pd.read_csv(f"{E}/{ACC}/{ACC}_expr.tsv", sep="\t", index_col=0); p = pd.read_csv(f"{E}/{ACC}/{ACC}_pheno.tsv", sep="\t", dtype={".sample_id": str})
X = e[p[p.condition == COND][".sample_id"]]; hv = X.var(axis=1).sort_values(ascending=False).index[:NGENES]
C = np.corrcoef(X.loc[hv].to_numpy()); A = np.abs(C) ** BETA; np.fill_diagonal(A, 0)
rng = np.random.default_rng(1); n = len(A); pos = rng.normal(size=(n, 2)); W = A / A.max()
for it in range(400):   # layout por fuerzas (atraccion ~ peso, repulsion ~ 1/d^2)
    d = pos[:, None, :] - pos[None, :, :]; dist = np.linalg.norm(d, axis=2) + 1e-3
    pos += ((d / dist[:, :, None] ** 3).sum(1) * 0.02 - (d * W[:, :, None]).sum(1) * 0.06) * (1 - it / 400); pos -= pos.mean(0)
D = 1 - np.abs(C); np.fill_diagonal(D, 0); mod = fcluster(linkage(squareform(D, checks=False), "average"), 4, "maxclust")
pickle.dump(dict(pos=pos, A=A, mod=mod, genes=list(hv), acc=ACC, cond=COND, beta=BETA), open(f"{OUT}/{ACC}_{COND}_network.pkl", "wb"))
src = os.path.join(os.environ.get("T2D_DISCOVERY", "results"), "summary", "all_metrics_combined.tsv")
if os.path.exists(src): pd.read_csv(src, sep="\t").to_csv(f"{OUT}/all_metrics_combined.tsv", sep="\t", index=False)
print(f"{ACC} {COND}: {n} nodes, {(A > 0.02).sum() // 2} edges above 0.02, {len(set(mod))} modules")
