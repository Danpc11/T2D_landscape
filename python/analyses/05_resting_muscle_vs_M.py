#!/usr/bin/env python3
"""Figura 1: el transcriptoma basal del musculo frente al valor M del clamp (GSE182120, HTA 2.0, 2 laboratorios).
Salidas: results/resting/GSE182120_genomewide_M.tsv, GSE182120_clock_basal.tsv, summary.txt"""
import sys, os, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..")); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib"))
import geo, response as R, landscape as L
OUT = os.environ.get("T2D_OUT", "results/resting"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
p, e = geo.read_series_matrix("GSE182120_series_matrix.txt.gz"); amap = geo.read_gpl_annot("GPL17586-45144.txt"); g = geo.probes_to_genes(e, amap)
p["M"] = pd.to_numeric(p["m-value"]); p["bmi"] = pd.to_numeric(p.bmi); p["age"] = pd.to_numeric(p.age); lab = (p.lab == "IRS1").astype(float).values; t2d = (p.disease == "T2D").astype(int).values
D = np.column_stack([np.ones(len(p)), lab, p.age.values, p.bmi.values])
Y = g[p.gsm].to_numpy().T; keep = (Y > np.percentile(Y, 30)).mean(0) > 0.5; Y = Y[:, keep]; genes = g.index[keep]
b, *_ = np.linalg.lstsq(D, Y, rcond=None); Rz = Y - D[:, 1:] @ b[1:]
# correlacion parcial: el valor M tambien se residualiza sobre las mismas covariables,
# en lugar de correlacionar expresion residualizada con M crudo
bM, *_ = np.linalg.lstsq(D, p.M.to_numpy(float), rcond=None); Mres = p.M.to_numpy(float) - D[:, 1:] @ bM[1:]
res = np.array([stats.spearmanr(Rz[:, j], Mres) for j in range(Rz.shape[1])]); tt = np.array([stats.ttest_ind(Rz[t2d == 1, j], Rz[t2d == 0, j], equal_var=False) for j in range(Rz.shape[1])])
gw = pd.DataFrame({"rho_M": res[:, 0], "p_M": res[:, 1], "fdr_M": R.bh(res[:, 1]), "t_T2D_vs_NGT": tt[:, 0], "p_T2D": tt[:, 1], "fdr_T2D": R.bh(tt[:, 1])}, index=genes).sort_values("p_M")
gw.to_csv(f"{OUT}/GSE182120_genomewide_M.tsv", sep="\t")
clock = ["DBP", "TEF", "HLF", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "BHLHE40", "ARNTL", "CLOCK", "CRY1", "NFIL3", "TXNIP", "KLF15", "PPP1R3B", "PPARGC1A", "PDK4"]
gw.loc[[c for c in clock if c in gw.index]].to_csv(f"{OUT}/GSE182120_clock_basal.tsv", sep="\t")
# geometria de grupo (dispersion, Ic, profundidad) con covariables centradas por grupo
X = (Y - Y.mean(0)) / (Y.std(0) + 1e-9); v = X.var(0); X = X[:, np.argsort(-v)[:800]]
Dg = [np.ones(len(p)), lab]
for c in ["age", "bmi"]:
    vv = p[c].to_numpy().astype(float).copy(); [vv.__setitem__(t2d == k, vv[t2d == k] - vv[t2d == k].mean()) for k in (0, 1)]; Dg.append(vv)
Dg = np.column_stack(Dg); b2, *_ = np.linalg.lstsq(Dg, X, rcond=None); Xr = X - Dg[:, 1:] @ b2[1:]; Z, _ = L.pca(Xr, 3)
h = L.silverman(Z[t2d == 0]); dep = np.array([-np.log(L.kde(Z[(t2d == 0) & (np.arange(len(Z)) != i)], Z[i:i + 1], h)[0][0]) for i in range(len(Z))])
with open(f"{OUT}/summary.txt", "w") as f:
    f.write(f"n={len(p)} (NGT {int((t2d==0).sum())}, T2D {int(t2d.sum())}); M NGT {p.M[t2d==0].mean():.1f}±{p.M[t2d==0].std():.1f}, T2D {p.M[t2d==1].mean():.1f}±{p.M[t2d==1].std():.1f}\n")
    f.write(f"genes tested {len(gw)}; FDR<0.1 vs M: {(gw.fdr_M<0.1).sum()}; p<0.01: {(gw.p_M<0.01).sum()} (expected by chance ~{int(0.01*len(gw))})\n")
    f.write(f"FDR<0.1 T2D vs NGT: {(gw.fdr_T2D<0.1).sum()}\n")
    for k, nm in [(0, "NGT"), (1, "T2D")]: f.write(f"{nm}: dispersion={np.trace(np.cov(Z[t2d==k].T)):.0f} Ic={L.critical_index(Xr[t2d==k]):.2f}\n")
    r, pv = stats.spearmanr(Mres, dep); f.write(f"depth ~ M (partial, same covariates): rho={r:.2f} p={pv:.3f}\n")
    f.write(f"M-value units: as deposited in GSE182120 (verify against the source publication before submission)\n")
print(open(f"{OUT}/summary.txt").read()); print(gw.loc[[c for c in clock if c in gw.index], ["rho_M", "p_M", "t_T2D_vs_NGT", "p_T2D"]].round(3).to_string())
