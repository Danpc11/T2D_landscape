#!/usr/bin/env python3
"""Expresion diferencial clasica (T2D vs control, en reposo) y enriquecimiento por conjuntos curados,
por tejido. Es la comparacion canonica contra la que se contrasta el metodo de respuesta.
Salidas: results/de/de_<acc>.tsv, de_summary.tsv, de_set_enrichment.tsv"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import response as R
E = os.environ.get("T2D_EXPORT", "data/export"); OUT = os.environ.get("T2D_OUT", "results/de"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
SETS = {"immediate-early TFs": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "ATF3", "KLF10", "ZFP36", "NR4A1", "NR4A2", "NR4A3", "BTG2", "DUSP1"],
        "metallothionein / HMOX1": ["MT1A", "MT1E", "MT1F", "MT1G", "MT1H", "MT1M", "MT1X", "MT2A", "HMOX1", "NFE2L2", "TXNRD1"],
        "chaperones (UPR/HSP)": ["DNAJA1", "DNAJB1", "DNAJB4", "DNAJB5", "HSPA1A", "HSPA1B", "HSPA8", "HSPH1", "HSP90AA1", "XBP1", "ATF4", "DDIT3"],
        "amino-acid transport": ["SLC7A5", "SLC3A2", "SLC38A2", "SLC7A1", "SLC1A5", "SLC7A11"],
        "clock output": ["DBP", "TEF", "HLF", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART", "BHLHE40", "BHLHE41"],
        "canonical insulin targets": ["TXNIP", "KLF15", "IRS2", "PPP1R3B", "HES1", "PIK3R1", "SREBF1", "INSIG1", "PDK4", "FOXO1"],
        "OXPHOS / PGC-1α": ["PPARGC1A", "NDUFA4", "NDUFB8", "SDHB", "UQCRC1", "COX5A", "COX7A1", "ATP5F1A", "ATP5F1B", "CYCS", "ESRRA", "TFAM"],
        "inflammation / complement": ["C1QA", "C1QB", "C3AR1", "ITGB2", "IL1R1", "IL1RL1", "TNF", "CCL2", "STAT3", "S100A4", "LSP1", "HLA-DOA"],
        "ECM": ["COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL6A1", "COL6A2", "LTBP4", "FN1", "SPARC", "LAMA2"],
        "lipid / adipogenesis": ["PPARG", "ADIPOQ", "LEP", "FASN", "SCD", "CD36", "LPL", "PLIN1", "PLIN2", "CIDEC"]}
rows, setrows = [], []
for acc, tis in [("GSE164416", "islet"), ("GSE50244", "islet"), ("GSE25462", "muscle"), ("METSIM", "adipose")]:
    d = f"{E}/{acc}"
    if not os.path.exists(f"{d}/{acc}_expr.tsv"): print("skip", acc); continue
    e = pd.read_csv(f"{d}/{acc}_expr.tsv", sep="\t", index_col=0); p = pd.read_csv(f"{d}/{acc}_pheno.tsv", sep="\t", dtype={".sample_id": str})
    e = e[p[".sample_id"]]; cond = p.condition.values; lv = list(pd.unique(cond))
    hi = "T2D" if "T2D" in lv else ("obese" if "obese" in lv else lv[-1]); lo = "ND" if "ND" in lv else ("normal" if "normal" in lv else lv[0])
    A = e.loc[:, cond == lo].to_numpy(); B = e.loc[:, cond == hi].to_numpy()
    t, pv = stats.ttest_ind(B, A, axis=1); lfc = B.mean(1) - A.mean(1); fdr = R.bh(pv)
    pd.DataFrame({"gene": e.index, "logFC": lfc, "t": t, "p": pv, "fdr": fdr}).to_csv(f"{OUT}/de_{acc}.tsv", sep="\t", index=False)
    rows.append(dict(acc=acc, tissue=tis, contrast=f"{hi} vs {lo}", n_lo=int((cond == lo).sum()), n_hi=int((cond == hi).sum()), n_genes=len(e), n_fdr10=int((fdr < 0.1).sum()), n_p01=int((pv < 0.01).sum()), expected_p01=int(0.01 * len(e))))
    tt = pd.Series(np.abs(t), index=e.index)
    for name, g in SETS.items():
        gg = [x for x in g if x in tt.index]
        if len(gg) < 4: continue
        obs = tt.loc[gg].mean(); null = np.array([tt.sample(len(gg), random_state=int(rng.integers(1e9))).mean() for _ in range(2000)])
        setrows.append(dict(acc=acc, tissue=tis, set=name, n=len(gg), mean_t=obs, null_mean=null.mean(), null_sd=null.std(), z=(obs - null.mean()) / null.std(), p=(np.sum(null >= obs) + 1) / 2001))
if not rows:
    print("no cohorts found under", E, "- run 01 first"); raise SystemExit(1)
pd.DataFrame(rows).to_csv(f"{OUT}/de_summary.tsv", sep="\t", index=False)
se = pd.DataFrame(setrows); se["z"] = (se.mean_t - se.null_mean) / se.null_sd; se.to_csv(f"{OUT}/de_set_enrichment.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).to_string(index=False))
