# Genera export_sheaf/ sintetico con la estructura real (3 tejidos, n reales, IGT incoherente)
import numpy as np, pandas as pd, os
rng = np.random.default_rng(3); P=900; K=6
genes = [f"G{i:04d}" for i in range(P)]
m0 = np.repeat(np.arange(K), P//K)
def sim(m, n, a=1.0): 
    z = rng.normal(size=(n,K)); return (a*z[:,m] + rng.normal(size=(n,P))).T + 8
os.makedirs("export_sheaf_FAKE", exist_ok=True)
designs = {"GSE76895": (["ND","IGT","T2D"], [32,15,36]),
           "GSE18732": (["ND","IGT","T2D"], [47,26,45]),
           "GSE27951": (["NGT","IGT","T2D"], [11,11,11])}
m_t2d = m0.copy(); idx = rng.choice(P, int(.35*P), replace=False); m_t2d[idx] = rng.integers(0,K,idx.size)
for acc,(sts,ns) in designs.items():
    m_igt = m0.copy(); idx = rng.choice(P, int(.35*P), replace=False); m_igt[idx] = rng.integers(0,K,idx.size)
    blocks, ids, conds = [], [], []
    for st,n,m,a in zip(sts, ns, [m0,m_igt,m_t2d],[1,1,.7]):
        X = sim(m,n,a); blocks.append(X); ids += [f"{acc}_{st}_{k}" for k in range(n)]; conds += [st]*n
    E = pd.DataFrame(np.hstack(blocks), index=genes, columns=ids); E.index.name="gene"
    E.to_csv(f"export_sheaf_FAKE/{acc}_expr.tsv", sep="\t")
    pd.DataFrame({".sample_id": ids, "condition": conds}).to_csv(f"export_sheaf_FAKE/{acc}_pheno.tsv", sep="\t", index=False)
pd.DataFrame({"gene": genes}).to_csv("export_sheaf_FAKE/high_variance_genes_ordered.tsv", sep="\t", index=False)
