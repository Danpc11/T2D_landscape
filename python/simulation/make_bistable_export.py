"""Datos sinteticos con dos atractores en un espacio latente 2D + IGT en la silla,
proyectados a 900 genes con ruido, mas covariable BMI confundida."""
import numpy as np, pandas as pd, os
rng = np.random.default_rng(11); P = 900; genes = [f"G{i:04d}" for i in range(P)]
os.makedirs("bistable_export", exist_ok=True)
W = 2.5 * rng.normal(size=(2, P)) / np.sqrt(2)          # carga latente -> genes
designs = {"GSE76895": (["ND","IGT","T2D"], [32,15,36]), "GSE18732": (["ND","IGT","T2D"], [47,26,45]),
           "GSE27951": (["NGT","IGT","T2D"], [11,10,12])}
for acc, (sts, ns) in designs.items():
    Z, cond, bmi = [], [], []
    centers = {0: np.array([-2., 0.]), 2: np.array([2., 0.5])}
    for k, (st, n) in enumerate(zip(sts, ns)):
        if k == 1:   # silla: entre los valles, con dispersion alta (huella critica)
            z = np.array([0., 0.25]) + rng.normal(scale=[0.9, 0.6], size=(n, 2))
        else:
            z = centers[k] + rng.normal(scale=0.5, size=(n, 2))
        Z.append(z); cond += [st]*n
        bmi += list(26 + 3*k + rng.normal(scale=3, size=n))      # BMI sube con estadio (confusor)
    Z = np.vstack(Z); bmi = np.array(bmi)
    X = Z @ W + 0.3 * (bmi - bmi.mean())[:, None] * rng.normal(size=P)[None, :] * 0.05 + rng.normal(scale=1.0, size=(len(Z), P))
    ids = [f"{acc}_{c}_{i}" for i, c in enumerate(cond)]
    pd.DataFrame(X.T + 8, index=genes, columns=ids).rename_axis("gene").to_csv(f"bistable_export/{acc}_expr.tsv", sep="\t")
    pd.DataFrame({".sample_id": ids, "condition": cond, "bmi": bmi.round(1), "age": rng.integers(40, 70, len(ids)),
                  "hba1c": (5.2 + 0.9*np.array([sts.index(c) for c in cond]) + rng.normal(scale=0.3, size=len(ids))).round(2)})\
      .to_csv(f"bistable_export/{acc}_pheno.tsv", sep="\t", index=False)
pd.DataFrame({"gene": genes}).to_csv("bistable_export/high_variance_genes_ordered.tsv", sep="\t", index=False)
print("ok")
