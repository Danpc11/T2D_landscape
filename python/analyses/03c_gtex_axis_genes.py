#!/usr/bin/env python3
"""Genes del estado sistemico individual (eje canonico musculo-adiposo, GTEx). Requiere 03a. Salida: results/gtex/systemic_axis_genes.tsv"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
U = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex")
ph = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
sa = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False); sa["donor"] = sa.SAMPID.str.split("-").str[:2].str.join("-")
TNAME = {"muscle": "Muscle - Skeletal", "adipose": "Adipose - Subcutaneous", "pancreas": "Pancreas"}
T = {k: pd.read_pickle(f"{OUT}/{k}.pkl") for k in TNAME}; donors = sorted(set.intersection(*[set(t.columns) for t in T.values()])); genes = sorted(set.intersection(*[set(t.index) for t in T.values()]))
def design(k, d):
    p = ph.loc[d]; age = p.AGE.str[:2].astype(float).values; sex = (p.SEX.values == 2).astype(float); h = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
    t = sa[(sa.SMTSD == TNAME[k]) & (sa.SMAFRZE == "RNASEQ")].drop_duplicates("donor").set_index("donor").reindex(d); rin = t.SMRIN.fillna(t.SMRIN.median()).values; isch = t.SMTSISCH.fillna(t.SMTSISCH.median()).values
    return np.column_stack([np.ones(len(d)), age, sex, h, rin, isch])
X = {}
for k, t in T.items():
    Y = t.loc[genes, donors].to_numpy().T; D = design(k, donors); b, *_ = np.linalg.lstsq(D, Y, rcond=None); X[k] = Y - D[:, 1:] @ b[1:]
def pcs(Yr, r=5): Z = (Yr - Yr.mean(0)) / (Yr.std(0) + 1e-9); U_, S, Vt = np.linalg.svd(Z, full_matrices=False); return U_[:, :r] * S[:r]
P = {k: pcs(x) for k, x in X.items()}
def cca(A, B): qa, _ = np.linalg.qr(A - A.mean(0)); qb, _ = np.linalg.qr(B - B.mean(0)); u, s, vt = np.linalg.svd(qa.T @ qb); return qa @ u[:, 0], qb @ vt[0], s[0]
ua, ub, rho = cca(P["muscle"], P["adipose"]); shared = (ua + ub) / 2; print("rho musculo-adiposo:", round(rho, 3))
out = {}
for k in TNAME:
    Z = (X[k] - X[k].mean(0)) / (X[k].std(0) + 1e-9); out[k] = pd.Series((Z.T @ ((shared - shared.mean()) / shared.std())) / len(shared), index=genes)
df = pd.DataFrame(out); df["min_abs"] = df[list(TNAME)].abs().min(axis=1); df["same_sign"] = (np.sign(df.muscle) == np.sign(df.adipose)) & (np.sign(df.muscle) == np.sign(df.pancreas))
df.sort_values("min_abs", ascending=False).to_csv(f"{OUT}/systemic_axis_genes.tsv", sep="\t"); print(df[df.same_sign].sort_values("min_abs", ascending=False).head(30).round(2).to_string())
