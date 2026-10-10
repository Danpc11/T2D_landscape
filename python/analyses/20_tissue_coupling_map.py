#!/usr/bin/env python3
"""Con que tejidos comparte posicion el musculo esqueletico (exploratorio, no entra en el articulo).

Extiende la coherencia individual de los tres organos metabolicos a 24 tejidos de GTEx, para
preguntar si el eje digestivo y el sistema nervioso participan del estado compartido entre tejidos
de una misma persona. Dos controles son imprescindibles porque ambos producen rankings espurios:
el numero de donantes, del que depende el nulo, y la similitud de composicion entre tejidos, que
hace que el musculo esqueletico se parezca a cualquier tejido con musculo liso.

Lectura de los resultados (octubre 2026, 24 tejidos, hasta 670 donantes):
  - En bruto, el nervio tibial es el primer acompanante del musculo, por encima del adiposo.
  - Tras ajustar por n y por similitud, el orden es nervio +0.09, metabolicos +0.04, musculo liso
    +0.03, digestivo 0.00, SNC -0.03, y NINGUNA diferencia alcanza significacion (SNC vs resto
    p = 0.15). El nervio es un unico tejido, de modo que su posicion no se distingue del ruido, y
    el nervio tibial de GTEx es rico en tejido adiposo y conectivo perineural.
  - El resultado solido es negativo: el SNC queda consistentemente por debajo, con el hipocampo el
    ultimo. La posicion individual compartida es un fenomeno de tejidos perifericos y no se
    extiende al sistema nervioso central.
Para poner a prueba una contribucion nerviosa o enterica harian falta disenos de intervencion
(bloqueo autonomico durante un clamp; comida mixta frente a insulina intravenosa en las mismas
personas con biopsias pareadas), no correlaciones entre tejidos de donantes post mortem.

Salida: results/gtex/muscle_partner_coupling.tsv
"""
import os, sys, itertools, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats
U = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0); FOCUS = os.environ.get("T2D_FOCUS", "muscle")
TN = {"muscle": "Muscle - Skeletal", "adipose": "Adipose - Subcutaneous", "adipose_visceral": "Adipose - Visceral (Omentum)",
      "pancreas": "Pancreas", "liver": "Liver", "adrenal": "Adrenal Gland", "kidney_cortex": "Kidney - Cortex", "stomach": "Stomach",
      "ileum": "Small Intestine - Terminal Ileum", "colon_transverse": "Colon - Transverse", "colon_sigmoid": "Colon - Sigmoid",
      "esophagus_mucosa": "Esophagus - Mucosa", "esophagus_muscularis": "Esophagus - Muscularis",
      "esophagus_gej": "Esophagus - Gastroesophageal Junction", "nerve_tibial": "Nerve - Tibial",
      "spinal_cord": "Brain - Spinal cord (cervical c-1)", "hypothalamus": "Brain - Hypothalamus",
      "nucleus_accumbens": "Brain - Nucleus accumbens (basal ganglia)", "caudate": "Brain - Caudate (basal ganglia)",
      "putamen": "Brain - Putamen (basal ganglia)", "cortex": "Brain - Cortex", "frontal_cortex": "Brain - Frontal Cortex (BA9)",
      "hippocampus": "Brain - Hippocampus", "breast": "Breast - Mammary Tissue"}
CLASS = {"nerve_tibial": "peripheral nerve", "spinal_cord": "CNS", "hypothalamus": "CNS", "cortex": "CNS",
         "frontal_cortex": "CNS", "hippocampus": "CNS", "caudate": "CNS", "putamen": "CNS", "nucleus_accumbens": "CNS",
         "adipose": "metabolic", "adipose_visceral": "metabolic", "pancreas": "metabolic", "liver": "metabolic",
         "stomach": "gut", "ileum": "gut", "colon_transverse": "gut", "colon_sigmoid": "gut", "esophagus_mucosa": "gut",
         "esophagus_muscularis": "smooth muscle", "esophagus_gej": "smooth muscle",
         "adrenal": "other", "kidney_cortex": "other", "breast": "other"}
T = {k: pd.read_pickle(f"{OUT}/{k}.pkl") for k in TN if os.path.exists(f"{OUT}/{k}.pkl")}
if FOCUS not in T: sys.exit(f"falta {OUT}/{FOCUS}.pkl: ejecutar antes 03a_gtex_prepare.py")
print(f"tejidos disponibles: {len(T)}")
ph = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
sa = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False)
sa["donor"] = sa.SAMPID.str.split("-").str[:2].str.join("-")
def cov(k, d):
    p = ph.loc[d]; age = p.AGE.str[:2].astype(float).values; sex = (p.SEX.values == 2).astype(float)
    h = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
    t = sa[(sa.SMTSD == TN[k]) & (sa.SMAFRZE == "RNASEQ")].drop_duplicates("donor").set_index("donor").reindex(d)
    return np.column_stack([np.ones(len(d)), age, sex, h, t.SMRIN.fillna(t.SMRIN.median()).values, t.SMTSISCH.fillna(t.SMTSISCH.median()).values])
def pcs(x, r=5):
    Z = (x.T - x.T.mean(0)) / (x.T.std(0) + 1e-9); U_, S, _ = np.linalg.svd(Z, full_matrices=False); return U_[:, :r] * S[:r]
def cca(A, B):
    qa, _ = np.linalg.qr(A - A.mean(0)); qb, _ = np.linalg.qr(B - B.mean(0)); return np.linalg.svd(qa.T @ qb, compute_uv=False)[0]
rows = []
for other in sorted(k for k in T if k != FOCUS):
    donors = sorted(set(T[FOCUS].columns) & set(T[other].columns))
    if len(donors) < 60: continue
    genes = sorted(set(T[FOCUS].index) & set(T[other].index))
    X = {}
    for k in (FOCUS, other):
        Y = T[k].loc[genes, donors].to_numpy().T; D = cov(k, donors)
        b, *_ = np.linalg.lstsq(D, Y, rcond=None); X[k] = (Y - D[:, 1:] @ b[1:]).T
    v = np.mean([x.var(1) / x.var(1).mean() for x in X.values()], 0); gi = np.argsort(-v)[:800]
    P = {k: pcs(x[gi]) for k, x in X.items()}; obs = cca(P[FOCUS], P[other])
    nl = [cca(P[FOCUS], P[other][rng.permutation(len(donors))]) for _ in range(200)]
    rows.append(dict(partner=other, partner_class=CLASS.get(other, "other"), n=len(donors),
                     rho=obs, null=float(np.mean(nl)), excess=obs - float(np.mean(nl))))
d = pd.DataFrame(rows)
# control 1: el nulo decrece con n, de modo que el exceso crece con n
Xn = np.column_stack([np.ones(len(d)), np.log(d.n)]); bn, *_ = np.linalg.lstsq(Xn, d.excess, rcond=None)
d["excess_adj_n"] = d.excess - Xn @ bn
# control 2: similitud de composicion, que hace que el musculo se parezca a cualquier tejido con musculo
mprof = T[FOCUS].mean(1)
d["similarity"] = [stats.spearmanr(mprof[sorted(set(mprof.index) & set(T[t].index))],
                                   T[t].loc[sorted(set(mprof.index) & set(T[t].index))].mean(1))[0] for t in d.partner]
Xs = np.column_stack([np.ones(len(d)), np.log(d.n), d.similarity]); bs, *_ = np.linalg.lstsq(Xs, d.excess, rcond=None)
d["excess_adj_n_similarity"] = d.excess - Xs @ bs
d = d.sort_values("excess", ascending=False); d.to_csv(f"{OUT}/{FOCUS}_partner_coupling.tsv", sep="\t", index=False)
print(d[["partner", "partner_class", "n", "rho", "excess", "excess_adj_n", "excess_adj_n_similarity"]].round(3).to_string(index=False))
print("\npor clase, tras ambos ajustes:")
print(d.groupby("partner_class").excess_adj_n_similarity.agg(["mean", "size"]).round(3).sort_values("mean", ascending=False).to_string())
cns = d[d.partner_class == "CNS"].excess_adj_n_similarity; rest = d[d.partner_class != "CNS"].excess_adj_n_similarity
print(f"\nSNC vs resto: {cns.mean():+.3f} vs {rest.mean():+.3f}   U p = {stats.mannwhitneyu(cns, rest).pvalue:.3f}")
