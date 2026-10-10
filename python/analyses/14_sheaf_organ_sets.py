#!/usr/bin/env python3
"""Extension y control de especificidad del haz.

Pregunta: la arquitectura compartida que observamos entre musculo, adiposo y pancreas, es propia de
los organos del metabolismo, o la comparten dos tejidos cualesquiera de la misma persona? Sin esta
comparacion, el resultado podria no decir nada sobre metabolismo.

Compara la energia del haz, siempre sobre los mismos donantes y el mismo nulo de correspondencia
de genes, para: (a) los organos metabolicos; (b) un conjunto de organos no metabolicos; (c) conjuntos
mixtos. Tambien calcula la correlacion canonica individual para todos los pares disponibles.
Salida: results/gtex/sheaf_organ_sets.tsv, canonical_all_pairs.tsv
"""
import os, sys, itertools, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); import sheaf_coherence as SC
U = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex"); rng = np.random.default_rng(0)
NPERM = int(os.environ.get("T2D_SHEAF_PERM", 100)); NGENES = int(os.environ.get("T2D_SHEAF_GENES", 800)); BETA, R = 6, 8
ph = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt", sep="\t").set_index("SUBJID")
sa = pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt", sep="\t", low_memory=False); sa["donor"] = sa.SAMPID.str.split("-").str[:2].str.join("-")
TNAME = {"muscle": "Muscle - Skeletal", "adipose": "Adipose - Subcutaneous", "adipose_visceral": "Adipose - Visceral (Omentum)",
         "pancreas": "Pancreas", "liver": "Liver", "adrenal": "Adrenal Gland", "kidney_cortex": "Kidney - Cortex",
         "stomach": "Stomach", "ileum": "Small Intestine - Terminal Ileum"}
METAB = ["muscle", "adipose", "adipose_visceral", "pancreas", "liver"]
OTHER = ["adrenal", "kidney_cortex", "stomach", "ileum"]
T = {k: pd.read_pickle(f"{OUT}/{k}.pkl") for k in TNAME if os.path.exists(f"{OUT}/{k}.pkl")}
def covariates(k, d):
    p = ph.loc[d]; age = p.AGE.str[:2].astype(float).values; sex = (p.SEX.values == 2).astype(float)
    h = pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int), drop_first=True).values.astype(float)
    t = sa[(sa.SMTSD == TNAME[k]) & (sa.SMAFRZE == "RNASEQ")].drop_duplicates("donor").set_index("donor").reindex(d)
    rin = t.SMRIN.fillna(t.SMRIN.median()).values; isch = t.SMTSISCH.fillna(t.SMTSISCH.median()).values
    return np.column_stack([np.ones(len(d)), age, sex, h, rin, isch])
def residualise(keys, donors, genes):
    X = {}
    for k in keys:
        Y = T[k].loc[genes, donors].to_numpy().T; D = covariates(k, donors)
        b, *_ = np.linalg.lstsq(D, Y, rcond=None); X[k] = (Y - D[:, 1:] @ b[1:]).T
    return X
def sheaf(keys, label):
    donors = sorted(set.intersection(*[set(T[k].columns) for k in keys]))
    genes = sorted(set.intersection(*[set(T[k].index) for k in keys]))
    if len(donors) < 40: return None
    X = residualise(keys, donors, genes)
    v = np.mean([x.var(1) / x.var(1).mean() for x in X.values()], 0)
    gi = np.argsort(-v)[:NGENES]; Xh = {k: x[gi] for k, x in X.items()}
    E = [SC.spectral_embedding(SC.adjacency(x, BETA), R) for x in Xh.values()]
    ref = E[0]
    for _ in range(3): ref = np.mean([SC.procrustes_align(e, ref) for e in E], 0)
    obs = SC.sheaf_energy([SC.procrustes_align(e, ref) for e in E])
    null = []
    for _ in range(NPERM):
        E2 = [E[0]] + [e[rng.permutation(NGENES)] for e in E[1:]]; rf = E2[0]
        for _ in range(3): rf = np.mean([SC.procrustes_align(e, rf) for e in E2], 0)
        null.append(SC.sheaf_energy([SC.procrustes_align(e, rf) for e in E2]))
    z = (obs - np.mean(null)) / np.std(null)
    print(f"  {label:34s} n={len(donors):4d} organs={len(keys)}  E={obs:.3f}  null={np.mean(null):.3f}±{np.std(null):.3f}  z={z:.1f}")
    return dict(set=label, organs=",".join(keys), n_organs=len(keys), n_donors=len(donors), n_genes=NGENES,
                observed=obs, null_mean=float(np.mean(null)), null_sd=float(np.std(null)), z=z, n_perm=NPERM)
rows = []
print("energia del haz por conjunto de organos (mismo nulo, mismos parametros):")
for keys, label in [(["muscle", "adipose", "pancreas"], "metabolic (original three)"),
                    ([k for k in METAB if k in T], "metabolic (five)"),
                    ([k for k in OTHER if k in T], "non-metabolic control"),
                    ([k for k in ["muscle", "adipose", "stomach"] if k in T], "mixed: muscle, adipose, stomach"),
                    ([k for k in ["adrenal", "kidney_cortex", "pancreas"] if k in T], "mixed: adrenal, kidney, pancreas")]:
    r = sheaf(keys, label)
    if r: rows.append(r)
pd.DataFrame(rows).to_csv(f"{OUT}/sheaf_organ_sets.tsv", sep="\t", index=False)
# correlacion canonica individual para todos los pares disponibles
def pcs(x, r=5):
    Z = (x.T - x.T.mean(0)) / (x.T.std(0) + 1e-9); U_, S, Vt = np.linalg.svd(Z, full_matrices=False); return U_[:, :r] * S[:r]
def cca_r(A, B):
    qa, _ = np.linalg.qr(A - A.mean(0)); qb, _ = np.linalg.qr(B - B.mean(0)); return np.linalg.svd(qa.T @ qb, compute_uv=False)[0]
pairs = []
print("\ncorrelacion canonica individual por par de organos:")
for a, b in itertools.combinations([k for k in TNAME if k in T], 2):
    donors = sorted(set(T[a].columns) & set(T[b].columns))
    if len(donors) < 60: continue
    genes = sorted(set(T[a].index) & set(T[b].index)); X = residualise([a, b], donors, genes)
    v = np.mean([x.var(1) / x.var(1).mean() for x in X.values()], 0); gi = np.argsort(-v)[:NGENES]
    P = {k: pcs(x[gi]) for k, x in X.items()}
    obs = cca_r(P[a], P[b]); nl = [cca_r(P[a], P[b][rng.permutation(len(donors))]) for _ in range(200)]
    both_metab = (a in METAB) and (b in METAB)
    pairs.append(dict(organ_a=a, organ_b=b, both_metabolic=both_metab, n_donors=len(donors), rho=obs,
                      null_mean=float(np.mean(nl)), null_sd=float(np.std(nl)), p=(np.sum(np.array(nl) >= obs) + 1) / 201))
    print(f"  {a:17s} {b:17s} n={len(donors):4d}  rho={obs:.3f}  null={np.mean(nl):.3f}  {'metabolic' if both_metab else ''}")
pd.DataFrame(pairs).to_csv(f"{OUT}/canonical_all_pairs.tsv", sep="\t", index=False)
