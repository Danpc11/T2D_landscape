#!/usr/bin/env python3
"""Figura 5 / R1: GTEx v11 tejidos emparejados. 03a prepara matrices por donante; 03b: coherencia del haz, coherencia individual (nulo de barajado de donantes), sangre -> tejido."""
import pandas as pd, numpy as np, warnings, sys; warnings.filterwarnings("ignore")
sys.path.insert(0,"/home/claude/fixed/python"); import sheaf_coherence as SC
import os; U=os.environ.get("T2D_RAW","data/raw"); OUT=os.environ.get("T2D_OUT","results/gtex"); os.makedirs(OUT, exist_ok=True); os.makedirs("gtex", exist_ok=True); rng=np.random.default_rng(0)
ph=pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt",sep="\t").set_index("SUBJID")
sa=pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",sep="\t",low_memory=False)
sa["donor"]=sa.SAMPID.str.split("-").str[:2].str.join("-")
TNAME={"muscle":"Muscle - Skeletal","adipose":"Adipose - Subcutaneous","pancreas":"Pancreas","blood":"Whole Blood"}
def sample_cov(k, d):
    """RIN, tiempo isquemico, lote de acido nucleico (one-hot de los lotes frecuentes) por muestra del tejido k"""
    t=sa[(sa.SMTSD==TNAME[k])&(sa.SMAFRZE=="RNASEQ")].drop_duplicates("donor").set_index("donor").reindex(d)
    rin=t.SMRIN.fillna(t.SMRIN.median()).values; isch=t.SMTSISCH.fillna(t.SMTSISCH.median()).values
    bt=t.SMNABTCH.fillna("NA"); top=bt.value_counts(); top=top[top>=8].index
    B=np.column_stack([(bt==b).astype(float).values for b in top]) if len(top)>1 else np.zeros((len(d),0))
    return np.column_stack([rin,isch,B])
T={k:pd.read_pickle(f"{OUT}/{k}.pkl") for k in ["muscle","adipose","pancreas"]}
donors=sorted(set.intersection(*[set(t.columns) for t in T.values()])); print("donantes emparejados:", len(donors))
genes=sorted(set.intersection(*[set(t.index) for t in T.values()]))
# residualizar sexo, edad (decada), Hardy
def design(d):
    p=ph.loc[d]; age=p.AGE.str[:2].astype(float).values; sex=(p.SEX.values==2).astype(float); h=pd.get_dummies(p.DTHHRDY.fillna(-1).astype(int),drop_first=True).values.astype(float)
    return np.column_stack([np.ones(len(d)),age,sex,h])
X={}
for k,t in T.items():
    Y=t.loc[genes,donors].to_numpy().T; D=np.column_stack([design(donors), sample_cov(k,donors)]); b,*_=np.linalg.lstsq(D,Y,rcond=None); X[k]=(Y-D[:,1:]@b[1:]).T   # genes x donantes
# genes HV comunes (varianza media)
v=np.mean([x.var(1)/x.var(1).mean() for x in X.values()],0); hv=np.array(genes)[np.argsort(-v)[:800]]; gi=[genes.index(g) for g in hv]
Xh={k:x[gi] for k,x in X.items()}
# ---------- (1) F1: haz de redes, 3 tejidos, mismos donantes ----------
beta=6; r=8
def embed(x): return SC.spectral_embedding(SC.adjacency(x,beta),r)
E=[embed(x) for x in Xh.values()]
ref=E[0]
for _ in range(3): ref=np.mean([SC.procrustes_align(e,ref) for e in E],0)
st=[SC.procrustes_align(e,ref) for e in E]; E_obs=SC.sheaf_energy(st)
null=[]
for b in range(100):
    E2=[E[0]]+[e[rng.permutation(800)] for e in E[1:]]; rf=E2[0]
    for _ in range(3): rf=np.mean([SC.procrustes_align(e,rf) for e in E2],0)
    null.append(SC.sheaf_energy([SC.procrustes_align(e,rf) for e in E2]))
print("[con RIN + tiempo isquemico + lote por muestra regresados]"); print(f"F1 GTEx (n={len(donors)}): E={E_obs:.4f}  nulo correspondencia={np.mean(null):.4f}±{np.std(null):.4f}  z={(E_obs-np.mean(null))/np.std(null):.1f}")
# ---------- (2) coherencia DENTRO de la persona: scores por donante y tejido ----------
def pcs(x,r=5):
    Z=(x.T-x.T.mean(0))/ (x.T.std(0)+1e-9); U_,S,Vt=np.linalg.svd(Z,full_matrices=False); return U_[:,:r]*S[:r]
P={k:pcs(x) for k,x in Xh.items()}
def cca_r(A,B):
    # primera correlacion canonica
    qa,_=np.linalg.qr(A-A.mean(0)); qb,_=np.linalg.qr(B-B.mean(0)); return np.linalg.svd(qa.T@qb,compute_uv=False)[0]
pairs=[("muscle","adipose"),("muscle","pancreas"),("adipose","pancreas")]
for a,b in pairs:
    obs=cca_r(P[a],P[b]); nl=[cca_r(P[a],P[b][rng.permutation(len(donors))]) for _ in range(500)]
    print(f"  coherencia individual {a}-{b}: rho_canonica={obs:.3f}  nulo(donantes barajados)={np.mean(nl):.3f}±{np.std(nl):.3f}  p={(np.sum(np.array(nl)>=obs)+1)/501:.3f}")
# ---------- (3) sangre -> tejidos (viabilidad clinica) ----------
B=pd.read_pickle(f"{OUT}/blood.pkl"); db=[d for d in donors if d in B.columns]; print("donantes con sangre ademas:", len(db))
gb=[g for g in genes if g in B.index]; Yb=B.loc[gb,db].to_numpy().T; Db=np.column_stack([design(db), sample_cov("blood",db)]); bb,*_=np.linalg.lstsq(Db,Yb,rcond=None); Yb=Yb-Db[:,1:]@bb[1:]
vb=Yb.var(0); Zb=Yb[:,np.argsort(-vb)[:800]]; Zb=(Zb-Zb.mean(0))/(Zb.std(0)+1e-9); Ub,Sb,_=np.linalg.svd(Zb,full_matrices=False); PB=Ub[:,:10]*Sb[:10]
idx=[donors.index(d) for d in db]
from numpy.linalg import lstsq
for k in ["muscle","adipose","pancreas"]:
    Y=P[k][idx]; # validacion cruzada 5-fold ridge-less (10 PCs de sangre -> 5 PCs tejido)
    n=len(db); folds=np.array_split(rng.permutation(n),5); r2=[]
    for f in folds:
        tr=np.setdiff1d(np.arange(n),f); A=np.column_stack([np.ones(len(tr)),PB[tr]]); w,*_=lstsq(A,Y[tr],rcond=None)
        pred=np.column_stack([np.ones(len(f)),PB[f]])@w; r2.append(1-((Y[f]-pred)**2).sum()/((Y[f]-Y[tr].mean(0))**2).sum())
    print(f"  sangre -> {k}: R2 (5-fold, 10 PCs sangre -> 5 PCs tejido) = {np.mean(r2):.3f}")
