#!/usr/bin/env python3
"""Figura 5 / R1: GTEx v11 tejidos emparejados. 03a prepara matrices por donante; 03b: coherencia del haz, coherencia individual (nulo de barajado de donantes), sangre -> tejido."""
import os, sys, pandas as pd, numpy as np, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")); import sheaf_coherence as SC
U=os.environ.get("T2D_RAW","data/raw"); OUT=os.environ.get("T2D_OUT","results/gtex"); os.makedirs(OUT, exist_ok=True); os.makedirs("gtex", exist_ok=True); rng=np.random.default_rng(0)
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
NPERM=int(os.environ.get("T2D_SHEAF_PERM",200))
null=[]
for b in range(NPERM):
    E2=[E[0]]+[e[rng.permutation(800)] for e in E[1:]]; rf=E2[0]
    for _ in range(3): rf=np.mean([SC.procrustes_align(e,rf) for e in E2],0)
    null.append(SC.sheaf_energy([SC.procrustes_align(e,rf) for e in E2]))
pd.DataFrame({"draw": np.arange(len(null)), "energy": null}).to_csv(f"{OUT}/sheaf_null_draws.tsv", sep="\t", index=False)
z = (E_obs - np.mean(null)) / np.std(null)
pd.DataFrame([dict(analysis="sheaf_gtex", n_donors=len(donors), n_genes=len(hv), r=r, beta=beta, n_perm=len(null),
                   observed=E_obs, null_mean=float(np.mean(null)), null_sd=float(np.std(null)), z=z,
                   p=(np.sum(np.array(null) <= E_obs) + 1) / (len(null) + 1))]).to_csv(f"{OUT}/sheaf_summary.tsv", sep="\t", index=False)
print("[con RIN + tiempo isquemico + lote por muestra regresados]"); print(f"F1 GTEx (n={len(donors)}): E={E_obs:.4f}  nulo correspondencia={np.mean(null):.4f}±{np.std(null):.4f}  z={z:.1f}  ({len(null)} permutaciones; draws en sheaf_null_draws.tsv)")
# ---------- (2) coherencia DENTRO de la persona: scores por donante y tejido ----------
def pcs(x,r=5):
    Z=(x.T-x.T.mean(0))/ (x.T.std(0)+1e-9); U_,S,Vt=np.linalg.svd(Z,full_matrices=False); return U_[:,:r]*S[:r]
P={k:pcs(x) for k,x in Xh.items()}
# variantes canonicas por donante (para la figura de coordenadas paralelas)
def _cca_vars(A,B): qa,_=np.linalg.qr(A-A.mean(0)); qb,_=np.linalg.qr(B-B.mean(0)); u,s,vt=np.linalg.svd(qa.T@qb); return qa@u[:,0], qb@vt[0]
ua,ub=_cca_vars(P["muscle"],P["adipose"]); um_,upan=_cca_vars(P["muscle"],P["pancreas"]); upan=-upan if np.corrcoef(um_,ua)[0,1]<0 else upan
pd.DataFrame({"muscle":ua,"adipose":ub,"pancreas":upan},index=donors).to_csv(f"{OUT}/donor_scores.tsv",sep="\t")
def cca_r(A,B):
    # primera correlacion canonica
    qa,_=np.linalg.qr(A-A.mean(0)); qb,_=np.linalg.qr(B-B.mean(0)); return np.linalg.svd(qa.T@qb,compute_uv=False)[0]
pairs=[("muscle","adipose"),("muscle","pancreas"),("adipose","pancreas")]
cca_rows=[]; cca_draws=[]
for a,b in pairs:
    obs=cca_r(P[a],P[b]); nl=[cca_r(P[a],P[b][rng.permutation(len(donors))]) for _ in range(500)]
    cca_rows.append(dict(pair=f"{a}-{b}", observed=obs, null_mean=float(np.mean(nl)), null_sd=float(np.std(nl)), n_perm=len(nl), p=(np.sum(np.array(nl)>=obs)+1)/(len(nl)+1)))
    cca_draws += [dict(pair=f"{a}-{b}", draw=i, rho=v) for i,v in enumerate(nl)]
    print(f"  coherencia individual {a}-{b}: rho_canonica={obs:.3f}  nulo(donantes barajados)={np.mean(nl):.3f}±{np.std(nl):.3f}  p={cca_rows[-1]['p']:.3f}")
pd.DataFrame(cca_rows).to_csv(f"{OUT}/canonical_summary.tsv",sep="\t",index=False); pd.DataFrame(cca_draws).to_csv(f"{OUT}/canonical_null_draws.tsv",sep="\t",index=False)
# ---------- (3) sangre -> tejidos (viabilidad clinica) ----------
B=pd.read_pickle(f"{OUT}/blood.pkl"); db=[d for d in donors if d in B.columns]; print("donantes con sangre ademas:", len(db))
gb=[g for g in genes if g in B.index]; Yb_raw=B.loc[gb,db].to_numpy().T; Db=np.column_stack([design(db), sample_cov("blood",db)])
idx=[donors.index(d) for d in db]; Xt={k: Xh[k][:, idx].T for k in ["muscle","adipose","pancreas"]}   # donantes x genes
from numpy.linalg import lstsq
rows=[]
n=len(db); folds=np.array_split(rng.permutation(n),5)
for k in ["muscle","adipose","pancreas"]:
    r2=[]
    for f in folds:
        tr=np.setdiff1d(np.arange(n),f)
        # ---- todo lo que sigue se aprende SOLO con el training fold ----
        bb,*_=lstsq(Db[tr],Yb_raw[tr],rcond=None); Rtr=Yb_raw[tr]-Db[tr][:,1:]@bb[1:]; Rte=Yb_raw[f]-Db[f][:,1:]@bb[1:]
        sel=np.argsort(-Rtr.var(0))[:800]; mu=Rtr[:,sel].mean(0); sd=Rtr[:,sel].std(0)+1e-9
        Ztr=(Rtr[:,sel]-mu)/sd; Zte=(Rte[:,sel]-mu)/sd
        U_,S_,Vt=np.linalg.svd(Ztr,full_matrices=False); PBtr=U_[:,:10]*S_[:10]; PBte=Zte@Vt[:10].T
        # el objetivo (componentes del tejido) tambien se define en el training fold
        Ytr_raw=Xt[k][tr]; muY=Ytr_raw.mean(0); sdY=Ytr_raw.std(0)+1e-9
        Uy,Sy,Vty=np.linalg.svd((Ytr_raw-muY)/sdY,full_matrices=False); Ytr=Uy[:,:5]*Sy[:5]; Yte=((Xt[k][f]-muY)/sdY)@Vty[:5].T
        A=np.column_stack([np.ones(len(tr)),PBtr]); w,*_=lstsq(A,Ytr,rcond=None)
        pred=np.column_stack([np.ones(len(f)),PBte])@w
        r2.append(1-((Yte-pred)**2).sum()/((Yte-Ytr.mean(0))**2).sum())
    rows.append(dict(tissue=k, R2_mean=float(np.mean(r2)), R2_sd=float(np.std(r2)), folds=len(folds), n=n))
    print(f"  sangre -> {k}: R2 (5-fold, preprocesamiento dentro del fold) = {np.mean(r2):.3f} ± {np.std(r2):.3f}")
pd.DataFrame(rows).to_csv(f"{OUT}/blood_to_tissue.tsv",sep="\t",index=False)
