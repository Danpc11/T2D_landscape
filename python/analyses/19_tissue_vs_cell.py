#!/usr/bin/env python3
"""Donde reside la coordinacion: tejido o miocito.

Tres analisis que localizan el fenomeno.
1. Tipo de fibra. Con los marcadores validados contra cadena pesada de miosina, el eje lento-rapido
   no difiere entre grupos y ajustar por el deja el efecto intacto.
2. Compartimento mononuclear. Las proporciones se estiman por minimos cuadrados no negativos contra
   una firma construida desde scRNA-seq de vasto lateral humano. El ajuste se hace con covariables
   CENTRADAS, de modo que la media de los desplazamientos (y por tanto la direccion compartida) se
   conserva; ajustar sin centrar restaria un vector casi constante y destruiria el alineamiento por
   construccion. Las proporciones se transforman a log-ratio centrado porque suman uno. El control
   es el mismo ajuste con covariables aleatorias.
3. Miotubos. Miocitos puros de 24 donantes con insulina aguda: si la coordinacion fuera celular,
   deberia aparecer ahi.
Salida: results/tissue/*.tsv
"""
import os, sys, gzip, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import stats, optimize, special
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo, response as R
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/tissue"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0)
TYPE_I = ["TPM3","TNNC1","TNNI1","TNNT1","MYH7","MYL2","MYL3","MYL6B","ANKRD2","MYOZ2","ATP2A2","CASQ2","PLN","CD36","CYB5R1","FABP3","LDHB","CA3","MYL12A","PDLIM1"]
TYPE_IIA = ["TPM1","TNNC2","TNNI2","TNNT3","MYH2","MYH1","MYL1","MYLPF","MYBPC2","ENO3","ATP2A1","SLN","ALDOA","GAPDH","LDHA","PKM","PFKM","PGM1","DDIT4L","G0S2"]

def adjust(D, C):
    """Regresa covariables centradas y escaladas, conservando la media de D."""
    C = np.asarray(C, float); C = C.reshape(len(D), -1)
    C = (C - C.mean(0)) / (C.std(0) + 1e-12)
    b, *_ = np.linalg.lstsq(np.column_stack([np.ones(len(C)), C]), D, rcond=None)
    return D - C @ b[1:]

rows = []
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); e = geo.probes_to_genes(e, geo.read_gpl_annot("GPL91.annot.gz"))
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"}); p["subj"] = np.arange(len(p)) // 2
rec = []
for s in p.subj.unique():
    a = p[(p.subj == s) & (p.agent == "untreated")]; b = p[(p.subj == s) & (p.agent == "insulin")]
    if len(a) and len(b): rec.append((a.gsm.iloc[0], b.gsm.iloc[0], a.grp.iloc[0]))
grp = np.array([r[2] for r in rec]); base = e[[a for a, _, _ in rec]]
Z = R.embed(e); D = np.array([Z.loc[b].values - Z.loc[a].values for a, b, _ in rec]); lb = np.where(grp == "IS", 0, 1)
al = R.alignment_per_person(D, grp == "IS")
z = lambda gs: ((lambda M: (M - M.mean(1, keepdims=True)) / (M.std(1, keepdims=True) + 1e-9))(base.loc[[g for g in gs if g in base.index]].to_numpy())).mean(0)
fibre = z(TYPE_I) - z(TYPE_IIA)
prop_file = f"{OUT}/../response/deconvolution_proportions.tsv"
P = pd.read_csv(prop_file, sep="\t").drop(columns=["group"], errors="ignore") if os.path.exists(prop_file) else None
settings = [("unadjusted", None), ("fibre type (slow - fast)", fibre.reshape(-1, 1))]
if P is not None:
    L = np.log(P.to_numpy() + 1e-6); clr = (L - L.mean(1, keepdims=True))[:, :-1]
    settings += [("mononuclear composition", clr), ("fibre + mononuclear", np.column_stack([fibre, clr])),
                 ("control: random covariates", rng.normal(size=(len(D), clr.shape[1] + 1)))]
print("efecto tras ajustar por composicion (GSE22309):")
for name, C in settings:
    Dx = D if C is None else adjust(D, C)
    o, pv, _ = R.perm_test_alignment(Dx, lb, rng, B=2000)
    rows.append(dict(analysis="composition adjustment", setting=name, n=len(D), diff_alignment=o, p=pv))
    print(f"  {name:28s} dif = {o:+.3f}  p = {pv:.4f}")
rows.append(dict(analysis="fibre axis", setting="group difference (Kruskal)", n=len(D),
                 diff_alignment=np.nan, p=stats.kruskal(*[fibre[grp == g] for g in ["IS", "IR", "T2D"]])[1]))
rows.append(dict(analysis="fibre axis", setting="correlation with alignment", n=len(D),
                 diff_alignment=stats.spearmanr(fibre, al)[0], p=stats.spearmanr(fibre, al)[1]))

# ---- miotubos: miocito puro, insulina aguda, 24 donantes ----
def meta(path):
    rws = {}
    with gzip.open(path, "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"): gsm = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_title"): ti = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in l.rstrip().split("\t")[1:]]; k = v[0].split(":")[0].strip()
                if k not in rws: rws[k] = [x.split(":", 1)[1].strip() if ":" in x else "" for x in v]
    d = pd.DataFrame(rws); d["gsm"] = gsm; d["title"] = ti; return d
try:
    an = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot.tsv.gz", sep="\t", usecols=["GeneID", "Symbol"], index_col=0)
    M = pd.concat([meta(f"{RAW}/{a}_series_matrix.txt.gz") for a in
                   ["GSE81965-GPL16791", "GSE81965-GPL11154", "GSE63887-GPL16791", "GSE63887-GPL11154"]], ignore_index=True)
    cnt = pd.concat([pd.read_csv(f"{RAW}/{a}_raw_counts_GRCh38_p13_NCBI.tsv.gz", sep="\t", index_col=0) for a in ["GSE81965", "GSE63887"]], axis=1)
    cnt = cnt.loc[:, ~cnt.columns.duplicated()]; cnt.index = an.Symbol.reindex(cnt.index).values
    cnt = cnt[pd.notna(cnt.index)].groupby(level=0).sum(); X = geo.log_cpm(cnt)
    M["t"] = M.title.str.extract(r"_([0-9.]+)h$")[0].astype(float); M["donor"] = M["cell id"]; M["grp"] = M.disease.str.upper()
    M = M[M.gsm.isin(X.columns)].drop_duplicates("gsm")
    Zm = R.embed(X[M.gsm.tolist()]); Zm.index = M.gsm.values
    print(f"\nmiotubos: {M.donor.nunique()} donantes, insulina 100 nM")
    for TP in [0.5, 1.0, 2.0]:
        Dm, gm = [], []
        for dn in M.donor.unique():
            m = M[M.donor == dn]; a = m[m.t == 0]; b = m[m.t == TP]
            if len(a) and len(b): Dm.append(Zm.loc[b.gsm.iloc[0]].values - Zm.loc[a.gsm.iloc[0]].values); gm.append(m.grp.iloc[0])
        if len(Dm) < 8: continue
        Dm = np.array(Dm); gm = np.array(gm); lbm = np.where(gm == "NGT", 0, 1)
        o, pv, _ = R.perm_test_alignment(Dm, lbm, rng, B=3000)
        rows.append(dict(analysis="myotubes", setting=f"{TP} h, NGT vs T2D", n=len(Dm), diff_alignment=o, p=pv))
        rows.append(dict(analysis="myotubes", setting=f"{TP} h, coherence NGT", n=int((lbm == 0).sum()),
                         diff_alignment=R.coherence_loo(Dm[gm == "NGT"]), p=np.nan))
        rows.append(dict(analysis="myotubes", setting=f"{TP} h, coherence T2D", n=int((lbm == 1).sum()),
                         diff_alignment=R.coherence_loo(Dm[gm == "T2D"]), p=np.nan))
        print(f"  {TP} h: coherencia NGT {R.coherence_loo(Dm[gm=='NGT']):+.3f}  T2D {R.coherence_loo(Dm[gm=='T2D']):+.3f}   dif = {o:+.3f}  p = {pv:.4f}")
except FileNotFoundError as ex:
    print("\n(miotubos omitidos:", ex, ")")
pd.DataFrame(rows).to_csv(f"{OUT}/tissue_vs_cell.tsv", sep="\t", index=False)
