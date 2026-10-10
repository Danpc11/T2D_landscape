#!/usr/bin/env python3
"""Firma de tipos celulares del musculo humano y deconvolucion de las biopsias en volumen.

Construye la firma desde el scRNA-seq de celulas mononucleares de vasto lateral de Rubenstein et al.
2020 (GSE130646; Sci Rep 10:229), agrupando y anotando por marcadores canonicos, y estima con ella
la composicion mononuclear de cada biopsia basal de GSE22309 por minimos cuadrados no negativos.
El paso 19 consume la salida.

Nota sobre el alcance: el scRNA-seq de musculo no captura los mionucleos, porque las fibras son
demasiado grandes para la suspension celular. Las proporciones estimadas son por tanto del
compartimento mononuclear ENTRE SI, no respecto a las fibras. El eje de tipo de fibra se trata
aparte en el paso 19, con los marcadores publicados de ese mismo trabajo.

Entrada: data/raw/GSE130646_RAW.tar (o los cuatro GSM*_Counts.csv.gz ya extraidos)
Salida: results/response/deconvolution_proportions.tsv, celltype_signature.tsv
"""
import os, sys, glob, tarfile, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import optimize
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo, response as R
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True)
WORK = os.environ.get("T2D_WORK", "data/work/GSE130646"); os.makedirs(WORK, exist_ok=True)
SEED = 0
MARKERS = {"FAP": ["PDGFRA", "DCN", "LUM", "COL1A1"],
           "endothelial": ["PECAM1", "VWF", "CDH5", "FLT1"],
           "satellite": ["PAX7", "MYF5", "CHRDL2", "CALCR"],
           "macrophage": ["PTPRC", "CD68", "C1QA", "LYZ"],
           "T/NK": ["CD3E", "IL7R", "NKG7", "CCL5"],
           "pericyte/SMC": ["MYH11", "ACTA2", "RGS5", "NOTCH3"],
           "myofibre": ["MYH1", "MYH2", "MYH7", "ACTA1", "TTN", "DES"],
           "B cell": ["MS4A1", "CD79A"], "lymphatic EC": ["PROX1", "LYVE1", "CCL21"], "neural": ["S100B", "PLP1", "MPZ"]}
try:
    import scanpy as sc
except ImportError:
    sys.exit("scanpy no esta instalado: pip install -r requirements.txt")
sc.settings.verbosity = 0
files = sorted(glob.glob(f"{WORK}/GSM*_Counts.csv.gz"))
if not files:
    tar = f"{RAW}/GSE130646_RAW.tar"
    if not os.path.exists(tar): sys.exit(f"{tar} no existe: descargar segun data/README.md")
    with tarfile.open(tar) as t: t.extractall(WORK)
    files = sorted(glob.glob(f"{WORK}/GSM*_Counts.csv.gz"))
print(f"muestras del atlas: {len(files)}")
ads = []
for f in files:
    X = pd.read_csv(f, index_col=0); a = sc.AnnData(X.T.astype("float32")); a.obs["donor"] = os.path.basename(f).split("_")[2]; ads.append(a)
A = sc.concat(ads, label="batch"); A.var_names_make_unique()
sc.pp.filter_cells(A, min_genes=200); sc.pp.filter_genes(A, min_cells=5)
A.var["mt"] = A.var_names.str.startswith("MT-"); sc.pp.calculate_qc_metrics(A, qc_vars=["mt"], inplace=True, percent_top=None)
A = A[A.obs.pct_counts_mt < 20].copy()
sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A); A.raw = A
sc.pp.highly_variable_genes(A, n_top_genes=2000, batch_key="donor"); Ah = A[:, A.var.highly_variable].copy()
sc.pp.scale(Ah, max_value=10); sc.tl.pca(Ah, n_comps=30, random_state=SEED); sc.pp.neighbors(Ah, n_neighbors=15, random_state=SEED)
sc.tl.leiden(Ah, resolution=0.6, flavor="igraph", n_iterations=2, directed=False, random_state=SEED)
A.obs["leiden"] = Ah.obs.leiden.values
E = pd.DataFrame(A.raw.X.toarray() if hasattr(A.raw.X, "toarray") else np.asarray(A.raw.X), columns=A.raw.var_names, index=A.obs_names)
S = pd.DataFrame({k: E[[g for g in v if g in E.columns]].mean(1) for k, v in MARKERS.items()}); S["leiden"] = A.obs.leiden.values
lab = S.groupby("leiden").mean().idxmax(1); A.obs["celltype"] = A.obs.leiden.map(lab).values
print(f"celulas {A.shape[0]}, clusters {A.obs.leiden.nunique()} -> tipos: {dict(A.obs.celltype.value_counts())}")
ct = A.obs.celltype.values
prof = pd.DataFrame({k: E[ct == k].mean(0) for k in pd.unique(ct)}); prof = prof[prof.max(1) > 0.3]
spec = prof.max(1) / (prof.sum(1) + 1e-9)                      # genes especificos de un tipo
top = pd.concat([prof.loc[spec > 0.5][c].nlargest(80) for c in prof.columns]).index.unique()
sig = prof.loc[top]; sig.to_csv(f"{OUT}/celltype_signature.tsv", sep="\t")
print(f"firma: {sig.shape[0]} genes x {sig.shape[1]} tipos")
# deconvolucion de las biopsias basales de GSE22309
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); e = geo.probes_to_genes(e, geo.read_gpl_annot("GPL91.annot.gz"))
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"}); p["subj"] = np.arange(len(p)) // 2
rec = []
for s in p.subj.unique():
    a = p[(p.subj == s) & (p.agent == "untreated")]; b = p[(p.subj == s) & (p.agent == "insulin")]
    if len(a) and len(b): rec.append((a.gsm.iloc[0], a.grp.iloc[0]))
base = e[[a for a, _ in rec]]; g = [x[1] for x in rec]
shared = [x for x in sig.index if x in base.index]; print(f"genes compartidos con el volumen: {len(shared)}")
M = sig.loc[shared].to_numpy(); B = base.loc[shared].to_numpy()
P = pd.DataFrame([(lambda w: w / (w.sum() + 1e-12))(optimize.nnls(M, B[:, j])[0]) for j in range(B.shape[1])], columns=sig.columns)
P["group"] = g; P.to_csv(f"{OUT}/deconvolution_proportions.tsv", sep="\t", index=False)
print(); print(P.groupby("group").mean().round(3).to_string())
