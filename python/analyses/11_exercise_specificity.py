#!/usr/bin/env python3
"""Prueba de especificidad: ¿pierde el músculo diabético la coordinación de CUALQUIER respuesta,
o solo la de la insulina? Compara la geometría de la respuesta al ejercicio agudo (GSE202295, músculo;
GSE198922, adiposo; biopsias basal/post/recuperación en NGT y T2D) con la de la insulina.
Salida: results/response/exercise_response_geometry.tsv"""
import os, sys, gzip, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import response as R, geo
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True)
an = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot.tsv.gz", sep="\t", index_col=0, usecols=["GeneID", "Symbol"])
def pheno(acc):
    rows = {}
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!Sample_geo_accession"): gsm = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_title"): ti = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
            if l.startswith("!Sample_characteristics_ch1"):
                v = [x.strip('"') for x in l.rstrip().split("\t")[1:]]; k = v[0].split(":")[0].strip()
                rows[k] = [x.split(":", 1)[1].strip() if ":" in x else "" for x in v]
    d = pd.DataFrame(rows); d["gsm"] = gsm; d["title"] = ti; return d
def analyse(acc, tissue, subj, timecol, basal, times):
    p = pheno(acc); cnt = pd.read_csv(f"{RAW}/{acc}_raw_counts_GRCh38_p13_NCBI.tsv.gz", sep="\t", index_col=0)
    cnt.index = an.Symbol.reindex(cnt.index).values; cnt = cnt[pd.notna(cnt.index)].groupby(level=0).sum()
    X = geo.log_cpm(cnt); p = p[p.gsm.isin(X.columns)].copy()
    p["subj"] = p.title.str.extract(r"((?:T2D|NGT)_\d+)")[0] if subj == "from_title" else p[subj]
    p["tp"] = p[timecol].str.lower(); Z = R.embed(X[p.gsm.tolist()]); Z.index = p.gsm.values
    rng = np.random.default_rng(0); out = []
    for grp in ["NGT", "T2D"]:
        q = p[p.diagnosis == grp]
        for t in times:
            prs = [(q[(q.subj == s) & (q.tp == basal)].gsm.iloc[0], q[(q.subj == s) & (q.tp == t)].gsm.iloc[0])
                   for s in q.subj.dropna().unique() if len(q[(q.subj == s) & (q.tp == basal)]) and len(q[(q.subj == s) & (q.tp == t)])]
            if len(prs) < 4: continue
            D = R.displacements(Z, prs)
            out.append(dict(tissue=tissue, dataset=acc, group=grp, timepoint=t, n=len(prs),
                            magnitude=float(np.linalg.norm(D, axis=1).mean()), coherence=R.coherence_loo(D)))
    return pd.DataFrame(out)
m = analyse("GSE202295", "muscle", "from_title", "timepoint", "basal", ["post", "recovery"])
a = analyse("GSE198922", "adipose", "participant id", "exercise", "pre", ["post", "rec"])
r = pd.concat([m, a]); r.to_csv(f"{OUT}/exercise_response_geometry.tsv", sep="\t", index=False); print(r.round(3).to_string(index=False))
