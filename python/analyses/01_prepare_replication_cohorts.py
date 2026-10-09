#!/usr/bin/env python3
"""Exporta cohortes de replicacion al formato export_sheaf (expr.tsv, pheno.tsv, high_variance_genes_ordered.tsv)
para landscape.py: GSE164416 (islote, lotes equilibrados), GSE50244 y GSE50398 (islote, HbA1c), GSE25462 (musculo),
GSE64567 y GSE59034 (adiposo), METSIM GSE135134 (adiposo, BMI). Salida: data/export/<acc>/"""
import os, sys, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib")); import geo
E = os.environ.get("T2D_EXPORT", "data/export")
def export(acc, expr, pheno):
    d = f"{E}/{acc}"; os.makedirs(d, exist_ok=True); pheno = pheno[pheno.condition.notna()]; expr = expr[pheno[".sample_id"]]
    expr.rename_axis("gene").to_csv(f"{d}/{acc}_expr.tsv", sep="\t"); pheno.to_csv(f"{d}/{acc}_pheno.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": expr.var(axis=1).sort_values(ascending=False).index}).to_csv(f"{d}/high_variance_genes_ordered.tsv", sep="\t", index=False)
    print(acc, expr.shape, pheno.condition.value_counts().to_dict())
m = geo.ensembl_map_from_gtex()
# GSE164416: donantes vivos; se restringe a lotes con los tres estadios y se centra por lote (batch confundido con estadio)
p, _ = geo.read_series_matrix("GSE164416_series_matrix.txt.gz"); p["sample"] = p.title.str.extract(r"(DP\d+)")[0]
cnt = pd.read_csv(geo.path("GSE164416_DP_htseq_counts.txt.gz"), sep="\t", index_col=0); cnt = cnt[~cnt.index.str.startswith("__")]
p = p[p["diabetes status"].isin(["ND", "IGT", "T2D"]) & p["sample"].isin(cnt.columns)]; X = geo.log_cpm(cnt[p["sample"]]); X.index = [m.get(i.split(".")[0], i) for i in X.index]; X = X[~X.index.str.startswith("ENSG")].groupby(level=0).mean()
keepb = p.library_prep_date.isin(["24/07/2015", "27/03/2014", "30/06/2014"]); p = p[keepb]; X = X[p["sample"]]
for b in p.library_prep_date.unique():
    cols = p.loc[p.library_prep_date == b, "sample"]; X[cols] = X[cols].sub(X[cols].mean(axis=1), axis=0).add(X.mean(axis=1), axis=0)
export("GSE164416", X, pd.DataFrame({".sample_id": p["sample"].values, "condition": p["diabetes status"].values}))
# GSE50244 (RNA-seq) y GSE50398 (array) : mismos donantes, estratos por HbA1c
p, _ = geo.read_series_matrix("GSE50244_series_matrix.txt.gz"); e = pd.read_csv(geo.path("GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz"), sep="\t", index_col=0)
p["sample"] = p.title.str.extract(r"(\d+)")[0]; p = p[p["sample"].isin(e.columns)]; p["hba1c"] = pd.to_numeric(p.hba1c, errors="coerce"); p = p.dropna(subset=["hba1c"])
X = np.log2(e[p["sample"]] * 1000 + 1); X = X[(e[p["sample"]] > 0).mean(axis=1) > 0.8]
cond = pd.cut(p.hba1c, [0, 5.69, 6.49, 99], labels=["ND", "PreD", "T2D"]).astype(str)
export("GSE50244", X, pd.DataFrame({".sample_id": p["sample"].values, "condition": cond.values, "hba1c": p.hba1c.values, "age": p.age.values, "bmi": p.bmi.values, "gender": p.gender.values}))
p, e = geo.read_series_matrix("GSE50398-GPL6244_series_matrix.txt.gz"); p["hba1c"] = pd.to_numeric(p.hba1c, errors="coerce"); p = p.dropna(subset=["hba1c"])
e = e[p.gsm]; e = e[(e > np.percentile(e.values, 30)).mean(axis=1) > 0.5]; cond = pd.cut(p.hba1c, [0, 5.69, 6.49, 99], labels=["ND", "PreD", "T2D"]).astype(str)
export("GSE50398", e, pd.DataFrame({".sample_id": p.gsm.values, "condition": cond.values, "hba1c": p.hba1c.values, "age": p.age.values, "bmi": p.bmi.values, "gender": p.gender.values}))
# GSE25462 musculo
p, e = geo.read_series_matrix("GSE25462_series_matrix.txt.gz"); e = np.log2(e.clip(lower=1)); e = e[(e > np.log2(100)).mean(axis=1) > 0.5]
fh = p["family history"]; cond = np.where(fh == "DM", "T2D", np.where(fh.str.contains("positive"), "ND_FH", "ND"))
export("GSE25462", e, pd.DataFrame({".sample_id": p.gsm.values, "condition": cond, "age": p["age (years)"].values, "sex": p.gender.values, "bmi": p["body mass index (kg/m2)"].values, "hba1c": p["hemoglobin a1c"].values}))
# GSE64567 adiposo (sin T2D por glucosa en ayunas)
p, e = geo.read_series_matrix("GSE64567_series_matrix.txt.gz"); e = e[(e > 7).mean(axis=1) > 0.5]; g_ = pd.to_numeric(p["fasting plasma glucose (mg/dl)"], errors="coerce")
cond = pd.cut(g_, [0, 99.9, 125.9, 999], labels=["NGT", "IFG", "T2D"]).astype(object).where(g_.notna(), None)
export("GSE64567", e, pd.DataFrame({".sample_id": p.gsm.values, "condition": cond.values, "age": p.age.values, "sex": p.gender.values, "bmi": p["bmi (kg/m2)"].values}))
# GSE59034 SAT antes/despues cirugia/nunca obesas
p, e = geo.read_series_matrix("GSE59034_series_matrix.txt.gz"); e = e[(e > 5).mean(axis=1) > 0.5]; st = p["obesity status"]
cond = np.where(st.str.contains("before"), "obese_before", np.where(st.str.contains("after"), "obese_after", "never_obese"))
export("GSE59034", e, pd.DataFrame({".sample_id": p.gsm.values, "condition": cond, "sex": p.gender.values, "title": p.title.values}))
# METSIM
p, _ = geo.read_series_matrix("GSE135134_series_matrix.txt.gz"); e = pd.read_csv(geo.path("GSE135134_METSIM_subcutaneousAdipose_RNAseq_TPMs_n434.txt"), sep="\t", index_col=0)
e.index = [m.get(i.split(".")[0], i) for i in e.index]; e = e[~e.index.str.startswith("ENSG")].groupby(level=0).mean(); X = np.log2(e + 1); X = X[(X > 1).mean(axis=1) > 0.5]; X = X[p.title]
p["bmi"] = pd.to_numeric(p.bmi); cond = pd.cut(p.bmi, [0, 24.99, 29.99, 99], labels=["normal", "overweight", "obese"]).astype(str)
export("METSIM", X, pd.DataFrame({".sample_id": p.title.values, "condition": cond.values, "age": p.age.values, "batch": p["sequencing batch"].values, "tin": p["transcript integrity number"].values, "bmi": p.bmi.values}))
