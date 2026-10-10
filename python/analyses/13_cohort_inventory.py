#!/usr/bin/env python3
"""Inventario de cohortes (Supplementary Table 1). Cuenta por separado estudios independientes,
accessions y participantes unicos, porque un estudio puede depositar varias series (superserie) y
una serie puede contener varias muestras por persona. Lee los series matrix presentes en data/raw
para extraer plataforma y numero de muestras, y completa diseno y uso desde la tabla curada de abajo.
Salida: results/inventory/supplementary_table1.tsv
"""
import os, sys, gzip, glob, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/inventory"); os.makedirs(OUT, exist_ok=True)

# study = publicacion independiente; accessions de la misma publicacion comparten study_id
CUR = [
 # accession, study, organ, design, perturbation, participants, groups, used_in
 ("GSE76895","Solimena 2018","islet","cross-sectional","none",83,"ND / IGT / T2D","networks, sheaf, landscape"),
 ("GSE164416","Wigger 2021","islet","cross-sectional","none",57,"ND / IGT / T2D","landscape, differential expression"),
 ("GSE50244","Fadista 2014","islet","cross-sectional","none",89,"by HbA1c","landscape, differential expression"),
 ("GSE50398","Fadista 2014","islet","cross-sectional (overlapping donors with GSE50244, array; two platforms in one series)","none",77,"by HbA1c","landscape (platform replicate)"),
 ("GSE159984","Marselli 2020","islet","cross-sectional + ex vivo paired","palmitate / high glucose / both, 4 d, plus 4 d washout",85,"ND / T2D; paired control per preparation","differential expression; response geometry (ED Fig. 1)"),
 ("GSE18732","Gallagher 2010","muscle","cross-sectional","none",118,"NGT / IGT / T2D","networks, sheaf, landscape"),
 ("GSE25462","Jin 2011","muscle","cross-sectional","none",50,"ND / ND family history / T2D","landscape, differential expression, example network"),
 ("GSE182120","Gabriel 2021","muscle","cross-sectional, clamp-phenotyped","none (resting biopsy)",49,"NGT / T2D, clamp M-value in all","resting association with insulin sensitivity (Fig. 3)"),
 ("GSE182117","Gabriel 2021","myotubes","in vitro time series","high glucose + insulin, 12–54 h",12,"NGT / T2D donors","clock amplitude and treatment effect (Fig. 5g,h)"),
 ("GSE22309","Wu 2007","muscle","paired, within person","hyperinsulinaemic clamp, 4 h",55,"insulin-sensitive / insulin-resistant / T2D","response geometry (Figs. 4, 5)"),
 ("GSE9105","Coletta 2008","muscle","paired, within person","hyperinsulinaemic clamp, 30 and 240 min",12,"healthy, family-history negative","response geometry, time course (Fig. 4)"),
 ("GSE7146","Rome 2007","muscle","paired, within person (two platforms, GPL96 used)","insulin infusion, 2 h",6,"healthy","response geometry (Fig. 4)"),
 ("GSE231509","2023","muscle","paired, within person","mixed meal, 1 h",21,"healthy / obese / T2D","response geometry (Fig. 4)"),
 ("GSE157988","Yoshino 2021","muscle","paired, within person, randomised trial","clamp before and after 10 weeks NMN or placebo",23,"prediabetic women, two arms","response geometry, intervention (Figs. 5, 6)"),
 ("GSE202295","Pillon 2022","muscle","paired, within person","acute exercise, post and 3 h recovery",37,"NGT / T2D","specificity (Fig. 5j)"),
 ("GSE224310","2023","muscle + adipose","paired, within person","acute exercise, before and after training",25,"healthy","characterisation of the measure (Supplementary Note)"),
 ("GSE106800","Wefers 2018","muscle","paired, within person","circadian misalignment, evening and morning",12,"healthy men","characterisation of the measure (Supplementary Note)"),
 ("GSE27951","2010","adipose","cross-sectional","none",33,"lean / obese","networks, sheaf, landscape"),
 ("GSE64567","2015","adipose","cross-sectional","none",62,"by fasting glucose","landscape"),
 ("GSE135134","METSIM","adipose","cross-sectional","none",434,"by BMI","landscape, differential expression"),
 ("GSE20950","Hardy 2011","adipose","paired depots, within person","none (subcutaneous vs omental)",19,"insulin-sensitive / resistant","individual position in vivo (Fig. 1g)"),
 ("GSE26637","Soronen 2012","adipose","paired, within person","hyperinsulinaemic clamp",10,"insulin-sensitive / resistant","response geometry (Fig. 6b)"),
 ("GSE59034","Mitchell 2015","adipose","paired, within person, longitudinal","bariatric surgery, 2 and 5 years",16,"obese / never obese","landscape (supplementary)"),
 ("Ryden2016","Rydén 2016","adipose","paired, within person + longitudinal","clamp, before and 2 years after bariatric surgery",46,"non-obese / obese / post-surgery","response geometry, intervention (Fig. 6a)"),
 ("GSE67297","Hanssen 2015","adipose","paired, within person","10 days cold acclimation",7,"T2D","characterisation of the measure (Supplementary Note)"),
 ("GSE66306","2015","monocytes","paired, within person, longitudinal","bariatric surgery, 3 months",19,"obese","supplementary"),
 ("GSE156993","2020","PBMC","cross-sectional","none",31,"healthy / T2D","supplementary (no signal)"),
 ("GSE21321","2010","blood","cross-sectional (two platforms)","none",25,"healthy / T2D","supplementary (no signal)"),
 ("GTEx v11","GTEx Consortium","muscle, adipose (SC and visceral), pancreas, liver, blood","cross-sectional, paired organs per donor","none",253,"post-mortem donors","sheaf, individual position, blood-to-tissue (Fig. 1)"),
]
df = pd.DataFrame(CUR, columns=["accession","study","organ","design","perturbation","participants","groups","used_in"])

def probe(acc):
    """plataforma y numero de muestras desde el series matrix, si esta disponible localmente"""
    hits = glob.glob(os.path.join(RAW, f"{acc}*series_matrix*"))
    if not hits: return pd.Series({"platform": "", "samples": np.nan, "series_matrix_found": False})
    plat, n = "", np.nan
    for h in hits:
        op = gzip.open if h.endswith(".gz") else open
        try:
            fh = op(h, "rt", errors="ignore")
        except OSError:
            continue
        with fh:
            for l in fh:
                if l.startswith("!Sample_platform_id") and not plat:
                    plat = sorted(set(x.strip('"') for x in l.rstrip().split("\t")[1:]))[0]
                if l.startswith("!Sample_geo_accession"):
                    k = len(l.rstrip().split("\t")) - 1; n = k if np.isnan(n) else n + k
                if l.startswith("!series_matrix_table_begin"): break
    return pd.Series({"platform": plat, "samples": n, "series_matrix_found": True})
df = pd.concat([df, df.accession.apply(probe)], axis=1)

n_studies = df.study.nunique(); n_acc = df.accession.nunique(); n_part = int(df.participants.sum())
# participantes unicos: GSE50398 comparte donantes con GSE50244, y GSE182117/182120/182121 son un
# superserie del mismo estudio; se descuentan para no contar personas dos veces
overlap = int(df[df.accession.isin(["GSE50398"])].participants.sum())
n_unique = n_part - overlap
paired = df[df.design.str.contains("paired")]
summary = pd.DataFrame([
    dict(quantity="independent studies", value=n_studies),
    dict(quantity="accessions (GEO series or portal)", value=n_acc),
    dict(quantity="participants, sum over accessions", value=n_part),
    dict(quantity="participants, after removing donors shared between accessions of the same study", value=n_unique),
    dict(quantity="cohorts with paired sampling around a perturbation", value=int((df.perturbation != "none").sum())),
    dict(quantity="participants in paired designs", value=int(paired.participants.sum())),
    dict(quantity="organs represented", value=df.organ.str.split(", ").explode().nunique()),
])
df.to_csv(f"{OUT}/supplementary_table1.tsv", sep="\t", index=False); summary.to_csv(f"{OUT}/inventory_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False)); print()
print(df[["accession","study","organ","perturbation","participants","platform","samples"]].to_string(index=False))
