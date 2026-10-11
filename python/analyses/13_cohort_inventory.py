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
 ("GSE50398","Fadista 2014","islet","SuperSeries grouping GSE50244 (RNA-seq) and GSE50397 (expression array); the array accession is GSE50397","none",77,"by HbA1c","landscape (platform replicate)"),
 ("GSE159984","Marselli 2020","islet","cross-sectional + ex vivo paired","palmitate / high glucose / both, 4 d, plus 4 d washout",85,"ND / T2D; paired control per preparation","differential expression; response geometry (ED Fig. 1)"),
 ("GSE18732","Gallagher 2010","muscle","cross-sectional","none",118,"NGT / IGT / T2D","networks, sheaf, landscape"),
 ("GSE25462","Jin 2011","muscle","cross-sectional","none",50,"ND / ND family history / T2D","landscape, differential expression, example network"),
 ("GSE182120","Gabriel 2021","muscle","cross-sectional, clamp-phenotyped","none (resting biopsy)",49,"NGT / T2D, clamp M-value in all","resting association with insulin sensitivity (Fig. 3)"),
 ("GSE182117","Gabriel 2021","myotubes","in vitro time series","high glucose + insulin, 12–54 h",12,"NGT / T2D donors","clock amplitude and treatment effect (Fig. 5g,h)"),
 ("GSE81965","Väremo 2017","myotubes","in vitro, paired within donor","insulin 100 nM, 0.5–2 h",24,"NGT / T2D, balanced for obesity and sex","response geometry in the isolated myocyte (Fig. 7c,d)"),
 ("GSE63887","Väremo 2017","myotubes","in vitro; baseline (0 h) samples of GSE81965, deposited separately","none (0 h)",0,"NGT, non-obese","completes the 24 donors of GSE81965; same study, no additional participants"),
 ("GSE130646","Rubenstein 2020","muscle (single cell)","reference; four preparations of ONE biopsy from a single young male donor","none",1,"one healthy donor","cell-type signature for deconvolution (Fig. 7a); fibre-type markers taken from Table 1 of the same publication, whose data (GSE130977) were not reanalysed"),
 ("GSE22309","Wu 2007","muscle","paired, within person","hyperinsulinaemic clamp, 4 h",55,"insulin-sensitive / insulin-resistant / T2D","response geometry (Figs. 4, 5)"),
 ("GSE9105","Coletta 2008","muscle","paired, within person","hyperinsulinaemic clamp, 30 and 240 min",12,"healthy, family-history negative","response geometry, time course (Fig. 4)"),
 ("GSE7146","Parikh, Groop and Mootha (depositors)","muscle","paired, within person (two platforms, GPL96 used)","insulin infusion, 2 h",6,"healthy","response geometry (Fig. 4)"),
 ("GSE231509","2023","muscle","paired, within person","mixed meal, 1 h",21,"healthy / obese / T2D","response geometry (Fig. 4)"),
 ("GSE157988","Yoshino 2021","muscle","paired, within person, randomised trial","clamp before and after 10 weeks NMN or placebo",23,"prediabetic women, two arms","response geometry, intervention (Figs. 5, 6)"),
 ("GSE202295","Pillon 2022","muscle","paired, within person","acute exercise, post and 3 h recovery",37,"NGT / T2D","specificity (Fig. 5j)"),
 ("GSE224310","Dreher 2023","muscle + adipose","paired, within person; SuperSeries of GSE161749, GSE161750 and GSE208032","acute exercise, before and after training",25,"healthy","characterisation of the measure (Supplementary Note)"),
 ("GSE106800","Wefers 2018","muscle","paired, within person","circadian misalignment, evening and morning",12,"healthy men","characterisation of the measure (Supplementary Note)"),
 ("GSE27951","Keller 2011","adipose","cross-sectional; human arm of a SuperSeries that also contains mouse data","none",33,"lean / obese","networks, sheaf, landscape"),
 ("GSE64567","Winnier 2015","adipose","cross-sectional","none",62,"by fasting glucose","landscape"),
 ("GSE135134","Raulerson 2019 (METSIM)","adipose","cross-sectional, men only","none",434,"by BMI","landscape, differential expression"),
 ("GSE20950","Hardy 2011","adipose","paired depots, within person","none (subcutaneous vs omental)",19,"insulin-sensitive / resistant","individual position in vivo (Fig. 1g)"),
 ("GSE26637","Soronen 2012","adipose","paired, within person","hyperinsulinaemic clamp",10,"insulin-sensitive / resistant","response geometry (Fig. 6b)"),
 ("GSE59034","Arner, Rydén and Mejhert (depositors)","adipose","paired, within person, longitudinal; same group as the 2016 clamp export, donor overlap not excluded","bariatric surgery, 2 and 5 years",16,"obese / never obese","landscape (supplementary)"),
 ("Ryden2016","Rydén 2016","adipose","paired, within person + longitudinal","clamp, before and 2 years after bariatric surgery",46,"non-obese / obese / post-surgery","response geometry, intervention (Fig. 6a)"),
 ("GSE67297","Hanssen 2015","adipose","paired, within person","10 days cold acclimation",7,"T2D","characterisation of the measure (Supplementary Note)"),
 ("GSE66306","Poitou 2015","monocytes","paired, within person, longitudinal","bariatric surgery, 3 months",19,"obese","supplementary"),
 ("GSE156993","no associated publication; cite the accession","PBMC","cross-sectional","none",31,"healthy / T2D","supplementary (no signal)"),
 ("GSE21321","Karolina 2011","blood","cross-sectional (two platforms)","none",25,"healthy / T2D","supplementary (no signal)"),
 ("GTEx v11 (dbGaP phs000424.v11)","GTEx Consortium 2020 (flagship describes v8; v11 is a dbGaP data release)","muscle, adipose (SC and visceral), pancreas, liver, blood and 19 further tissues","cross-sectional, paired organs per donor","none",253,"post-mortem donors","sheaf, individual position, blood-to-tissue (Fig. 1)"),
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
# GSE50398 comparte donantes con GSE50244. GSE63887 ya se cuenta con 0 participantes porque son
# muestras basales de los donantes de GSE81965.
overlap = int(df[df.accession.isin(["GSE50398"])].participants.sum())
n_unique = n_part - overlap
paired = df[df.design.str.contains("paired")]
summary = pd.DataFrame([
    dict(quantity="independent studies", value=n_studies),
    dict(quantity="accessions (GEO series or portal)", value=n_acc),
    dict(quantity="participants, sum over accessions", value=n_part),
    dict(quantity="participants, sum after discounting accessions known to redeposit the same donors (not an identifier-level count of unique individuals; donor identifiers are not public for most cohorts)", value=n_unique),
    dict(quantity="cohorts with paired sampling around a perturbation", value=int((df.perturbation != "none").sum())),
    dict(quantity="participants in paired designs", value=int(paired.participants.sum())),
    dict(quantity="organs represented", value=df.organ.str.split(", ").explode().nunique()),
])
df.to_csv(f"{OUT}/supplementary_table1.tsv", sep="\t", index=False); summary.to_csv(f"{OUT}/inventory_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False)); print()
print(df[["accession","study","organ","perturbation","participants","platform","samples"]].to_string(index=False))
