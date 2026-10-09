#!/usr/bin/env python3
"""
Figuras 2-3 (y 4A): geometria de la respuesta a la insulina dentro de persona.
Datasets: GSE22309 (musculo IS/IR/T2D, GPL91), GSE9105 (sanos 0/30/240 min, GPL96),
GSE7146 (sanos 2 h), GSE231509 (comida mixta H/Ob/T2D, RNA-seq), GSE157988 (prediabetes,
NMN/placebo, RNA-seq), GSE26637 (adiposo IS/IR), Ryden 2016 CAGE (adiposo NO/OB/POB).
Salidas: results/response/*.tsv
"""
import sys, os, re, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..")); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../lib"))
import geo, response as R
OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True); rng = np.random.default_rng(0)
CLOCK = ["DBP", "TEF", "HLF", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "BHLHE40"]
rows = []

def report(name, groups, D, extra=""):
    for g, d in D.items():
        rows.append(dict(dataset=name, group=g, n=len(d), magnitude=np.linalg.norm(d, axis=1).mean(), coherence_loo=R.coherence_loo(d),
                         coherence_naive=R.coherence_naive(d), net=np.linalg.norm(d.mean(0)), note=extra))

# ---------------- GSE22309 musculo ----------------
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); a91 = geo.read_gpl_annot("GPL91.annot.gz")
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"}); p["subj"] = np.arange(len(p)) // 2
p["run"] = p.title.str.extract(r"(Run\d+)")[0].fillna("none")
Z = R.embed(e)
def pairs22309(subjs): return [(p[(p.subj == s) & (p.agent == "untreated")].gsm.iloc[0], p[(p.subj == s) & (p.agent == "insulin")].gsm.iloc[0]) for s in subjs]
D = {g: R.displacements(Z, pairs22309(p[p.grp == g].subj.unique())) for g in ["IS", "IR", "T2D"]}
report("GSE22309_muscle_4h", None, D)
same = p.groupby("subj").run.nunique().eq(1); keep = same[same].index
Dsr = {g: R.displacements(Z, pairs22309(p[(p.grp == g) & (p.subj.isin(keep))].subj.unique())) for g in ["IS", "IR", "T2D"]}
report("GSE22309_muscle_4h_sameRun", None, Dsr, "pares basal/insulina en el mismo lote (analisis primario)")
tests = []
for name, DD in [("full", D), ("sameRun", Dsr)]:
    for g in ["IR", "T2D"]:
        obs, pv = R.perm_test_coherence(DD["IS"], DD[g], rng); tests.append(dict(dataset="GSE22309", subset=name, test=f"coherence IS - {g}", observed=obs, p=pv))
    obs, nm, pv = R.perm_test_direction(DD["IS"], DD["T2D"], rng); tests.append(dict(dataset="GSE22309", subset=name, test="cos(dir IS, dir T2D)", observed=obs, null_mean=nm, p=pv))
# reloj: interaccion grupo x insulina por gen, FDR sobre el set
g22 = geo.probes_to_genes(e, a91); clock_rows = []
for other in ["IR", "T2D"]:
    tmp = []
    for gn in CLOCK:
        if gn not in g22.index: continue
        dA = np.array([g22.loc[gn, b_] - g22.loc[gn, a_] for a_, b_ in pairs22309(p[p.grp == "IS"].subj.unique())])
        dB = np.array([g22.loc[gn, b_] - g22.loc[gn, a_] for a_, b_ in pairs22309(p[p.grp == other].subj.unique())])
        obs, pv = R.interaction_perm(dA, dB, rng); tmp.append(dict(gene=gn, vs=other, resp_IS=dA.mean(), resp_other=dB.mean(), p=pv))
    t = pd.DataFrame(tmp); t["fdr"] = R.bh(t.p); clock_rows.append(t)
pd.concat(clock_rows).to_csv(f"{OUT}/GSE22309_clock_interaction.tsv", sep="\t", index=False)
# programa sano gen a gen y su destino
tabs = {g: R.gene_response_table(g22[g22.mean(axis=1) > 7], pairs22309(p[p.grp == g].subj.unique()))[0] for g in ["IS", "IR", "T2D"]}
prog = tabs["IS"].join(tabs["IR"], rsuffix="_IR").join(tabs["T2D"], rsuffix="_T2D"); prog.to_csv(f"{OUT}/GSE22309_gene_response_by_group.tsv", sep="\t")

# ---------------- GSE9105 sanos, 30 y 240 min ----------------
p9, e9 = geo.read_series_matrix("GSE9105_series_matrix.txt.gz"); a96 = geo.read_gpl_annot("GPL96.annot.gz")
p9["subj"] = p9.title.str.extract(r"^(S\d+)")[0]; p9["tp"] = p9.title.str.extract(r"(\d+min)$")[0]; Z9 = R.embed(e9)
for tp in ["30min", "240min"]:
    prs = [(p9[(p9.subj == s) & (p9.tp == "0min")].gsm.iloc[0], p9[(p9.subj == s) & (p9.tp == tp)].gsm.iloc[0]) for s in p9.subj.dropna().unique() if ((p9.subj == s) & (p9.tp == tp)).any()]
    report("GSE9105_muscle_healthy", None, {tp: R.displacements(Z9, prs)})
g9 = geo.probes_to_genes(e9, a96); g9 = g9[g9.mean(axis=1) > 7]
for tp in ["30min", "240min"]:
    prs = [(p9[(p9.subj == s) & (p9.tp == "0min")].gsm.iloc[0], p9[(p9.subj == s) & (p9.tp == tp)].gsm.iloc[0]) for s in p9.subj.dropna().unique() if ((p9.subj == s) & (p9.tp == tp)).any()]
    R.gene_response_table(g9, prs)[0].to_csv(f"{OUT}/GSE9105_gene_response_{tp}.tsv", sep="\t")
# programa sano replicado (IS 22309 |t|>4 y 9105 240min |t|>3 mismo signo)
t240 = pd.read_csv(f"{OUT}/GSE9105_gene_response_240min.tsv", sep="\t", index_col=0); common = tabs["IS"].index.intersection(t240.index)
core = tabs["IS"].loc[common]; core = core[core.t.abs() > 4]; rep = t240.loc[core.index]
prog55 = core[(np.sign(core["mean"]) == np.sign(rep["mean"])) & (rep.t.abs() > 3)].join(rep, rsuffix="_9105").join(tabs["IR"], rsuffix="_IR").join(tabs["T2D"], rsuffix="_T2D")
prog55.to_csv(f"{OUT}/healthy_program_replicated.tsv", sep="\t")

# ---------------- GSE7146 sanos 2 h ----------------
p7, e7 = geo.read_series_matrix("GSE7146-GPL96_series_matrix.txt.gz"); p7["subj"] = p7.title.str.extract(r"sample (\d+)")[0]; p7["tp"] = np.where(p7.title.str.contains("Pre"), "pre", "post")
Z7 = R.embed(e7); prs = [(p7[(p7.subj == s) & (p7.tp == "pre")].gsm.iloc[0], p7[(p7.subj == s) & (p7.tp == "post")].gsm.iloc[0]) for s in p7.subj.dropna().unique()]
report("GSE7146_muscle_healthy_2h", None, {"healthy": R.displacements(Z7, prs)})

# ---------------- GSE231509 comida mixta ----------------
c = pd.read_csv(geo.path("GSE231509_Read_counts.txt.gz"), sep="\t"); c = c[c.gene_biotype == "protein_coding"].set_index("gene_name").drop(columns=["gene_id", "gene_biotype"]).groupby(level=0).sum()
X = geo.log_cpm(c); cols = pd.Series(X.columns); m = pd.DataFrame({"col": cols, "grp": cols.str.extract(r"^(\w+?)_(?:before|1h)")[0], "tp": np.where(cols.str.contains("before"), "b", "a"), "subj": cols.str.extract(r"_(\d+)$")[0]})
Zm = R.embed(X)
Dm = {g: R.displacements(Zm, [(m[(m.grp == g) & (m.subj == s) & (m.tp == "b")].col.iloc[0], m[(m.grp == g) & (m.subj == s) & (m.tp == "a")].col.iloc[0]) for s in m[m.grp == g].subj.unique()]) for g in m.grp.unique()}
report("GSE231509_muscle_meal_1h", None, Dm)

# ---------------- GSE157988 prediabetes clamp, NMN/placebo ----------------
p15, _ = geo.read_series_matrix("GSE157988_series_matrix.txt.gz"); cc = pd.read_excel(geo.path("GSE157988_NMN_all_gene_counts.xlsx")).set_index("Symbol").groupby(level=0).sum()
Y = geo.log_cpm(cc); cols = pd.Series(Y.columns); mm = pd.DataFrame({"col": cols, "subj": cols.str.extract(r"_(\d+)[AC]_")[0], "tp": cols.str.extract(r"_\d+([AC])_")[0], "st": cols.str.extract(r"_(Basal|Insulin)_")[0]})
ag = p15.assign(subj=p15.title.str.extract(r"(\d+)")[0]).drop_duplicates("subj")[["subj", "agent"]]; mm = mm.merge(ag, on="subj", how="left"); Zp = R.embed(Y)
Dp = {}
for agn in mm.agent.dropna().unique():
    for tp, lab in [("A", "before"), ("C", "after")]:
        prs = [(mm[(mm.subj == s) & (mm.tp == tp) & (mm.st == "Basal")].col.iloc[0], mm[(mm.subj == s) & (mm.tp == tp) & (mm.st == "Insulin")].col.iloc[0]) for s in mm[(mm.agent == agn) & (mm.tp == tp)].subj.unique() if ((mm.subj == s) & (mm.tp == tp) & (mm.st == "Insulin")).any()]
        Dp[f"{agn}_{lab}"] = R.displacements(Zp, prs)
report("GSE157988_muscle_prediabetes_clamp", None, Dp)
for agn in mm.agent.dropna().unique():
    obs, pv = R.paired_swap_test(Dp[f"{agn}_before"], Dp[f"{agn}_after"], R.coherence_loo, rng); tests.append(dict(dataset="GSE157988", subset=agn, test="coherence after - before (paired)", observed=obs, p=pv))

# ---------------- GSE26637 adiposo IS/IR ----------------
p26, e26 = geo.read_series_matrix("GSE26637_series_matrix.txt.gz"); p26["subj"] = p26.title.str.extract(r"(\d+)")[0] + p26["resistance status"]; Z26 = R.embed(e26)
D26 = {}
for st in ["sensitive", "resistant"]:
    prs = [(p26[(p26.subj == s) & (p26.stimulation == "fasting")].gsm.iloc[0], p26[(p26.subj == s) & (p26.stimulation == "hyperinsulinemia")].gsm.iloc[0]) for s in p26[p26["resistance status"] == st].subj.unique()]
    D26[st] = R.displacements(Z26, prs)
report("GSE26637_adipose_clamp", None, D26)

# ---------------- Ryden 2016 CAGE adiposo NO / OB / POB ----------------
try:
    e = pd.read_csv(geo.path("ExpressionTable.txt"), sep="\t", index_col=0); ch = pd.read_csv(geo.path("Cohort.txt"), sep="\t")
    e = e.loc[:, e.columns.isin(ch.Cond)]; ch = ch.set_index("Cond").loc[e.columns]; X = geo.log_cpm(e); Zr = R.embed(X)
    Dr = {}
    for grp, f_, h_ in [("non_obese", "f0", "h0"), ("obese_before", "f0", "h0"), ("obese_2y_after", "f2", "h2")]:
        subs = ch[ch.OB_NO_POB == {"non_obese": "NO", "obese_before": "OB", "obese_2y_after": "POB"}[grp]].Patnr.unique()
        prs = [(f"{s}{f_}", f"{s}{h_}") for s in subs if f"{s}{f_}" in Zr.index and f"{s}{h_}" in Zr.index]; Dr[grp] = R.displacements(Zr, prs)
    report("Ryden2016_adipose_clamp", None, Dr)
    obs, pv = R.perm_test_coherence(Dr["non_obese"], Dr["obese_before"], rng); tests.append(dict(dataset="Ryden2016", subset="full", test="coherence NO - OB", observed=obs, p=pv))
    obs, pv = R.paired_swap_test(Dr["obese_before"], Dr["obese_2y_after"], R.coherence_loo, rng); tests.append(dict(dataset="Ryden2016", subset="paired", test="coherence after - before surgery", observed=obs, p=pv))
except FileNotFoundError as ex: print("Ryden 2016 omitido:", ex)

pd.DataFrame(rows).to_csv(f"{OUT}/response_geometry_by_group.tsv", sep="\t", index=False)
pd.DataFrame(tests).to_csv(f"{OUT}/response_tests.tsv", sep="\t", index=False)
print(pd.DataFrame(rows).round(2).to_string(index=False)); print(pd.DataFrame(tests).round(3).to_string(index=False))
