#!/usr/bin/env python3
"""Modelo jerarquico von Mises-Fisher sobre todas las cohortes pareadas.

Sustituye al meta-analisis de estadisticos de resumen. Cada participante aporta una direccion de
respuesta unitaria; cada cohorte tiene su propia direccion media (porque tejido, estimulo y
plataforma difieren), y la concentracion se modela como log kappa = b0 + b_dataset + b1*deterioro,
de modo que el efecto del deterioro metabolico se estima conjuntamente con todos los participantes
en lugar de combinar un numero por grupo. El parametro de interes es b1: cuanto cae la
concentracion de la respuesta en los grupos con sensibilidad reducida.

Tambien se ajusta la version continua, con el valor M o el estado como covariable numerica cuando
esta disponible, y se compara por razon de verosimilitudes contra el modelo sin efecto.
Salida: results/response/hierarchical_vmf.tsv, pooled_directions.tsv
"""
import os, sys, gzip, numpy as np, pandas as pd, warnings; warnings.filterwarnings("ignore")
from scipy import optimize, special, stats
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../lib")); import geo, response as R
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/response"); os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(0); rows = []

def collect(ds, tissue, stim, D, groups, impaired, extra=None):
    Un = D / np.linalg.norm(D, axis=1, keepdims=True)
    for i, (u, g) in enumerate(zip(Un, groups)):
        r = dict(dataset=ds, tissue=tissue, stimulus=stim, group=g, impaired=int(g in impaired),
                 magnitude=float(np.linalg.norm(D[i])))
        r.update({f"u{j}": float(v) for j, v in enumerate(u)})
        if extra is not None: r.update(extra[i])
        rows.append(r)

def paired(acc, annot, grpmap, basal_key, stim_key, col_agent, col_group, tissue, stim, impaired, subj_col=None):
    p, e = geo.read_series_matrix(f"{acc}_series_matrix.txt.gz")
    if e is None: print("  (omitido)", acc, "sin matriz de expresion"); return
    if annot: e = geo.probes_to_genes(e, geo.read_gpl_annot(annot))
    p["grp"] = p[col_group].map(grpmap)
    if subj_col and subj_col in p.columns: p["subj"] = p[subj_col].astype(str) + "|" + p[col_group].astype(str)
    else:
        key = p[col_group].astype(str)
        p["subj"] = [f"{k}_{i}" for k, i in zip(key, key.groupby(key).cumcount() // 1)]
        p["subj"] = p.groupby([col_group, col_agent]).cumcount().astype(str) + "|" + p[col_group].astype(str)
    rec = []
    for s in p.subj.unique():
        a = p[(p.subj == s) & (p[col_agent] == basal_key)]; b = p[(p.subj == s) & (p[col_agent] == stim_key)]
        if len(a) and len(b) and pd.notna(a.grp.iloc[0]): rec.append((a.gsm.iloc[0], b.gsm.iloc[0], a.grp.iloc[0]))
    if len(rec) < 5: print("  (omitido)", acc, len(rec)); return
    Z = R.embed(e); D = np.array([Z.loc[b].values - Z.loc[a].values for a, b, _ in rec])
    collect(acc, tissue, stim, D, [r[2] for r in rec], impaired); print(f"  {acc}: {len(rec)} participantes")

print("reuniendo direcciones por persona:")
# --- musculo, insulina ---
p, e = geo.read_series_matrix("GSE22309_series_matrix.txt.gz"); e = geo.probes_to_genes(e, geo.read_gpl_annot("GPL91.annot.gz"))
p["grp"] = p.status.map({"insulin sensitive": "IS", "insulin resistant": "IR", "diabetic": "T2D"})
p["subj"] = np.arange(len(p)) // 2; p["run"] = p.title.str.extract(r"(Run\d+)")[0].fillna("none")
rec = []
for s in p.subj.unique():
    a = p[(p.subj == s) & (p.agent == "untreated")]; b = p[(p.subj == s) & (p.agent == "insulin")]
    if len(a) and len(b): rec.append((a.gsm.iloc[0], b.gsm.iloc[0], a.grp.iloc[0], a.run.iloc[0] == b.run.iloc[0]))
Z = R.embed(e); D = np.array([Z.loc[b].values - Z.loc[a].values for a, b, _, _ in rec])
collect("GSE22309", "muscle", "insulin", D, [r[2] for r in rec], {"IR", "T2D"}, [{"same_batch": int(r[3])} for r in rec])
print(f"  GSE22309: {len(rec)} participantes")
# --- adiposo, clamp ---
paired("GSE26637", None, {"sensitive": "sensitive", "resistant": "resistant"},
       "fasting", "hyperinsulinemia", "stimulation", "resistance status", "adipose", "insulin", {"resistant"})
# --- musculo, comida mixta ---
paired("GSE231509", None, {"Metabolically healthy": "healthy", "Obesity": "obese", "T2D, Obesity": "T2D"},
       "before treatment", "1 hour after treatment", "time point", "diagnosis", "muscle", "meal", {"obese", "T2D"},
       subj_col="individual")
# --- ejercicio agudo: permite probar, dentro del modelo, si el efecto es propio de la insulina ---
an = pd.read_csv(f"{RAW}/Human_GRCh38_p13_annot.tsv.gz", sep="\t", usecols=["GeneID", "Symbol"], index_col=0)
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
for acc, tissue, subjcol, tcol, basal, tp in [("GSE202295", "muscle", "from_title", "timepoint", "basal", "recovery"),
                                              ("GSE198922", "adipose", "participant id", "exercise", "pre", "rec")]:
    try:
        q = pheno(acc); cnt = pd.read_csv(f"{RAW}/{acc}_raw_counts_GRCh38_p13_NCBI.tsv.gz", sep="\t", index_col=0)
        cnt.index = an.Symbol.reindex(cnt.index).values; cnt = cnt[pd.notna(cnt.index)].groupby(level=0).sum(); X = geo.log_cpm(cnt)
        q = q[q.gsm.isin(X.columns)].copy()
        q["subj"] = q.title.str.extract(r"((?:T2D|NGT)_\d+)")[0] if subjcol == "from_title" else q[subjcol]
        q["tp"] = q[tcol].str.lower(); Zx = R.embed(X[q.gsm.tolist()]); Zx.index = q.gsm.values
        Dx, gx = [], []
        for grp in ["NGT", "T2D"]:
            w = q[q.diagnosis == grp]
            for s2 in w.subj.dropna().unique():
                a2 = w[(w.subj == s2) & (w.tp == basal)]; b2 = w[(w.subj == s2) & (w.tp == tp)]
                if len(a2) and len(b2): Dx.append(Zx.loc[b2.gsm.iloc[0]].values - Zx.loc[a2.gsm.iloc[0]].values); gx.append(grp)
        if len(Dx) >= 5:
            collect(acc, tissue, "exercise", np.array(Dx), gx, {"T2D"}); print(f"  {acc}: {len(Dx)} participantes")
    except FileNotFoundError as ex: print("  (omitido)", acc, ex)
df = pd.DataFrame(rows); df.to_csv(f"{OUT}/pooled_directions.tsv", sep="\t", index=False)
print(f"\ntotal: {len(df)} participantes, {df.dataset.nunique()} cohortes")

# ---------------- modelo jerarquico ----------------
UCOLS = [c for c in df.columns if c.startswith("u")]
def logC(k, d): return (d / 2 - 1) * np.log(k) - (d / 2) * np.log(2 * np.pi) - np.log(special.ive(d / 2 - 1, k)) - k
def fit(df, with_effect=True, with_interaction=False, joint=True, iters=25):
    """Ajuste conjunto de las direcciones medias por cohorte y de la concentracion.

    log kappa = b0 + b_cohorte + b1*deterioro + b2*(deterioro x ejercicio).

    La direccion media de maxima verosimilitud de una von Mises-Fisher es el vector resultante
    normalizado SOLO si todos los individuos comparten kappa. Aqui kappa varia dentro de cada
    cohorte segun el deterioro, de modo que el estimador conjunto de mu es la media PONDERADA por
    kappa. Se resuelve alternando: dadas las kappa actuales se recalculan las mu ponderadas, y
    dadas las mu se reoptimizan los coeficientes, hasta convergencia. joint=False reproduce el
    estimador en dos pasos (mu fijada al resultante no ponderado) y se conserva para comparacion.
    """
    ds = sorted(df.dataset.unique()); d = len(UCOLS)
    U = df[UCOLS].to_numpy(); imp = df.impaired.to_numpy(float); exo = (df.stimulus == "exercise").to_numpy(float)
    dsv = np.array([ds.index(x) for x in df.dataset])
    mus = np.array([(lambda m: m / np.linalg.norm(m))(U[dsv == i].sum(0)) for i in range(len(ds))])
    def nll(par, mus=None):
        cos = np.einsum("ij,ij->i", U, mus[dsv])
        b0 = par[0]; bd = np.r_[0.0, par[1:len(ds)]]
        b1 = par[len(ds)] if with_effect else 0.0
        b2 = par[len(ds) + 1] if with_interaction else 0.0
        k = np.exp(np.clip(b0 + bd[dsv] + b1 * imp + b2 * imp * exo, -5, 6))
        return -np.sum(logC(k, d) + k * cos)
    x0 = np.r_[np.log(4.0), np.zeros(len(ds) - 1)]
    if with_effect: x0 = np.r_[x0, -0.5]
    if with_interaction: x0 = np.r_[x0, 0.5]
    opt = dict(maxiter=40000, xatol=1e-9, fatol=1e-9)
    r = optimize.minimize(lambda p: nll(p, mus), x0, method="Nelder-Mead", options=opt)
    if joint:
        for _ in range(iters):
            b0 = r.x[0]; bd = np.r_[0.0, r.x[1:len(ds)]]
            b1 = r.x[len(ds)] if with_effect else 0.0
            b2 = r.x[len(ds) + 1] if with_interaction else 0.0
            k = np.exp(np.clip(b0 + bd[dsv] + b1 * imp + b2 * imp * exo, -5, 6))
            mus = np.array([(lambda m: m / np.linalg.norm(m))((U[dsv == i] * k[dsv == i, None]).sum(0)) for i in range(len(ds))])
            r2 = optimize.minimize(lambda p: nll(p, mus), r.x, method="Nelder-Mead", options=opt)
            if abs(r2.fun - r.fun) < 1e-8: r = r2; break
            r = r2
    return r, ds
r0, ds = fit(df, False, False); r1, _ = fit(df, True, False)
has_ex = bool((df.stimulus == "exercise").any())
r2 = fit(df, True, True)[0] if has_ex else None
main = r2 if has_ex else r1
b1 = float(main.x[len(ds)]); b2 = float(main.x[len(ds) + 1]) if has_ex else np.nan
LR1 = 2 * (r0.fun - main.fun); p1 = stats.chi2.sf(LR1, 2 if has_ex else 1)
print(f"\nmodelo jerarquico: {len(df)} participantes, {len(ds)} cohortes, dimension {len(UCOLS)}")
print(f"  b1 (deterioro metabolico, bajo insulina) = {b1:+.3f}   kappa x {np.exp(b1):.2f}")
if has_ex:
    LRi = 2 * (r1.fun - r2.fun); pi_ = stats.chi2.sf(LRi, 1)
    print(f"  b2 (deterioro x ejercicio)               = {b2:+.3f}   LR = {LRi:.2f}  p = {pi_:.4f}")
    print(f"     es decir: bajo ejercicio el efecto del deterioro se cancela ({np.exp(b1 + b2):.2f} veces kappa)")
else: LRi, pi_ = np.nan, np.nan
bs1, bs2 = [], []
for _ in range(int(os.environ.get("T2D_BOOT", 60))):
    take = np.concatenate([rng.choice(np.where(df.dataset == k)[0], int((df.dataset == k).sum()), replace=True) for k in ds])
    try:
        rr = fit(df.iloc[take].reset_index(drop=True), True, has_ex)[0]
        bs1.append(rr.x[len(ds)]);  bs2.append(rr.x[len(ds) + 1] if has_ex else np.nan)
    except Exception: pass
lo, hi = np.percentile(bs1, [2.5, 97.5])
print(f"  IC95% bootstrap de b1: [{lo:+.3f}, {hi:+.3f}]  ->  kappa x [{np.exp(lo):.2f}, {np.exp(hi):.2f}]")
res = [dict(parameter="b1 (metabolic impairment, insulin)", estimate=b1, ci_lo=lo, ci_hi=hi, kappa_ratio=float(np.exp(b1)),
            LR=float(LR1), p=float(p1), n_participants=len(df), n_datasets=len(ds), dim=len(UCOLS))]
if has_ex:
    l2, h2 = np.percentile([x for x in bs2 if np.isfinite(x)], [2.5, 97.5])
    res.append(dict(parameter="b2 (impairment x exercise)", estimate=b2, ci_lo=l2, ci_hi=h2, kappa_ratio=float(np.exp(b2)),
                    LR=float(LRi), p=float(pi_), n_participants=len(df), n_datasets=len(ds), dim=len(UCOLS)))
for i_, k in enumerate(ds):
    lk = main.x[0] + (0.0 if i_ == 0 else main.x[i_])
    res.append(dict(parameter=f"kappa, {k} (reference group)", estimate=float(np.exp(lk)), ci_lo=np.nan, ci_hi=np.nan,
                    kappa_ratio=np.nan, LR=np.nan, p=np.nan, n_participants=int((df.dataset == k).sum()), n_datasets=1, dim=len(UCOLS)))
pd.DataFrame(res).to_csv(f"{OUT}/hierarchical_vmf.tsv", sep="\t", index=False)
print(); print(pd.DataFrame(res).round(4).to_string(index=False))
