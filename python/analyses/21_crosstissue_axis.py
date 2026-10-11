#!/usr/bin/env python3
"""Eje individual compartido entre musculo y tejido adiposo: concordancia, metaanalisis y redes.

Analisis principal del articulo. Pregunta que genes situan a una persona de forma consistente
respecto a las demas en SUS DOS tejidos a la vez, cantidad distinta de la correlacion entre genes
dentro de un tejido.

Pasos:
  1. Concordancia entre tejidos: para cada gen, correlacion de Spearman ENTRE INDIVIDUOS de su
     expresion en musculo y en adiposo. El fondo es la media de esa cantidad sobre genes al azar.
  2. Significacion por conjunto frente a conjuntos aleatorios del mismo tamano y, donde se indica,
     equiparados ademas en expresion media y desviacion tipica (los 80 genes mas proximos en ese
     espacio bidimensional).
  3. Intervalos por bootstrap SOBRE INDIVIDUOS, no sobre genes: el individuo es la unidad de
     replicacion.
  4. Metaanalisis de efectos aleatorios (DerSimonian-Laird) sobre la z de Fisher, con I2 y Q.
  5. Dos redes: una sobre la concordancia entre tejidos, otra DENTRO de cada tejido con el mismo
     exponente y los mismos participantes, de modo que la conectividad sea comparable.

Los conjuntos de genes proceden de la literatura y NUNCA se derivan de los datos que se ponen a
prueba; se definen aqui de forma explicita para que queden congelados.

RESULTADO DE LOS DOS NULOS (AAGMEx, n = 232). El control equiparado es notablemente mas severo:

  conjunto          r      p (tamano)   p (expresion y variabilidad)
  Genetic        +0.808      0.001              0.002
  Clock output   +0.269      0.001              0.014
  Glucocorticoid +0.168      0.001              0.100
  Insulin core   +0.190      0.007              0.116
  Oscillator     +0.127      0.055              0.260
  Oxidative      +0.030      0.533              0.810

Solo el componente genetico y la salida del reloj superan ambos controles. Los glucocorticoides y
el nucleo de insulina superan el de tamano pero no el de expresion. Caben dos lecturas y los datos
no permiten elegir: que parte de su concordancia se deba a ser genes muy expresados y variables,
que se miden mejor y por tanto correlacionan mejor entre tejidos; o que el control sea demasiado
conservador, porque los genes equiparados a FKBP5 en expresion y variabilidad son ellos mismos
candidatos a estar bajo influencia hormonal. Debe reportarse el par de valores, no solo el primero.

Entradas: GSE95674 (adiposo) y GSE95675 (musculo) de AAGMEx; GSE41168 como replicacion.
Salidas: results/gtex/live_crosstissue_r.tsv, ci_*.tsv, meta_random_effects.tsv,
         crosstissue_modules.tsv, connectivity_AAGMEx.tsv
"""
import os, sys, gzip, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
RAW = os.environ.get("T2D_RAW", "data/raw"); OUT = os.environ.get("T2D_OUT", "results/gtex"); os.makedirs(OUT, exist_ok=True)
NBOOT = int(os.environ.get("T2D_BOOT", "400")); NPERM = int(os.environ.get("T2D_PERM", "3000"))

SETS = {
    "Genetic": ["GSTM1", "GSTT1", "ERAP2"],
    "Clock output": ["DBP", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART"],
    "Glucocorticoid": ["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "DDIT4", "SESN1", "GLUL", "SGK1", "CEBPD"],
    "Insulin core": ["TXNIP", "PDK4", "HK2", "IRS2", "ANGPTL4"],
    "Oscillator": ["ARNTL", "CLOCK", "CRY1", "CRY2"],
    "Oxidative": ["PPARGC1A", "NDUFB8", "COX5A", "SDHB", "UQCRC1"],
}

def read_series(acc):
    meta, tab, intab = {}, [], False
    with gzip.open(f"{RAW}/{acc}_series_matrix.txt.gz", "rt", errors="ignore") as fh:
        for l in fh:
            if l.startswith("!series_matrix_table_begin"): intab = True; continue
            if l.startswith("!series_matrix_table_end"): break
            if intab: tab.append(l.rstrip("\n").split("\t")); continue
            for k, t in [("gsm", "!Sample_geo_accession"), ("title", "!Sample_title")]:
                if l.startswith(t) and k not in meta: meta[k] = [x.strip('"') for x in l.rstrip().split("\t")[1:]]
    e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
    e.index = [i.strip('"') for i in e.index]
    return pd.DataFrame(meta), e.apply(pd.to_numeric, errors="coerce").dropna()

def paired_aagmex():
    pa, ea = read_series("GSE95674"); pm, em = read_series("GSE95675")
    pa["s"] = pa.title.str.replace("-M", "", regex=False); pm["s"] = pm.title.str.replace("-M", "", regex=False)
    ids = sorted(set(pa.s) & set(pm.s))
    A = np.log2(ea[[pa.set_index("s").gsm[s] for s in ids]].clip(lower=1)); A.columns = ids
    M = np.log2(em[[pm.set_index("s").gsm[s] for s in ids]].clip(lower=1)); M.columns = ids
    g = sorted(set(A.index) & set(M.index))
    return M.loc[g], A.loc[g], ids

def concordance(M, A):
    return pd.Series([stats.spearmanr(M.loc[x], A.loc[x])[0] for x in M.index], index=M.index).dropna()

def matched_null(target, mu, sd, rng, pool):
    """Conjunto nulo equiparado en expresion media y variabilidad, no solo en tamano."""
    out = []
    for g in target:
        d = np.abs(mu - mu[g]) / mu.std() + np.abs(sd - sd[g]) / sd.std()
        out.append(rng.choice(d.sort_values().index[1:81]))
    return out

def boot_ci(M, A, genes, ids, rng, nboot=NBOOT):
    n = len(ids); vals = []
    for _ in range(nboot):
        i = rng.integers(0, n, n)
        vals.append(np.mean([stats.spearmanr(M.loc[x].values[i], A.loc[x].values[i])[0] for x in genes]))
    return np.percentile(vals, 2.5), np.percentile(vals, 97.5), float(np.std(vals))

def dersimonian_laird(r, se):
    """Efectos aleatorios sobre z de Fisher; devuelve r combinada, IC, tau2, I2, Q."""
    z = np.arctanh(r); sez = se / (1 - r ** 2)
    w = 1 / sez ** 2; zbar = np.sum(w * z) / np.sum(w)
    Q = np.sum(w * (z - zbar) ** 2); k = len(z)
    Cc = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    tau2 = max(0.0, (Q - (k - 1)) / Cc) if Cc > 0 else 0.0
    wr = 1 / (sez ** 2 + tau2); zm = np.sum(wr * z) / np.sum(wr); sem = np.sqrt(1 / np.sum(wr))
    I2 = max(0.0, (Q - (k - 1)) / Q * 100) if Q > 0 else 0.0
    return np.tanh(zm), np.tanh(zm - 1.96 * sem), np.tanh(zm + 1.96 * sem), tau2, I2, Q, k

def crosstissue_network(M, A, con, thr=0.30, k=8):
    """Red SOBRE la concordancia entre tejidos: agrupa genes que ademas covarian entre individuos."""
    hi = con[con > thr].index.tolist()
    Zm = ((M.loc[hi].T - M.loc[hi].mean(1)) / (M.loc[hi].std(1) + 1e-9)).T
    Za = ((A.loc[hi].T - A.loc[hi].mean(1)) / (A.loc[hi].std(1) + 1e-9)).T
    Z = (Zm + Za) / 2
    Cm = np.corrcoef(Z.values); D = 1 - np.abs(Cm); np.fill_diagonal(D, 0); D = (D + D.T) / 2
    lab = fcluster(linkage(squareform(D, checks=False), method="average"), t=k, criterion="maxclust")
    return pd.Series(lab, index=hi)

def within_tissue_connectivity(M, A, beta=6, ntop=4000):
    """Conectividad dentro de cada tejido con el MISMO exponente y los MISMOS participantes."""
    v = (A.std(1).rank(pct=True) + M.std(1).rank(pct=True)) / 2
    sel = list(v.sort_values(ascending=False).index[:ntop])
    for s in set(sum(SETS.values(), [])):
        if s in A.index and s not in sel: sel.append(s)
    out = {}
    for nm, E in [("adipose", A), ("muscle", M)]:
        X = E.loc[sel]; Cm = np.corrcoef((X.T - X.mean(1)).T.values)
        Adj = np.abs(Cm) ** beta; np.fill_diagonal(Adj, 0)
        out[nm] = pd.Series(Adj.sum(1), index=sel)
    return pd.DataFrame(out)

if __name__ == "__main__":
    rng = np.random.default_rng(0)
    M, A, ids = paired_aagmex()
    print(f"AAGMEx: {len(ids)} individuos con ambos tejidos, {len(M)} genes comunes")
    con = concordance(M, A); con.to_frame("r").rename_axis("probe").to_csv(f"{OUT}/live_crosstissue_r.tsv", sep="\t")
    bg = con.sample(min(3000, len(con)), random_state=0).mean()
    print(f"fondo: {bg:+.3f} | genes con r>0.5: {(con > 0.5).mean()*100:.1f}% | >0.7: {(con > 0.7).mean()*100:.2f}%")
    mu, sd = A.mean(1), A.std(1)
    rows = []
    for nm, gs in SETS.items():
        g = [x for x in gs if x in con.index]
        if len(g) < 3: continue
        obs = con[g].mean()
        null_sz = [con.sample(len(g), random_state=int(rng.integers(1e9))).mean() for _ in range(NPERM)]
        null_mt = [con[matched_null(g, mu, sd, rng, con.index)].mean() for _ in range(min(NPERM, 500))]
        lo, hi, se = boot_ci(M, A, g, ids, rng)
        rows.append(dict(cohort="AAGMEx", set=nm, n=len(ids), n_genes=len(g), r=obs, lo=lo, hi=hi, se=se,
                         p_size=(np.sum(np.array(null_sz) >= obs) + 1) / (NPERM + 1),
                         p_matched=(np.sum(np.array(null_mt) >= obs) + 1) / (len(null_mt) + 1)))
        print(f"  {nm:16s} r={obs:+.3f} [{lo:+.3f}, {hi:+.3f}]  p={rows[-1]['p_size']:.4f} (matched {rows[-1]['p_matched']:.4f})")
    ci = pd.DataFrame(rows); ci.to_csv(f"{OUT}/ci_AAGMEx.tsv", sep="\t", index=False)
    mods = crosstissue_network(M, A, con)
    mods.to_frame("module").rename_axis("gene").to_csv(f"{OUT}/crosstissue_modules.tsv", sep="\t")
    print(f"\nred de concordancia: {mods.nunique()} modulos sobre {len(mods)} genes")
    K = within_tissue_connectivity(M, A); K.to_csv(f"{OUT}/connectivity_AAGMEx.tsv", sep="\t")
    print("conectividad dentro de tejido (percentil medio por conjunto):")
    for nm, gs in SETS.items():
        g = [x for x in gs if x in K.index]
        if len(g) < 3: continue
        pa_, pm_ = K.adipose.rank(pct=True)[g].mean(), K.muscle.rank(pct=True)[g].mean()
        print(f"  {nm:16s} adiposo {pa_:.3f} | musculo {pm_:.3f}")
