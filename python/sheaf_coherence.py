#!/usr/bin/env python3
"""
sheaf_coherence.py
==================
Coherencia cross-tejido de la reorganizacion de redes de coexpresion mediante
un haz celular sobre el grafo de tejidos.

Objeto
------
Para cada estadio s (0 = sano, 1 = intermedio, 2 = T2D) y tejido t se construye
la red de coexpresion W_ts (bicor^beta, beta FIJO por tejido) sobre un conjunto
COMUN de genes, con n IGUAL entre estadios dentro de cada tejido (submuestreo).
El stalk del gen i en el tejido t es su posicion en el embedding espectral de
W_ts (r autovectores del laplaciano normalizado). Los mapas de restriccion son
rotaciones O(r) (Procrustes) que alinean cada embedding a una referencia comun.

Energia del haz por estadio:
    E_s      = sum_{t<u} ||X_ts - X_us||_F^2 / sum_t ||X_ts||_F^2
    E^Delta_s= igual, sobre los CAMBIOS respecto al sano: D_ts = X_ts - X_t0
Score por gen (incoherencia cross-tejido):
    c_i(s)   = sum_t ||D_ts[i] - mean_t D_ts[i]||^2   (varianza entre tejidos)

Nulos
-----
(a) permutacion de etiquetas de estadio dentro de cada tejido (H0: estadios
    intercambiables; controla tamano de muestra porque n se iguala siempre).
(b) permutacion de la correspondencia gen<->gen entre tejidos (H0: no hay
    estructura cross-tejido; controla la topologia de cada red por separado).

Entrada: carpeta export_sheaf/ producida por 01_download_qc_preprocess.R
Salida : tablas TSV listas para el articulo (ver --out).

Uso:
  python sheaf_coherence.py --export_dir export_sheaf --out results/sheaf \
      --tissues GSE76895 GSE18732 GSE27951 --n_genes 800 --reps 30 --B 500
"""
import argparse, json, os, sys, time, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Configuracion: orden fisiologico de los estados por dataset
# ----------------------------------------------------------------------------
STATE_ORDERS = {
    "GSE76895": ["ND", "IGT", "T2D"],
    "GSE18732": ["ND", "IGT", "T2D"],
    "GSE15653": ["Lean", "Obese_noT2D", "Obese_T2D"],
    "GSE27951": ["NGT", "IGT", "T2D"],
}
STAGE_NAMES = ["healthy", "intermediate", "T2D"]

# ----------------------------------------------------------------------------
# Correlacion robusta (biweight midcorrelation, como WGCNA::bicor sin maxPOutliers)
# ----------------------------------------------------------------------------
def bicor(X):
    """X: genes x muestras -> matriz genes x genes."""
    med = np.median(X, axis=1, keepdims=True)
    mad = np.median(np.abs(X - med), axis=1, keepdims=True)
    mad[mad == 0] = 1e-12
    u = (X - med) / (9.0 * mad)
    w = (1 - u ** 2) ** 2 * (np.abs(u) < 1)
    Xw = (X - med) * w
    Xw = Xw / np.sqrt((Xw ** 2).sum(1, keepdims=True) + 1e-300)
    C = Xw @ Xw.T
    C[np.isnan(C)] = 0
    np.fill_diagonal(C, 1.0)
    return np.clip(C, -1, 1)

def adjacency(X, beta, method="bicor"):
    C = bicor(X) if method == "bicor" else np.corrcoef(X)
    C[np.isnan(C)] = 0
    W = np.abs(C) ** beta
    np.fill_diagonal(W, 0)
    return W

def pick_beta(X, powers=range(1, 21), r2_cut=0.80, method="bicor"):
    """Menor beta con ajuste libre de escala R^2 >= r2_cut y pendiente negativa
    (convencion WGCNA). Fallback 6."""
    C = np.abs(bicor(X) if method == "bicor" else np.corrcoef(X))
    np.fill_diagonal(C, 0)
    for b in powers:
        k = (C ** b).sum(1)
        k = k[k > 0.1]
        if k.size < 10:
            continue
        lk = np.log10(k)
        hist, edges = np.histogram(lk, bins=10)
        mids = (edges[:-1] + edges[1:]) / 2
        sel = hist > 0
        if sel.sum() < 3:
            continue
        y = np.log10(hist[sel] / hist.sum())
        A = np.vstack([mids[sel], np.ones(sel.sum())]).T
        coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
        ss_tot = ((y - y.mean()) ** 2).sum()
        r2 = 1 - (res[0] / ss_tot if res.size and ss_tot > 0 else 1)
        if r2 >= r2_cut and coef[0] < 0:
            return int(b)
    return 6

# ----------------------------------------------------------------------------
# Embedding espectral y mapas de restriccion (Procrustes)
# ----------------------------------------------------------------------------
def spectral_embedding(W, r):
    d = W.sum(1) + 1e-12
    Ln = np.eye(W.shape[0]) - W / np.sqrt(np.outer(d, d))
    vals, vecs = np.linalg.eigh(Ln)
    return vecs[:, 1:r + 1]

def procrustes_align(X, R):
    """Rotacion O(r) que minimiza ||X Q - R||_F."""
    U, _, Vt = np.linalg.svd(X.T @ R)
    return X @ (U @ Vt)

def sheaf_energy(stalks):
    """stalks: lista de matrices P x r (una por tejido)."""
    T = len(stalks)
    num = sum(((stalks[t] - stalks[u]) ** 2).sum()
              for t in range(T) for u in range(t + 1, T))
    den = sum((S ** 2).sum() for S in stalks) + 1e-300
    return num / den

def gene_incoherence(stalks):
    A = np.stack(stalks, 0)                  # T x P x r
    return ((A - A.mean(0, keepdims=True)) ** 2).sum((0, 2))

# ----------------------------------------------------------------------------
# Carga de datos
# ----------------------------------------------------------------------------
def load_tissue(export_dir, acc):
    expr = pd.read_csv(os.path.join(export_dir, f"{acc}_expr.tsv"), sep="\t", index_col=0)
    pheno = pd.read_csv(os.path.join(export_dir, f"{acc}_pheno.tsv"), sep="\t")
    pheno = pheno.dropna(subset=["condition"])
    pheno = pheno[pheno[".sample_id"].isin(expr.columns)]
    order = STATE_ORDERS[acc]
    stage = {st: k for k, st in enumerate(order)}
    pheno["stage"] = pheno["condition"].map(stage)
    pheno = pheno.dropna(subset=["stage"])
    pheno["stage"] = pheno["stage"].astype(int)
    return expr, pheno

# ----------------------------------------------------------------------------
# Nucleo del analisis
# ----------------------------------------------------------------------------
def build_stalks(data, genes, betas, rng, r, n_by_tissue, label_perm=False,
                 corr_perm=False):
    """
    Devuelve stalks[s][t] (P x r) alineados a una referencia comun, para cada
    estadio s. Submuestrea n_by_tissue[t] muestras por estadio dentro del tejido.
    label_perm: permuta las etiquetas de estadio dentro de cada tejido (nulo a).
    corr_perm : permuta la correspondencia de genes en tejidos t>=1 (nulo b).
    """
    tissues = list(data)
    T = len(tissues)
    emb = {}
    for t, acc in enumerate(tissues):
        expr, pheno = data[acc]
        X_all = expr.loc[genes].to_numpy()
        stages = pheno["stage"].to_numpy().copy()
        if label_perm:
            stages = rng.permutation(stages)
        for s in range(3):
            idx = np.where(stages == s)[0]
            if idx.size < n_by_tissue[acc]:
                return None
            sel = rng.choice(idx, n_by_tissue[acc], replace=False)
            W = adjacency(X_all[:, sel], betas[acc])
            E = spectral_embedding(W, r)
            if corr_perm and t > 0:
                E = E[rng.permutation(E.shape[0])]
            emb[(s, t)] = E
    # referencia: media (tras alineacion iterativa) de los embeddings sanos
    ref = emb[(0, 0)]
    for _ in range(2):
        aligned = [procrustes_align(emb[(0, t)], ref) for t in range(T)]
        ref = np.mean(aligned, 0)
    stalks = {s: [procrustes_align(emb[(s, t)], ref) for t in range(T)] for s in range(3)}
    return stalks

def analyse_once(data, genes, betas, rng, r, n_by_tissue, **kw):
    st = build_stalks(data, genes, betas, rng, r, n_by_tissue, **kw)
    if st is None:
        return None
    E = np.array([sheaf_energy(st[s]) for s in range(3)])
    delta = {s: [st[s][t] - st[0][t] for t in range(len(st[s]))] for s in (1, 2)}
    E_delta = np.array([np.nan] + [sheaf_energy(delta[s]) for s in (1, 2)])
    gene_c = {s: gene_incoherence(delta[s]) for s in (1, 2)}
    return E, E_delta, gene_c

def run(data, genes, betas, n_by_tissue, r, reps, B, seed, tag, out_dir):
    rng = np.random.default_rng(seed)
    tissues = list(data)
    P = len(genes)
    t0 = time.time()

    # --- observado: media sobre submuestras ---
    Es, Eds, gcs = [], [], {1: [], 2: []}
    for _ in range(reps):
        res = analyse_once(data, genes, betas, rng, r, n_by_tissue)
        if res is None:
            continue
        Es.append(res[0]); Eds.append(res[1])
        for s in (1, 2): gcs[s].append(res[2][s])
    Es, Eds = np.array(Es), np.array(Eds)
    E_obs, Ed_obs = Es.mean(0), np.nanmean(Eds, 0)
    print(f"[{tag}] observado ({len(Es)} submuestras, {time.time()-t0:.0f}s): "
          + " ".join(f"E_{STAGE_NAMES[s]}={E_obs[s]:.4f}" for s in range(3)))

    # --- nulo (a): etiquetas ---
    stat_obs = {"int_max": E_obs[1] - max(E_obs[0], E_obs[2]),
                "T2D_gt_healthy": E_obs[2] - E_obs[0],
                "int_gt_healthy": E_obs[1] - E_obs[0],
                "delta_int_gt_T2D": Ed_obs[1] - Ed_obs[2]}
    null = {k: [] for k in stat_obs}
    gene_null = {1: [], 2: []}
    for b in range(B):
        res = analyse_once(data, genes, betas, rng, r, n_by_tissue, label_perm=True)
        if res is None:
            continue
        E, Ed, gc = res
        null["int_max"].append(E[1] - max(E[0], E[2]))
        null["T2D_gt_healthy"].append(E[2] - E[0])
        null["int_gt_healthy"].append(E[1] - E[0])
        null["delta_int_gt_T2D"].append(Ed[1] - Ed[2])
        for s in (1, 2): gene_null[s].append(gc[s])
    p_lab = {k: (np.sum(np.array(v) >= stat_obs[k]) + 1) / (len(v) + 1) for k, v in null.items()}

    # --- nulo (b): correspondencia de genes ---
    Eb = []
    for b in range(min(B, 200)):
        res = analyse_once(data, genes, betas, rng, r, n_by_tissue, corr_perm=True)
        if res is not None: Eb.append(res[0])
    Eb = np.array(Eb)
    z_corr = (E_obs - Eb.mean(0)) / (Eb.std(0) + 1e-12)
    p_corr = (np.sum(Eb <= E_obs, 0) + 1) / (len(Eb) + 1)   # coherencia > azar => E menor

    # --- tablas ---
    rows = []
    for s in range(3):
        rows.append(dict(analysis=tag, stage=STAGE_NAMES[s], n_tissues=len(tissues),
                         tissues=";".join(tissues), n_genes=P, r=r,
                         n_per_tissue=";".join(f"{a}={n_by_tissue[a]}" for a in tissues),
                         E=E_obs[s], E_sd=Es.std(0)[s],
                         E_delta=Ed_obs[s], E_delta_sd=np.nanstd(Eds, 0)[s],
                         E_null_corr_mean=Eb.mean(0)[s], z_vs_corr_null=z_corr[s],
                         p_vs_corr_null=p_corr[s]))
    stage_df = pd.DataFrame(rows)
    test_df = pd.DataFrame([dict(analysis=tag, statistic=k, observed=stat_obs[k],
                                 null_mean=np.mean(null[k]), null_sd=np.std(null[k]),
                                 p_perm=p_lab[k], B=len(null[k])) for k in stat_obs])

    gene_rows = {"gene": genes}
    for s in (1, 2):
        obs = np.mean(gcs[s], 0)
        nul = np.array(gene_null[s])                 # B x P
        pz = (obs - nul.mean(0)) / (nul.std(0) + 1e-12)
        pp = (np.sum(nul >= obs, 0) + 1) / (nul.shape[0] + 1)
        gene_rows[f"incoherence_{STAGE_NAMES[s]}"] = obs
        gene_rows[f"z_{STAGE_NAMES[s]}"] = pz
        gene_rows[f"p_{STAGE_NAMES[s]}"] = pp
    gene_df = pd.DataFrame(gene_rows).sort_values(f"z_{STAGE_NAMES[1]}", ascending=False)

    os.makedirs(out_dir, exist_ok=True)
    stage_df.to_csv(os.path.join(out_dir, f"{tag}_coherence_by_stage.tsv"), sep="\t", index=False)
    test_df.to_csv(os.path.join(out_dir, f"{tag}_permutation_tests.tsv"), sep="\t", index=False)
    gene_df.to_csv(os.path.join(out_dir, f"{tag}_gene_incoherence.tsv"), sep="\t", index=False)
    print(f"[{tag}] p(int max)={p_lab['int_max']:.3f}  p(T2D>healthy)={p_lab['T2D_gt_healthy']:.3f}  "
          f"z vs nulo-corr={np.round(z_corr, 2)}")
    return stage_df, test_df

# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export_dir", default="export_sheaf")
    ap.add_argument("--out", default="results/sheaf")
    ap.add_argument("--tissues", nargs="+", default=["GSE76895", "GSE18732", "GSE27951"],
                    help="analisis principal; el higado (GSE15653, n=4 en IGT) solo en sensibilidad")
    ap.add_argument("--n_genes", type=int, default=800)
    ap.add_argument("--r", type=int, default=8, help="dimension del embedding espectral")
    ap.add_argument("--beta", type=int, default=None, help="beta fijo global; por defecto se estima por tejido en el estado sano")
    ap.add_argument("--reps", type=int, default=30, help="submuestras a n igual")
    ap.add_argument("--B", type=int, default=500, help="permutaciones nulo (a)")
    ap.add_argument("--min_n", type=int, default=8, help="n minimo por estadio para incluir un tejido")
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--no_loto", action="store_true", help="omitir leave-one-tissue-out")
    ap.add_argument("--no_sens", action="store_true", help="omitir sensibilidad a r (dimension del embedding)")
    args = ap.parse_args()

    # --- cargar ---
    data = {}
    for acc in args.tissues:
        expr, pheno = load_tissue(args.export_dir, acc)
        counts = pheno["stage"].value_counts().reindex(range(3), fill_value=0)
        print(f"{acc}: genes={expr.shape[0]} muestras por estadio={counts.tolist()}")
        if counts.min() < args.min_n:
            print(f"  -> excluido del analisis principal (n_min={counts.min()} < {args.min_n})")
            continue
        data[acc] = (expr, pheno)
    if len(data) < 2:
        sys.exit("Se necesitan al menos 2 tejidos con n suficiente.")

    hv = pd.read_csv(os.path.join(args.export_dir, "high_variance_genes_ordered.tsv"), sep="\t")["gene"].tolist()
    common = set.intersection(*[set(d[0].index) for d in data.values()])
    genes = [g for g in hv if g in common][:args.n_genes]
    print(f"Genes comunes usados: {len(genes)}")

    # n igual por tejido (minimo entre estadios) y beta por tejido (estado sano)
    n_by_tissue, betas = {}, {}
    for acc, (expr, pheno) in data.items():
        n_by_tissue[acc] = int(pheno["stage"].value_counts().min())
        if args.beta is not None:
            betas[acc] = args.beta
        else:
            healthy = pheno.loc[pheno["stage"] == 0, ".sample_id"]
            betas[acc] = pick_beta(expr.loc[genes, healthy].to_numpy())
        print(f"{acc}: n igualado={n_by_tissue[acc]}  beta={betas[acc]}")

    os.makedirs(args.out, exist_ok=True)
    json.dump(dict(vars(args), n_by_tissue=n_by_tissue, betas=betas, genes_used=len(genes)),
              open(os.path.join(args.out, "config.json"), "w"), indent=2)

    # --- principal ---
    run(data, genes, betas, n_by_tissue, args.r, args.reps, args.B, args.seed, "main", args.out)

    # --- sensibilidad: leave-one-tissue-out ---
    if not args.no_loto and len(data) >= 3:
        for drop in list(data):
            sub = {a: v for a, v in data.items() if a != drop}
            run(sub, genes, betas, n_by_tissue, args.r, args.reps, max(args.B // 2, 50),
                args.seed + 1, f"loto_without_{drop}", args.out)

    # --- sensibilidad: dimension del embedding ---
    if not args.no_sens:
        for r_alt in sorted({4, 12} - {args.r}):
            run(data, genes, betas, n_by_tissue, r_alt, args.reps, max(args.B // 4, 50),
                args.seed + 2, f"sens_r{r_alt}", args.out)

    print("Listo. Tablas en", args.out)

if __name__ == "__main__":
    main()
