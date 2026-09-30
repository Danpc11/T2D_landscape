#!/usr/bin/env python3
"""
landscape.py
============
Paisaje cuasi-potencial (U, J) de la progresion sano -> intermedio -> T2D por
tejido, a partir de cohortes transversales.

Objeto central: U(x) = -ln P(x) sobre un embedding de baja dimension del
transcriptoma, condicionado a covariables (BMI, edad, sexo) si estan en pheno.

ETAPA 1 (paisaje):
  - Potencial condicional: residuos tras regresar covariables, PCA a r dims.
  - Densidad por KDE con barrido de ancho de banda; U = -ln p.
  - Puntos criticos por SCORE (grad ln p) y su jacobiano: minimos vs sillas.
  - Persistencia de subnivel (H0) de U en malla 2D: numero de cuencas
    persistentes y barreras = muerte del par minimo-silla, robustas a escala.
  - GMM: dBIC(2 comp - 1 comp); bimodalidad.
  - Asimetria de barreras dU(sano->T2D) vs dU(T2D->sano), IC bootstrap.
  - Indice de transicion critica I_c (Mojtahedi/Huang 2016) por estadio con
    n igualado y nulo por permutacion.
  - Informacion de Fisher a lo largo de una covariable continua (HbA1c/glucosa)
    si existe: g(theta) ~ 2*KL(P_theta || P_theta+d)/d^2 por ventanas.

ETAPA 2 (flujo):
  - Puente de Schrodinger (OT entropico, Sinkhorn) entre P_sano y P_T2D en el
    embedding: plan pi (matriz de atencion doblemente estocastica) -> deriva
    b(x_i) = sum_j pi_ij x_j / a_i - x_i.
  - Descomposicion de Hodge de la deriva sobre el grafo kNN de pacientes:
    gradiente (= -grad U dinamico) + residual (rotacional+armonico = J).
    ratio_J = ||J|| / ||grad U||  (0 = equilibrio / dinamica de gradiente).
  - Comprobacion: la ruta del puente pasa por la region de los pacientes IGT.

Entrada: export_sheaf/ (01) ; columnas opcionales en *_pheno.tsv: bmi, age,
sex, hba1c, glucose (nombres insensibles a mayusculas).
Salida: results/landscape/<acc>_*.tsv y landscape_summary.tsv

Uso:
  python landscape.py --export_dir export_sheaf --out results/landscape \
      --n_genes 800 --r 3 --reps 20 --B 200
"""
import argparse, json, os, sys, warnings
import numpy as np, pandas as pd
from scipy import linalg
from scipy.spatial.distance import cdist
warnings.filterwarnings("ignore", category=RuntimeWarning)

STATE_ORDERS = {"GSE76895": ["ND", "IGT", "T2D"], "GSE18732": ["ND", "IGT", "T2D"],
                "GSE15653": ["Lean", "Obese_noT2D", "Obese_T2D"], "GSE27951": ["NGT", "IGT", "T2D"]}
STAGE = ["healthy", "intermediate", "T2D"]
NO_COVAR = False
NO_BALANCE = False
COVARS = {"bmi": ["bmi", "body mass index", "body_mass_index", "bodymassindex", "body mass"],
          "age": ["age", "edad", "years"], "sex": ["sex", "gender"]}
# grupos de controles continuos: se calcula Fisher para cada grupo disponible
CONTROL_GROUPS = {"hba1c":   ["hba1c", "hb a1c", "a1c", "glycated", "hemoglobin a1c"],
                  "glucose": ["fasting glucose", "fasting plasma glucose", "fpg", "glucose 0h", "glucose", "glucemia"]}
CONTROLS = sum(CONTROL_GROUPS.values(), [])
CONTROL_OVERRIDE = None   # --control nombre_de_columna

def expand_characteristics(pheno):
    """GEO guarda covariables como texto 'clave: valor' en characteristics_ch1.N.
    Expande cada par en una columna 'clave' (si no existe ya) para que find_col
    las encuentre por nombre. Tambien parte celdas con varios pares ' | '."""
    pheno = pheno.copy(); new = {}
    for c in list(pheno.columns):
        if pd.api.types.is_numeric_dtype(pheno[c]): continue
        vals = pheno[c].astype(str)
        if not (vals.str.contains(r"^[^:|]{1,40}:\s*\S", regex=True).mean() > 0.5): continue
        for i, cell in vals.items():
            for part in str(cell).split(" | "):
                if ":" not in part: continue
                k, v = part.split(":", 1); k = k.strip().lower().replace(" ", "_"); v = v.strip()
                if not k or k in pheno.columns: continue
                new.setdefault(k, {})[i] = v
    for k, d in new.items():
        pheno[k] = pd.Series(d, dtype=object).reindex(pheno.index)
    return pheno

# ----------------------------------------------------------------------------
# utilidades
# ----------------------------------------------------------------------------
def find_col(pheno, names):
    low = {c.lower().replace("_", " "): c for c in pheno.columns}
    for n in names:                       # 1) nombre exacto
        if n in low: return low[n]
    for n in names:                       # 2) contenido, como palabra (evita 'age' en 'percentage')
        for lc, c in low.items():
            import re
            if re.search(r"(^|[^a-z])" + re.escape(n) + r"([^a-z]|$)", lc): return c
    return None

def residualize(Y, pheno, stage=None):
    """Y: muestras x genes. Regresa covariables disponibles (bmi, age, sex).
    [within-stage] Las covariables se centran DENTRO de cada estadio: asi se
    elimina su efecto entre pacientes del mismo grupo sin borrar la diferencia
    entre grupos (BMI es colineal con el estadio; regresarlo crudo destruye
    exactamente la senal que se quiere estudiar)."""
    cols, X = [], [np.ones(Y.shape[0])]
    for key, names in COVARS.items():
        c = find_col(pheno, names)
        if c is None: continue
        v = pheno[c]
        if key == "sex": v = pd.Series(pd.factorize(v.astype(str).str.strip().str.lower())[0], index=v.index).astype(float)
        v = pd.to_numeric(v.astype(str).str.replace(",", "."), errors="coerce").to_numpy() if key != "sex" else v.to_numpy()
        if np.isnan(v).mean() > 0.2: continue
        v = np.where(np.isnan(v), np.nanmean(v), v)
        if stage is not None:
            for g in np.unique(stage): v[stage == g] -= v[stage == g].mean()
        X.append(v); cols.append(key)
    X = np.column_stack(X)
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X[:, 1:] @ beta[1:] if len(cols) else Y, cols

def pca(Y, r):
    Yc = Y - Y.mean(0)
    U, S, Vt = np.linalg.svd(Yc, full_matrices=False)
    Z = U[:, :r] * S[:r]
    return Z, (S[:r] ** 2 / (S ** 2).sum())

# ----------------------------------------------------------------------------
# densidad, score, puntos criticos
# ----------------------------------------------------------------------------
def kde(Z, X, h, w=None):
    """KDE gaussiano de los datos Z evaluado en X. w: pesos por muestra (suman 1).
    [FIX-design] Sin pesos, la densidad refleja cuantos pacientes se reclutaron por
    estadio (diseno del estudio), no la frecuencia de los estados. Por defecto el
    pipeline pondera cada estadio a masa igual (1/3)."""
    D2 = cdist(X, Z, "sqeuclidean")
    K = np.exp(-D2 / (2 * h * h))
    if w is None: w = np.full(Z.shape[0], 1.0 / Z.shape[0])
    return (K @ w) / (2 * np.pi * h * h) ** (Z.shape[1] / 2) + 1e-300, K * w[None, :]

def score(Z, X, h, w=None):
    p, Kw = kde(Z, X, h, w)
    num = (Kw[:, :, None] * (Z[None, :, :] - X[:, None, :])).sum(1) / (h * h)
    return num / (Kw.sum(1)[:, None] + 1e-300)      # grad ln p

def stage_weights(st, balance=True):
    if not balance: return np.full(len(st), 1.0 / len(st))
    w = np.zeros(len(st)); S = np.unique(st)
    for g in S: w[st == g] = 1.0 / (len(S) * (st == g).sum())
    return w

def silverman(Z):
    n, d = Z.shape
    return (4 / (d + 2)) ** (1 / (d + 4)) * n ** (-1 / (d + 4)) * Z.std(0).mean()

def critical_points(Z, h, n_starts=200, seed=0, w=None):
    """Ascenso por el score desde puntos aleatorios -> modos; jacobiano numerico
    para clasificar (minimo de U = modo de p). Sillas: puntos donde el score se
    anula con jacobiano indefinido, buscados como minimos de |score| en la
    trayectoria entre modos."""
    rng = np.random.default_rng(seed)
    lo, hi = Z.min(0), Z.max(0)
    X = rng.uniform(lo, hi, size=(n_starts, Z.shape[1]))
    for _ in range(300):
        X = X + 0.1 * h * h * score(Z, X, h, w)
    # agrupar modos
    modes = []
    for x in X:
        if not any(np.linalg.norm(x - m) < 1.0 * h for m in modes): modes.append(x)
    modes = np.array(modes)
    dens, _ = kde(Z, modes, h, w)
    keep = dens > 0.05 * dens.max()          # descartar modos espurios de baja densidad
    return modes[keep], -np.log(dens[keep])

# ----------------------------------------------------------------------------
# persistencia de subnivel H0 sobre malla (2 primeras dims)
# ----------------------------------------------------------------------------
def sublevel_persistence_2d(Z, h, grid=80, w=None):
    z = Z[:, :2]
    lo, hi = z.min(0) - 2 * h, z.max(0) + 2 * h
    gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], grid), np.linspace(lo[1], hi[1], grid))
    G = np.column_stack([gx.ravel(), gy.ravel()])
    p, _ = kde(z, G, h, w)
    U = -np.log(p).reshape(grid, grid)
    order = np.argsort(U.ravel())
    parent = -np.ones(U.size, int); birth = {}; pairs = []
    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for idx in order:
        i, j = divmod(idx, grid)
        nbrs = [(i + di) * grid + (j + dj) for di in (-1, 0, 1) for dj in (-1, 0, 1)
                if (di or dj) and 0 <= i + di < grid and 0 <= j + dj < grid and parent[(i + di) * grid + (j + dj)] >= 0]
        parent[idx] = idx
        roots = {find(n) for n in nbrs}
        if not roots:
            birth[idx] = U.ravel()[idx]
        else:
            roots = sorted(roots, key=lambda r: birth[r])
            elder = roots[0]
            parent[idx] = elder
            for r in roots[1:]:
                pairs.append((birth[r], U.ravel()[idx], r))     # (nacimiento, muerte=silla)
                parent[r] = elder
    pairs.sort(key=lambda t: -(t[1] - t[0]))
    # (nacimiento, muerte, persistencia, coordenadas del minimo joven)
    pers = [(b, d, d - b, G[r]) for b, d, r in pairs]
    return U, (gx, gy), pers, birth

def persistence_threshold(Z2, h, rng, B=30, alpha=0.90, w=None):
    """[Stability theorem] d_B(Dgm(U_hat), Dgm(U)) <= ||U_hat - U||_inf. Estimamos
    eps = cuantil alpha de ||U_hat_b - U_hat||_inf sobre bootstrap del KDE (Fasy et
    al. 2014) restringido al soporte de los datos. Toda barra con persistencia
    > 2 eps corresponde a una barra real con confianza alpha. Umbral final:
    max(ln 2, 2 eps)  (interpretabilidad + garantia estadistica)."""
    n = len(Z2); base = -np.log(kde(Z2, Z2, h, w)[0]); sup = []
    if w is None: w = np.full(n, 1.0 / n)
    # sup restringido a la region de alta densidad (80% de puntos mas densos): es
    # donde estan minimos y sillas; en las colas log p es inestable y no informa.
    core = base <= np.quantile(base, 0.8)
    for _ in range(B):
        idx = rng.choice(n, n, p=w); sup.append(np.abs(-np.log(kde(Z2[idx], Z2[core], h)[0]) - base[core]).max())
    eps = np.quantile(sup, alpha)
    return max(np.log(2), 2 * eps), eps

# ----------------------------------------------------------------------------
# GMM BIC
# ----------------------------------------------------------------------------
def gmm_dbic(Z, seed=0):
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        return np.nan, np.nan
    b = [GaussianMixture(k, covariance_type="full", n_init=5, random_state=seed).fit(Z).bic(Z) for k in (1, 2, 3)]
    return b[0] - b[1], b[0] - b[2]           # >0 favorece 2 (3) componentes

# ----------------------------------------------------------------------------
# indice de transicion critica (Mojtahedi et al. 2016)
# ----------------------------------------------------------------------------
def critical_index(Yg):
    """Yg: muestras x genes (un estadio). I_c = <|corr genes|> / <|corr muestras|>"""
    Cg = np.corrcoef(Yg.T); Cs = np.corrcoef(Yg)
    ug = np.abs(Cg[np.triu_indices_from(Cg, 1)]); us = np.abs(Cs[np.triu_indices_from(Cs, 1)])
    return np.nanmean(ug) / (np.nanmean(us) + 1e-12)

# ----------------------------------------------------------------------------
# Fisher a lo largo de covariable continua
# ----------------------------------------------------------------------------
def fisher_along(Z, theta, h, n_win=6, min_dtheta_frac=0.08):
    """g(theta) ~ 2 KL(p_w || p_w') / (dtheta)^2 entre ventanas deslizantes.
    [FIX-edge] Se exige dtheta >= min_dtheta_frac * rango(theta): en zonas donde
    muchos pacientes comparten casi el mismo theta (extremo sano), dtheta -> 0 y
    g explotaba artificialmente. Las ventanas se separan hasta cumplir el minimo."""
    o = np.argsort(theta); Z, theta = Z[o], theta[o]
    n = len(theta); w = max(8, n // n_win); rng_th = theta[-1] - theta[0]
    if rng_th <= 0: return np.array([]), np.array([])
    centers, g = [], []
    for s0 in range(0, n - w, max(1, w // 2)):
        A = Z[s0:s0 + w]; thA = theta[s0:s0 + w].mean()
        s1 = s0 + w
        while s1 + w <= n and theta[s1:s1 + w].mean() - thA < min_dtheta_frac * rng_th: s1 += 1
        if s1 + w > n: break
        B = Z[s1:s1 + w]; dth = theta[s1:s1 + w].mean() - thA
        pa, _ = kde(A, A, h); pb, _ = kde(B, A, h)
        kl = np.mean(np.log(pa) - np.log(pb))
        centers.append((thA + theta[s1:s1 + w].mean()) / 2); g.append(2 * max(kl, 0) / dth ** 2)
    return np.array(centers), np.array(g)

# ----------------------------------------------------------------------------
# ETAPA 2: puente de Schrodinger (Sinkhorn) + Hodge de la deriva
# ----------------------------------------------------------------------------
def sinkhorn(A, B, eps, iters=500):
    C = cdist(A, B, "sqeuclidean"); K = np.exp(-C / eps)
    a, b = np.ones(len(A)) / len(A), np.ones(len(B)) / len(B)
    u, v = np.ones(len(A)), np.ones(len(B))
    for _ in range(iters):
        u = a / (K @ v + 1e-300); v = b / (K.T @ u + 1e-300)
    return u[:, None] * K * v[None, :]         # plan pi (atencion doblemente estocastica)

def knn_graph(X, k):
    D = cdist(X, X); np.fill_diagonal(D, np.inf)
    edges = set()
    for i in range(len(X)):
        for j in np.argsort(D[i])[:k]: edges.add((min(i, j), max(i, j)))
    return np.array(sorted(edges))

def hodge_decompose_drift(X, drift, k=8):
    """Proyecta la deriva sobre las aristas del grafo kNN (flujo f_ij), ajusta
    potencial phi por minimos cuadrados (f ~ B1 phi) y separa gradiente/residual."""
    E = knn_graph(X, k); n = len(X)
    d = X[E[:, 1]] - X[E[:, 0]]; L = np.linalg.norm(d, 1) + 1e-12
    f = ((drift[E[:, 0]] + drift[E[:, 1]]) / 2 * d).sum(1) / L      # flujo por arista
    B1 = np.zeros((len(E), n)); B1[np.arange(len(E)), E[:, 0]] = -1; B1[np.arange(len(E)), E[:, 1]] = 1
    phi, *_ = np.linalg.lstsq(B1, f, rcond=None)
    grad = B1 @ phi; res = f - grad
    return phi, np.linalg.norm(res) / (np.linalg.norm(grad) + 1e-12), E, f

# ----------------------------------------------------------------------------
def analyse_tissue(acc, expr, pheno, genes, r, reps, B, rng, out, eps_rel=0.5):
    order = STATE_ORDERS[acc]
    pheno = expand_characteristics(pheno); pheno["stage"] = pheno["condition"].map({s: k for k, s in enumerate(order)})
    pheno = pheno.dropna(subset=["stage"]); pheno["stage"] = pheno["stage"].astype(int)
    pheno = pheno[pheno[".sample_id"].isin(expr.columns)]
    Y = expr.loc[genes, pheno[".sample_id"]].to_numpy().T          # muestras x genes
    Y = (Y - Y.mean(0)) / (Y.std(0) + 1e-12)
    st = pheno["stage"].to_numpy()
    Yr, covs = (Y, []) if NO_COVAR else residualize(Y, pheno, stage=st)
    Z, varexp = pca(Yr, r)
    wst = stage_weights(st, balance=not NO_BALANCE)   # [FIX-design] masa igual por estadio
    h0 = silverman(Z)
    res = dict(accession=acc, n=len(Z), r=r, var_explained=float(varexp.sum()), covariates=";".join(covs) or "none",
               n_by_stage=";".join(f"{STAGE[s]}={int((st == s).sum())}" for s in range(3)))

    # --- cuencas: barrido de ancho de banda, persistencia H0 ---
    n_basins, barriers = [], []
    tau, eps_hat = persistence_threshold(Z[:, :2], h0, rng, w=wst)
    res.update(persistence_threshold_guaranteed=tau, stability_eps=eps_hat, stage_balanced=int(not NO_BALANCE))
    # Dos niveles de evidencia (no se mezclan):
    #  - nominal: persistencia > ln 2 (densidad en la silla < 1/2 del valle menor),
    #    con estabilidad en escala -> decision two_attractors (validado en sinteticos)
    #  - garantizado: persistencia > tau = max(ln2, 2 eps_hat) [teorema de estabilidad,
    #    Cohen-Steiner et al.; eps por bootstrap, Fasy et al. 2014] -> la cuenca no
    #    puede ser un artefacto de estimacion con confianza 90%. Condicion suficiente,
    #    conservadora: se reporta como fuerza de la evidencia.
    n_basins_g = []
    for h in h0 * np.array([0.7, 1.0, 1.4]):
        U, grid, pers, _ = sublevel_persistence_2d(Z, h, w=wst)
        ok_sup = lambda p: np.linalg.norm(Z[:, :2] - p[3], axis=1).min() < 1.5 * h
        sig = [p for p in pers if p[2] > np.log(2) and ok_sup(p)]
        n_basins_g.append(1 + len([p for p in pers if p[2] > tau and ok_sup(p)]))
        n_basins.append(1 + len(sig)); barriers.append(sig[0][2] if sig else 0.0)
    res.update(n_basins_h07=n_basins[0], n_basins_h10=n_basins[1], n_basins_h14=n_basins[2],
               barrier_persistence=barriers[1], basins_stable=int(len(set(n_basins)) == 1),
               n_basins_guaranteed=int(np.median(n_basins_g)))
    d2, d3 = gmm_dbic(Z); res.update(dBIC_2vs1=d2, dBIC_3vs1=d3)
    # Persistencia en el plano (eje entre medias de la GMM de 2 comp., PC ortogonal
    # principal): la bimodalidad puede vivir fuera de PC1-PC2. No usa etiquetas.
    n_basins_gp = np.nan
    try:
        from sklearn.mixture import GaussianMixture
        gm = GaussianMixture(2, covariance_type="full", n_init=5, random_state=0).fit(Z)
        ax = gm.means_[1] - gm.means_[0]; ax /= np.linalg.norm(ax) + 1e-12
        Zp = Z - np.outer(Z @ ax, ax); u2 = np.linalg.svd(Zp - Zp.mean(0), full_matrices=False)[2][0]
        Z2 = np.column_stack([Z @ ax, Z @ u2]); hp = silverman(Z2)
        nb = []
        for hh in hp * np.array([0.7, 1.0, 1.4]):
            _, _, pers, _ = sublevel_persistence_2d(Z2, hh, w=wst)
            sig = [p for p in pers if p[2] > np.log(2) and np.linalg.norm(Z2 - p[3], axis=1).min() < 1.5 * hh]
            nb.append(1 + len(sig))
        n_basins_gp = int(np.median(nb))   # mediana sobre escalas (robusta a una sola escala)
    except Exception:
        pass
    res.update(n_basins_gmmplane=n_basins_gp)
    # Regla de decision (prefijada): dos atractores <=> cuencas persistentes
    # estables a traves de escalas Y mezcla de 2 componentes preferida por BIC.
    # Ninguna prueba sola basta (un continuo puede dar 2 cuencas a un solo h).
    # Dos atractores <=> dBIC>0 Y cuencas persistentes estables en PC1-PC2 O en el plano GMM
    stable_pc = res["basins_stable"] and n_basins[1] >= 2
    stable_gp = (not np.isnan(n_basins_gp)) and n_basins_gp >= 2
    two = bool((not np.isnan(d2)) and d2 > 0 and (stable_pc or stable_gp))
    res.update(two_attractors=int(two))

    # --- puntos criticos por score y a que estadio pertenece cada valle ---
    modes, Umodes = critical_points(Z[:, :2], 1.3 * h0, w=wst)      # 2D, ancho mayor: menos fragmentacion
    mode_stage = []
    for m in modes:
        d = np.linalg.norm(Z[:, :2] - m, axis=1); near = st[np.argsort(d)[:max(5, len(Z) // 8)]]
        mode_stage.append(STAGE[int(np.bincount(near, minlength=3).argmax())])
    res.update(n_modes_score=len(modes), modes_stage=";".join(mode_stage))
    pd.DataFrame(np.column_stack([modes, Umodes]), columns=["pc1", "pc2", "U"]).assign(stage=mode_stage)\
      .to_csv(f"{out}/{acc}_critical_points.tsv", sep="\t", index=False)

    # --- asimetria de barreras entre valle sano y valle T2D (bootstrap) ---
    def barrier_asym(Zb, stb):
        h = silverman(Zb); wb = stage_weights(stb, balance=not NO_BALANCE)
        cH = Zb[stb == 0].mean(0); cD = Zb[stb == 2].mean(0)
        path = np.linspace(cH, cD, 60); p, _ = kde(Zb, path, h, wb); Up = -np.log(p)
        UH, UD, Umax = Up[0], Up[-1], Up.max()
        return Umax - UH, Umax - UD, np.argmax(Up) / 59
    dH, dD, pos = barrier_asym(Z, st)
    boot = np.array([barrier_asym(Z[i], st[i]) for i in (rng.choice(len(Z), len(Z)) for _ in range(reps))])
    nan = np.nan
    res.update(barrier_healthy_to_T2D=dH if two else nan, barrier_T2D_to_healthy=dD if two else nan,
               barrier_asymmetry=(dD - dH) if two else nan,
               barrier_asym_CI_lo=np.percentile(boot[:, 1] - boot[:, 0], 2.5) if two else nan,
               barrier_asym_CI_hi=np.percentile(boot[:, 1] - boot[:, 0], 97.5) if two else nan,
               saddle_position_on_path=pos if two else nan)
    # ¿los pacientes IGT caen cerca de la silla?
    if (st == 1).any():
        cH, cD = Z[st == 0].mean(0), Z[st == 2].mean(0); v = cD - cH
        t_igt = ((Z[st == 1] - cH) @ v) / (v @ v)
        res.update(igt_position_mean=t_igt.mean(), igt_position_sd=t_igt.std())

    # --- indice de transicion critica, n igualado, nulo por etiquetas ---
    n_min = min((st == s).sum() for s in range(3))
    def ic_by_stage(stb):
        return np.array([np.mean([critical_index(Yr[rng.choice(np.where(stb == s)[0], n_min, replace=False)])
                                  for _ in range(max(3, reps // 4))]) for s in range(3)])
    ic = ic_by_stage(st)
    null = np.array([ic_by_stage(rng.permutation(st)) for _ in range(max(20, B // 5))])
    stat = ic[1] - max(ic[0], ic[2]); nstat = null[:, 1] - np.maximum(null[:, 0], null[:, 2])
    res.update(Ic_healthy=ic[0], Ic_intermediate=ic[1], Ic_T2D=ic[2],
               p_Ic_intermediate_max=(np.sum(nstat >= stat) + 1) / (len(nstat) + 1))

    # --- Fisher a lo largo de covariables continuas (todas las disponibles) ---
    ctrl_cols = {}
    if CONTROL_OVERRIDE is not None:
        c = find_col(pheno, [CONTROL_OVERRIDE.lower()]); 
        if c is not None: ctrl_cols["control"] = c
    else:
        for grp, names in CONTROL_GROUPS.items():
            c = find_col(pheno, names)
            if c is not None: ctrl_cols[grp] = c
    fisher_summary = []
    for grp, c in ctrl_cols.items():
        th = pd.to_numeric(pheno[c].astype(str).str.replace(",", "."), errors="coerce").to_numpy(); ok = ~np.isnan(th)
        if ok.sum() < 30: continue
        cen, g = fisher_along(Z[ok], th[ok], h0)
        if len(g) == 0: continue
        pd.DataFrame(dict(theta=cen, fisher=g)).to_csv(f"{out}/{acc}_fisher_{grp}.tsv", sep="\t", index=False)
        # ¿en que estadio cae el pico? (estadio modal de los pacientes con theta cercano al pico)
        pk = cen[np.argmax(g)]; near = st[ok][np.argsort(np.abs(th[ok] - pk))[:max(8, ok.sum() // 6)]]
        res[f"fisher_{grp}_col"] = c; res[f"fisher_{grp}_peak_theta"] = pk
        res[f"fisher_{grp}_peak_stage"] = STAGE[int(np.bincount(near, minlength=3).argmax())]
        fisher_summary.append(f"{grp}:{pk:.2f}->{res[f'fisher_{grp}_peak_stage']}")
    res["fisher_controls"] = ";".join(ctrl_cols.values()) or "none"
    res["fisher_peaks"] = ";".join(fisher_summary) or "none"

    # --- ETAPA 2: puente sano -> T2D y Hodge de la deriva ---
    A, Bz = Z[st == 0], Z[st == 2]
    eps = eps_rel * np.median(cdist(A, Bz, "sqeuclidean"))
    pi = sinkhorn(A, Bz, eps)
    drift_A = (pi @ Bz) / pi.sum(1)[:, None] - A
    drift = np.zeros_like(Z); drift[st == 0] = drift_A
    drift[st == 2] = -((pi.T @ A) / pi.sum(0)[:, None] - Bz)      # deriva inversa en T2D
    phi, ratio_J, E, f = hodge_decompose_drift(Z, drift)
    # ruta media del puente y distancia de los IGT a ella
    mid = (A[:, None, :] + Bz[None, :, :]) / 2
    route = (pi[:, :, None] * mid).sum((0, 1)) / pi.sum()
    res.update(bridge_eps=eps, flux_ratio_J_over_gradU=ratio_J,
               igt_dist_to_bridge=np.linalg.norm(Z[st == 1] - route, axis=1).mean() / np.linalg.norm(Z.std(0)) if (st == 1).any() else np.nan)
    # nulo del ratio: permutar etiquetas sano/T2D
    rnull = []
    for _ in range(max(20, B // 10)):
        p = rng.permutation(np.where(st != 1)[0]); nA = (st == 0).sum()
        A2, B2 = Z[p[:nA]], Z[p[nA:]]; pi2 = sinkhorn(A2, B2, eps)
        dr = np.zeros_like(Z); dr[p[:nA]] = (pi2 @ B2) / pi2.sum(1)[:, None] - A2
        rnull.append(hodge_decompose_drift(Z, dr)[1])
    # p pequeno = el flujo no-gradiente observado es MAYOR que bajo etiquetas
    # permutadas (evidencia de dinamica fuera de equilibrio); p ~ 1 = compatible
    # con dinamica de gradiente (la deriva inferida es mas 'recta' que el azar).
    res.update(p_flux_ratio_gt_null=(np.sum(np.array(rnull) >= ratio_J) + 1) / (len(rnull) + 1),
               flux_ratio_null_mean=float(np.mean(rnull)))

    pd.DataFrame(Z, columns=[f"pc{i+1}" for i in range(r)]).assign(sample=pheno[".sample_id"].to_numpy(), stage=[STAGE[s] for s in st],
        U=-np.log(kde(Z, Z, h0, wst)[0]), stage_weight=wst, phi_dynamic=phi).to_csv(f"{out}/{acc}_embedding_potential.tsv", sep="\t", index=False)
    return res

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--export_dir", default="export_sheaf"); ap.add_argument("--out", default="results/landscape")
    ap.add_argument("--tissues", nargs="+", default=["GSE76895", "GSE18732", "GSE27951", "GSE15653"])
    ap.add_argument("--n_genes", type=int, default=800); ap.add_argument("--r", type=int, default=3)
    ap.add_argument("--reps", type=int, default=20); ap.add_argument("--B", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1234)
    ap.add_argument("--no_covar", action="store_true", help="no regresar covariables (analisis de sensibilidad)")
    ap.add_argument("--no_balance", action="store_true", help="no reponderar estadios a masa igual (sensibilidad)")
    ap.add_argument("--control", default=None, help="columna de pheno a usar como control continuo para Fisher (por defecto: hba1c y glucosa si existen)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    global NO_COVAR, CONTROL_OVERRIDE, NO_BALANCE; NO_COVAR = a.no_covar; CONTROL_OVERRIDE = a.control; NO_BALANCE = a.no_balance; rng = np.random.default_rng(a.seed)
    hv = pd.read_csv(os.path.join(a.export_dir, "high_variance_genes_ordered.tsv"), sep="\t")["gene"].tolist()[:a.n_genes]
    rows = []
    for acc in a.tissues:
        expr = pd.read_csv(os.path.join(a.export_dir, f"{acc}_expr.tsv"), sep="\t", index_col=0)
        pheno = pd.read_csv(os.path.join(a.export_dir, f"{acc}_pheno.tsv"), sep="\t").dropna(subset=["condition"])
        genes = [g for g in hv if g in expr.index]
        print(f"[{acc}] n={len(pheno)} genes={len(genes)}")
        r = analyse_tissue(acc, expr, pheno, genes, a.r, a.reps, a.B, rng, a.out); rows.append(r)
        print(f"   cuencas PC={r['n_basins_h10']} garantizadas={r['n_basins_guaranteed']}(tau={r['persistence_threshold_guaranteed']:.2f}) estables={r['basins_stable']} planoGMM={r['n_basins_gmmplane']} dBIC2={r['dBIC_2vs1']:.1f} "
              f"modos={r['modes_stage']} | barrera asim={r['barrier_asymmetry']:.2f} "
              f"[{r['barrier_asym_CI_lo']:.2f},{r['barrier_asym_CI_hi']:.2f}] 2attr={r['two_attractors']} | Ic={r['Ic_healthy']:.2f}/{r['Ic_intermediate']:.2f}/{r['Ic_T2D']:.2f} "
              f"p={r['p_Ic_intermediate_max']:.3f} | Fisher {r['fisher_peaks']} | J/gradU={r['flux_ratio_J_over_gradU']:.2f} p_gt_null={r['p_flux_ratio_gt_null']:.3f}")
    pd.DataFrame(rows).to_csv(os.path.join(a.out, "landscape_summary.tsv"), sep="\t", index=False)
    json.dump(vars(a), open(os.path.join(a.out, "config.json"), "w"), indent=2)

if __name__ == "__main__":
    main()
