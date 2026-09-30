"""
Prototipo: haz celular sobre redes de coexpresion multi-tejido.
Pregunta: ¿la energia del haz distingue un estado con reorganizacion
tejido-ESPECIFICA (incoherente, 'IGT') de estados coherentes ('ND', 'T2D')?
"""
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh

rng = np.random.default_rng(7)

# ---------------- parametros ----------------
P, T, N = 300, 4, 30          # genes, tejidos, muestras por estado
K = 6                          # modulos
BETA = 4                       # potencia fija (misma para todos los estados)
EDGE_FRAC = 0.02               # aristas del grafo base (consenso)
F_REORG = 0.35                 # fraccion de genes reorganizados
B_NULL = 300                   # permutaciones

# ---------------- generador ----------------
def simulate_tissue(assign, n, a=1.0, noise=1.0):
    z = rng.normal(size=(n, K))
    X = a * z[:, assign] + noise * rng.normal(size=(n, P))
    return X.T                                   # genes x muestras

def make_state(kind, m0, a_common=1.0):
    """Devuelve lista de T matrices de expresion (genes x muestras)."""
    tissues = []
    if kind == "ND":
        for t in range(T):
            tissues.append(simulate_tissue(m0, N, a_common))
    elif kind == "IGT":                          # reorganizacion independiente por tejido
        for t in range(T):
            m = m0.copy()
            idx = rng.choice(P, int(F_REORG * P), replace=False)
            m[idx] = rng.integers(0, K, size=idx.size)
            tissues.append(simulate_tissue(m, N, a_common))
    elif kind == "T2D":                          # reorganizacion COMUN a los tejidos
        m = m0.copy()
        idx = rng.choice(P, int(F_REORG * P), replace=False)
        m[idx] = rng.integers(0, K, size=idx.size)
        for t in range(T):
            tissues.append(simulate_tissue(m, N, a_common * 0.7))  # menos modular
    return tissues

# ---------------- redes ----------------
def adjacency(X, beta=BETA):
    C = np.corrcoef(X)
    C[np.isnan(C)] = 0
    W = np.abs(C) ** beta
    np.fill_diagonal(W, 0)
    return W

def strength_profile(Ws):
    """Stalk por gen: fuerza estandarizada en cada tejido -> P x T."""
    S = np.column_stack([W.sum(1) for W in Ws])
    S = (S - S.mean(0)) / S.std(0)
    return S

def base_graph(Ws, frac=EDGE_FRAC):
    Wc = np.mean(Ws, axis=0)
    iu = np.triu_indices(P, 1)
    vals = Wc[iu]
    k = int(frac * vals.size)
    thr = np.partition(vals, -k)[-k]
    mask = vals >= thr
    return iu[0][mask], iu[1][mask], vals[mask]

# ---------------- haz ----------------
def sheaf_laplacian(i, j, w, maps=None):
    """
    Haz constante (maps=None): F_{i->ij} = sqrt(w) I_T.
    Laplaciano de bloques PT x PT, sparse.
    """
    E = i.size
    rows, cols, data = [], [], []
    for e in range(E):
        for t in range(T):
            a, b = i[e] * T + t, j[e] * T + t
            rows += [a, b, a, b]; cols += [a, b, b, a]
            data += [w[e], w[e], -w[e], -w[e]]
    L = sparse.coo_matrix((data, (rows, cols)), shape=(P * T, P * T)).tocsr()
    return L

def sheaf_energy(L, S):
    x = S.reshape(-1)                            # orden gen-mayor: (i,t) -> i*T+t
    return float(x @ (L @ x)) / float(x @ x)     # cociente de Rayleigh (escala-libre)

def local_energy(i, j, w, S):
    """Energia por gen: suma de ||x_i - x_j||^2 w_ij sobre sus aristas."""
    d = ((S[i] - S[j]) ** 2).sum(1) * w
    e = np.zeros(P)
    np.add.at(e, i, d); np.add.at(e, j, d)
    return e

def sheaf_gap(L):
    vals = eigsh(L.asfptype(), k=T + 2, sigma=-1e-6, which="LM",
                 return_eigenvectors=False)
    vals = np.sort(vals)
    # los primeros T (por componente conexa) son ~0 -> reportar el siguiente
    return vals[T]

# ---------------- baseline sin grafo ----------------
def profile_coherence(S):
    """Correlacion media entre tejidos de los perfiles de fuerza (sin topologia)."""
    C = np.corrcoef(S.T)
    return C[np.triu_indices(T, 1)].mean()

# ---------------- pipeline por estado ----------------
def analyse(tissues):
    Ws = [adjacency(X) for X in tissues]
    S = strength_profile(Ws)
    i, j, w = base_graph(Ws)
    L = sheaf_laplacian(i, j, w)
    return dict(E=sheaf_energy(L, S), gap=sheaf_gap(L),
                coh=profile_coherence(S), S=S, edges=(i, j, w), Ws=Ws)

# ---------------- nulo (b): romper correspondencia gen<->gen entre tejidos ----------------
def null_correspondence(res, B=B_NULL):
    i, j, w = res["edges"]; S0 = res["S"]
    L = sheaf_laplacian(i, j, w)
    out = []
    for _ in range(B):
        S = S0.copy()
        for t in range(1, T):
            S[:, t] = S[rng.permutation(P), t]
        out.append(sheaf_energy(L, S))
    out = np.array(out)
    return (res["E"] - out.mean()) / out.std(), (out <= res["E"]).mean()

# ---------------- nulo (a): permutar etiquetas de estado dentro de tejido ----------------
def null_labels(states, B=200):
    """Distribucion nula de E(IGT) - max(E(ND), E(T2D)) bajo H0 de intercambiabilidad."""
    names = list(states)
    pooled = [np.concatenate([states[s][t] for s in names], axis=1) for t in range(T)]
    n_tot = pooled[0].shape[1]
    obs = analyse(states["IGT"])["E"] - max(analyse(states["ND"])["E"],
                                            analyse(states["T2D"])["E"])
    null = []
    for _ in range(B):
        perm = rng.permutation(n_tot)
        splits = {s: perm[k * N:(k + 1) * N] for k, s in enumerate(names)}
        Es = {s: analyse([X[:, splits[s]] for X in pooled])["E"] for s in names}
        null.append(Es["IGT"] - max(Es["ND"], Es["T2D"]))
    null = np.array(null)
    return obs, (null >= obs).mean()

# ============================================================
if __name__ == "__main__":
    m0 = np.repeat(np.arange(K), P // K)
    states = {s: make_state(s, m0) for s in ["ND", "IGT", "T2D"]}

    print(f"P={P} genes, T={T} tejidos, n={N}/estado, beta={BETA} fijo, "
          f"grafo base top-{EDGE_FRAC*100:.0f}% aristas\n")
    print(f"{'estado':6s} {'E_haz':>8s} {'gap':>8s} {'coh_perfil':>11s} "
          f"{'z_vs_nulo_b':>12s} {'p_b':>6s}")
    results = {}
    for s in ["ND", "IGT", "T2D"]:
        r = analyse(states[s]); results[s] = r
        z, p = null_correspondence(r)
        print(f"{s:6s} {r['E']:8.4f} {r['gap']:8.4f} {r['coh']:11.3f} {z:12.2f} {p:6.3f}")

    print("\nNulo (a): permutacion de etiquetas de estado dentro de cada tejido")
    obs, p = null_labels(states)
    print(f"  E(IGT) - max(E(ND),E(T2D)) observado = {obs:.4f}   p_perm = {p:.3f}")

    # Robustez a n: ¿se sostiene con n pequeno?
    print("\nRobustez al tamano muestral (E_haz, media de 10 replicas):")
    for n in [10, 20, 30, 45]:
        acc = {s: [] for s in ["ND", "IGT", "T2D"]}
        for _ in range(10):
            N_bak = N; globals()["N"] = n
            st = {s: make_state(s, m0) for s in acc}
            for s in acc: acc[s].append(analyse(st[s])["E"])
            globals()["N"] = N_bak
        print(f"  n={n:2d}: " + "  ".join(f"{s}={np.mean(v):.4f}" for s, v in acc.items()))

    # Energia local: ¿los genes reorganizados por tejido cargan la incoherencia?
    r = results["IGT"]
    e_loc = local_energy(*r["edges"], r["S"])
    print("\nEnergia local IGT: top-10 genes mas incoherentes:", np.argsort(-e_loc)[:10])
