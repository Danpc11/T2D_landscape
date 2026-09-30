import numpy as np
from sheaf_prototype_synthetic import *
import sheaf_prototype_synthetic as sp

def spectral_embed(W, r=6):
    d = W.sum(1); Ln = np.eye(P) - (W / np.sqrt(np.outer(d, d)))
    vals, vecs = np.linalg.eigh(Ln)
    return vecs[:, 1:r+1]                      # P x r

def procrustes(A, B):
    U, _, Vt = np.linalg.svd(A.T @ B); return U @ Vt

def stalks(Ws, kind):
    if kind == "strength":
        return strength_profile(Ws)[:, :, None]               # P x T x 1
    if kind == "neighborhood":                                # fila de W (perfil de vecinos)
        X = np.stack([W / W.sum(1, keepdims=True) for W in Ws], 1)   # P x T x P
        return X
    if kind == "embedding":                                   # posicion espectral alineada (Procrustes)
        E = [spectral_embed(W) for W in Ws]
        E = [E[0]] + [e @ procrustes(e, E[0]) for e in E[1:]] # mapas de restriccion O(r)
        return np.stack(E, 1)                                 # P x T x r

def energy_gene_sheaf(Ws, S):
    """haz sobre el grafo de genes (consenso); energia sum_ij w_ij sum_t ||x_it - x_jt||^2 / ||x||^2"""
    i, j, w = base_graph(Ws)
    d = ((S[i] - S[j]) ** 2).sum((1, 2)) * w
    return d.sum() / (S ** 2).sum() * (P / w.sum())

def energy_tissue_sheaf(S):
    """haz sobre el grafo completo de tejidos; energia sum_tt' ||x_.t - x_.t'||^2 / ||x||^2 (sin topologia genica)"""
    e = 0
    for t in range(T):
        for u in range(t+1, T):
            e += ((S[:, t] - S[:, u]) ** 2).sum()
    return e / (S ** 2).sum()

def run(kind, reps=8):
    out = {s: {"gene": [], "tissue": []} for s in ["ND", "IGT", "T2D"]}
    for _ in range(reps):
        m0 = np.repeat(np.arange(K), P // K)
        for s in out:
            Ws = [adjacency(X) for X in make_state(s, m0)]
            S = stalks(Ws, kind)
            out[s]["gene"].append(energy_gene_sheaf(Ws, S))
            out[s]["tissue"].append(energy_tissue_sheaf(S))
    print(f"\nstalk = {kind}")
    for s in out:
        g, t_ = np.array(out[s]["gene"]), np.array(out[s]["tissue"])
        print(f"  {s:4s}  haz-genes: {g.mean():.4f} ± {g.std():.4f}   haz-tejidos: {t_.mean():.4f} ± {t_.std():.4f}")

for kind in ["strength", "neighborhood", "embedding"]:
    run(kind)
