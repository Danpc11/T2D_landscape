"""Geometry of within-person responses: embedding, displacement, leave-one-out coherence, nulls."""
import numpy as np, pandas as pd

def embed(expr, n_genes=800, r=5, expr_floor_pct=30):
    """expr: genes x muestras (log). Devuelve DataFrame muestras x r (PCA de los n_genes mas variables)."""
    X = expr.to_numpy().T
    if X.max() > 50: X = np.log2(np.clip(X, 1, None))
    X = X[:, (X > np.percentile(X, expr_floor_pct)).mean(0) > 0.5]
    v = X.var(0); X = X[:, np.argsort(-v)[:n_genes]]; X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    U, S, Vt = np.linalg.svd(X - X.mean(0), full_matrices=False)
    return pd.DataFrame(U[:, :r] * S[:r], index=expr.columns)

def displacements(Z, pairs):
    """pairs: lista de (id_basal, id_estimulo). Devuelve matriz n x r."""
    return np.array([Z.loc[b_].values - Z.loc[a_].values for a_, b_ in pairs])

def coherence_loo(D):
    """Media del coseno entre cada respuesta y la media del grupo SIN esa respuesta."""
    out = []
    for i in range(len(D)):
        m = np.delete(D, i, 0).mean(0); out.append((D[i] @ m) / (np.linalg.norm(D[i]) * np.linalg.norm(m) + 1e-9))
    return float(np.mean(out))

def coherence_naive(D):
    m = D.mean(0); return float(np.mean([(x @ m) / (np.linalg.norm(x) * np.linalg.norm(m) + 1e-9) for x in D]))

def cosine(a, b): return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))

def perm_test_coherence(DA, DB, rng, B=3000):
    """p para coherencia(A) - coherencia(B) > 0, permutando etiquetas de individuo."""
    allD = np.vstack([DA, DB]); nA = len(DA); obs = coherence_loo(DA) - coherence_loo(DB); null = []
    for _ in range(B):
        idx = rng.permutation(len(allD)); null.append(coherence_loo(allD[idx[:nA]]) - coherence_loo(allD[idx[nA:]]))
    return obs, (np.sum(np.array(null) >= obs) + 1) / (B + 1)

def perm_test_direction(DA, DB, rng, B=3000):
    """p para cos(dir A, dir B) < nulo (direcciones mas distintas de lo esperado)."""
    allD = np.vstack([DA, DB]); nA = len(DA); obs = cosine(DA.mean(0), DB.mean(0)); null = []
    for _ in range(B):
        idx = rng.permutation(len(allD)); null.append(cosine(allD[idx[:nA]].mean(0), allD[idx[nA:]].mean(0)))
    return obs, float(np.mean(null)), (np.sum(np.array(null) <= obs) + 1) / (B + 1)

def paired_swap_test(DA, DB, stat, rng, B=3000):
    """A y B: respuestas de LAS MISMAS personas en dos condiciones (mismo orden). Intercambia A/B dentro de persona."""
    n = min(len(DA), len(DB)); obs = stat(DB[:n]) - stat(DA[:n]); null = []
    for _ in range(B):
        sw = rng.random(n) < 0.5; a = np.where(sw[:, None], DB[:n], DA[:n]); b = np.where(sw[:, None], DA[:n], DB[:n]); null.append(stat(b) - stat(a))
    return obs, (np.sum(np.abs(null) >= abs(obs)) + 1) / (B + 1)

def gene_response_table(g, pairs):
    """g: genes x muestras. Respuesta por gen: media, t pareado, consistencia de signo."""
    d = np.array([g[b_].values - g[a_].values for a_, b_ in pairs])
    return pd.DataFrame({"mean": d.mean(0), "t": d.mean(0) / (d.std(0, ddof=1) / np.sqrt(len(d)) + 1e-9),
                         "consist": (np.sign(d) == np.sign(d.mean(0))).mean(0), "n": len(d)}, index=g.index), d

def interaction_perm(dA, dB, rng, B=5000):
    """dA, dB: respuestas por individuo (vector) en dos grupos. p de |mean A - mean B| por permutacion."""
    obs = abs(dA.mean() - dB.mean()); pool = np.concatenate([dA, dB]); null = []
    for _ in range(B):
        r = rng.permutation(pool); null.append(abs(r[:len(dA)].mean() - r[len(dA):].mean()))
    return obs, (np.sum(np.array(null) >= obs) + 1) / (B + 1)

def bh(p):
    p = np.asarray(p); o = np.argsort(p); q = np.minimum.accumulate((p[o] * len(p) / (np.arange(len(p)) + 1))[::-1])[::-1]
    out = np.empty(len(p)); out[o] = np.minimum(q, 1); return out
