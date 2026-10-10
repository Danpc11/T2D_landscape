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

def perm_test_coherence(DA, DB, rng, B=2000):
    """p para coherencia(A) - coherencia(B) > 0, permutando etiquetas de individuo."""
    allD = np.vstack([DA, DB]); nA = len(DA); obs = coherence_loo(DA) - coherence_loo(DB); null = []
    for _ in range(B):
        idx = rng.permutation(len(allD)); null.append(coherence_loo(allD[idx[:nA]]) - coherence_loo(allD[idx[nA:]]))
    return obs, (np.sum(np.array(null) >= obs) + 1) / (B + 1)

def perm_test_direction(DA, DB, rng, B=2000):
    """p para cos(dir A, dir B) < nulo (direcciones mas distintas de lo esperado)."""
    allD = np.vstack([DA, DB]); nA = len(DA); obs = cosine(DA.mean(0), DB.mean(0)); null = []
    for _ in range(B):
        idx = rng.permutation(len(allD)); null.append(cosine(allD[idx[:nA]].mean(0), allD[idx[nA:]].mean(0)))
    return obs, float(np.mean(null)), (np.sum(np.array(null) <= obs) + 1) / (B + 1)

def paired_swap_test(DA, DB, stat, rng, B=3000):   # pareado: 3000, como dice la leyenda
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

# ---------------- visual helpers (no dependen de matplotlib hasta llamarse) ----------------
def rose(ax, cosines, color, n_bins=18, label=None):
    """Rosa de los vientos de direcciones: angulo = arccos(cos) con signo aleatorio-estable (simetrica), norte = direccion sana."""
    import numpy as np
    ang = np.arccos(np.clip(np.asarray(cosines), -1, 1)); ang = np.concatenate([ang, -ang])   # simetrica: solo importa el angulo a la direccion sana
    edges = np.linspace(-np.pi, np.pi, n_bins + 1); h, _ = np.histogram(ang, bins=edges); h = h / 2
    ax.bar((edges[:-1] + edges[1:]) / 2, h, width=edges[1] - edges[0], bottom=0, color=color, alpha=0.75, edgecolor="white", lw=0.3, label=label)
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1); ax.set_xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2]); ax.set_xticklabels(["healthy\ndirection", "90°", "opposite", "270°"], fontsize=8); ax.set_yticks([])
    ax.spines["polar"].set_linewidth(0.4)

def chord(ax, names, weights, colors, null_weights=None, gap=0.08):
    """Diagrama de cuerdas simple: nodos en circulo, cuerdas con grosor ~ peso (0-1); null_weights dibuja la cuerda del nulo en gris."""
    import numpy as np; from matplotlib.patches import PathPatch; from matplotlib.path import Path
    n = len(names); th = np.linspace(0, 2 * np.pi, n, endpoint=False) + np.pi / 2
    for i, nm in enumerate(names):
        ax.add_patch(__import__("matplotlib").patches.Wedge((0, 0), 1.0, np.degrees(th[i]) - 20, np.degrees(th[i]) + 20, width=0.08, fc=colors[i], ec="none"))
        ax.text(1.22 * np.cos(th[i]), 1.22 * np.sin(th[i]), nm, ha="center", va="center", fontsize=6)
    def band(i, j, w, col, alpha, z):
        a0, a1 = th[i] - 0.25 * w, th[i] + 0.25 * w; b0, b1 = th[j] - 0.25 * w, th[j] + 0.25 * w; r = 0.9
        p = [(r * np.cos(a0), r * np.sin(a0))]; codes = [Path.MOVETO]
        p += [(r * np.cos(a1), r * np.sin(a1)), (0, 0), (r * np.cos(b0), r * np.sin(b0))]; codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
        p += [(r * np.cos(b1), r * np.sin(b1)), (0, 0), (r * np.cos(a0), r * np.sin(a0))]; codes += [Path.LINETO, Path.CURVE3, Path.CURVE3]
        ax.add_patch(PathPatch(Path(p, codes), fc=col, ec="none", alpha=alpha, zorder=z))
    for (i, j), w in weights.items():
        if null_weights: band(i, j, null_weights[(i, j)], "#bbbbbb", 0.9, 1)
        band(i, j, w, colors[i], 0.55, 2)
    ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45); ax.set_aspect("equal"); ax.axis("off")

def alluvial(ax, flows, left_label, right_labels, colors, title=None):
    """flows: dict destino -> fraccion (suman 1). Una fuente a la izquierda, destinos apilados a la derecha, bandas Bezier."""
    import numpy as np; from matplotlib.patches import PathPatch; from matplotlib.path import Path
    y = 0.0; ax.add_patch(__import__("matplotlib").patches.Rectangle((0, 0), 0.08, 1, fc="#444444", ec="none")); ax.text(-0.04, 0.5, left_label, ha="right", va="center", fontsize=8, rotation=90)
    gap = 0.03; tot = sum(flows.values()); yr = 0.0
    for k, lab in enumerate(right_labels):
        f = flows[lab] / tot; h_r = f * (1 - gap * (len(right_labels) - 1))
        ax.add_patch(__import__("matplotlib").patches.Rectangle((0.92, yr), 0.08, h_r, fc=colors[lab], ec="none")); ax.text(1.03, yr + h_r / 2, f"{lab}\n{f:.0%}", va="center", fontsize=8)
        p = [(0.08, y), (0.5, y), (0.5, yr), (0.92, yr), (0.92, yr + h_r), (0.5, yr + h_r), (0.5, y + f), (0.08, y + f), (0.08, y)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(p, codes), fc=colors[lab], ec="none", alpha=0.45)); y += f; yr += h_r + gap
    ax.set_xlim(-0.1, 1.4); ax.set_ylim(0, 1); ax.axis("off")
    if title: ax.set_title(title)

def random_rotation(d, rng):
    """Matriz ortogonal uniforme (Haar) de dimension d, para nulos isotropos."""
    A = rng.normal(size=(d, d)); Q, Rm = np.linalg.qr(A); return Q * np.sign(np.diag(Rm))

def isotropic_null(D, rng):
    """Desplazamientos con la misma longitud pero direcciones isotropas: destruye cualquier
    direccion compartida sin conservar los ejes originales (a diferencia de cambiar el signo)."""
    n, d = D.shape; L = np.linalg.norm(D, axis=1, keepdims=True)
    V = rng.normal(size=(n, d)); V /= (np.linalg.norm(V, axis=1, keepdims=True) + 1e-12)
    return V * L
