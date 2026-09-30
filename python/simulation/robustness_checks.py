import numpy as np
import sheaf_prototype_synthetic as sp
from sheaf_prototype_synthetic import *
from stalk_comparison import stalks, energy_gene_sheaf, energy_tissue_sheaf

def E_emb(tissues):
    Ws = [adjacency(X) for X in tissues]
    return energy_gene_sheaf(Ws, stalks(Ws, "embedding"))

print("Robustez a n (stalk embedding, haz-genes; media ± sd de 6 replicas)")
for n in [8, 15, 30, 45]:
    sp.N = n
    acc = {s: [] for s in ["ND","IGT","T2D"]}
    for _ in range(6):
        m0 = np.repeat(np.arange(K), P // K)
        for s in acc: acc[s].append(E_emb(make_state(s, m0)))
    print(f"  n={n:2d}: " + "  ".join(f"{s}={np.mean(v):.3f}±{np.std(v):.3f}" for s,v in acc.items()))
sp.N = 30

print("\nNulo (a): permutar etiquetas de estado dentro de cada tejido (n=30, B=150)")
m0 = np.repeat(np.arange(K), P // K)
states = {s: make_state(s, m0) for s in ["ND","IGT","T2D"]}
names = list(states)
obs = E_emb(states["IGT"]) - max(E_emb(states["ND"]), E_emb(states["T2D"]))
pooled = [np.concatenate([states[s][t] for s in names], axis=1) for t in range(T)]
null = []
for _ in range(150):
    perm = rng.permutation(pooled[0].shape[1])
    Es = {s: E_emb([X[:, perm[k*N:(k+1)*N]] for X in pooled]) for k,s in enumerate(names)}
    null.append(Es["IGT"] - max(Es["ND"], Es["T2D"]))
null = np.array(null)
print(f"  observado = {obs:.3f}   nulo media={null.mean():.3f} sd={null.std():.3f}   p = {(null>=obs).mean():.3f}")

print("\nControl negativo: IGT sin reorganizacion (misma red que ND, muestras distintas)")
acc = []
for _ in range(6):
    m0 = np.repeat(np.arange(K), P // K)
    acc.append((E_emb(make_state("ND", m0)), E_emb(make_state("ND", m0))))
acc = np.array(acc)
print(f"  E(ND_a)={acc[:,0].mean():.3f}  E(ND_b)={acc[:,1].mean():.3f}  (deben coincidir)")

print("\nLocalizacion: ¿los genes reorganizados cargan la energia local en IGT?")
# reconstruir IGT con registro de genes reorganizados
def make_igt_tracked(m0):
    tissues, reorg = [], np.zeros(P, bool)
    for t in range(T):
        m = m0.copy(); idx = rng.choice(P, int(F_REORG*P), replace=False)
        m[idx] = rng.integers(0, K, size=idx.size); reorg[idx] = True
        tissues.append(sp.simulate_tissue(m, N))
    return tissues, reorg
from sklearn.metrics import roc_auc_score
aucs = []
for _ in range(6):
    m0 = np.repeat(np.arange(K), P // K)
    tissues, reorg = make_igt_tracked(m0)
    Ws = [adjacency(X) for X in tissues]; S = stalks(Ws, "embedding")
    i, j, w = base_graph(Ws)
    d = ((S[i]-S[j])**2).sum((1,2))*w
    e = np.zeros(P); np.add.at(e, i, d); np.add.at(e, j, d)
    # energia por gen sin grafo: varianza entre tejidos del embedding
    e_t = S.var(1).sum(1)
    aucs.append((roc_auc_score(reorg, e), roc_auc_score(reorg, e_t)))
aucs = np.array(aucs)
print(f"  AUC gen reorganizado ~ energia local:  haz-genes={aucs[:,0].mean():.3f}   sin grafo={aucs[:,1].mean():.3f}")
