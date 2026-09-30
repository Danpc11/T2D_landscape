import numpy as np
import sheaf_prototype_synthetic as sp
from sheaf_prototype_synthetic import *
from stalk_comparison import stalks, energy_gene_sheaf, energy_tissue_sheaf

def make_state_n(kind, m0, n_per_tissue, T_):
    sp.T = T_
    tissues = []
    if kind == "ND":
        for t in range(T_): tissues.append(sp.simulate_tissue(m0, n_per_tissue[t]))
    elif kind == "IGT":
        for t in range(T_):
            m = m0.copy(); idx = rng.choice(P, int(sp.F_REORG*P), replace=False)
            m[idx] = rng.integers(0, K, size=idx.size)
            tissues.append(sp.simulate_tissue(m, n_per_tissue[t]))
    else:
        m = m0.copy(); idx = rng.choice(P, int(sp.F_REORG*P), replace=False)
        m[idx] = rng.integers(0, K, size=idx.size)
        for t in range(T_): tissues.append(sp.simulate_tissue(m, n_per_tissue[t], 0.7))
    return tissues

def E_tissue_emb(tissues):
    Ws = [adjacency(X) for X in tissues]
    return energy_tissue_sheaf(stalks(Ws, "embedding"))

if __name__ == "__main__":
  designs = {
   "4 tejidos, n igualado por tejido (15,26,4,11)": [15,26,4,11],
   "3 tejidos sin higado (15,26,11)":               [15,26,11],
   "2 tejidos pancreas+musculo (15,26)":            [15,26],
   "3 tejidos, n igualado GLOBAL a 11":             [11,11,11],
   "referencia: 3 tejidos n=30":                    [30,30,30],
  }
  print("Haz sobre tejidos, stalk embedding. 12 replicas. Se reporta E y fraccion de replicas con IGT maximo.\n")
  for name, ns in designs.items():
      T_ = len(ns)
      import stalk_comparison; stalk_comparison.T = T_
      res = {s: [] for s in ["ND","IGT","T2D"]}
      for _ in range(12):
          m0 = np.repeat(np.arange(K), P // K)
          for s in res: res[s].append(E_tissue_emb(make_state_n(s, m0, ns, T_)))
      R = np.array([res[s] for s in ["ND","IGT","T2D"]])   # 3 x reps
      frac_igt_max = (R.argmax(0) == 1).mean()
      frac_t2d_gt_nd = (R[2] > R[0]).mean()
      print(f"{name:48s} ND={R[0].mean():.2f}  IGT={R[1].mean():.2f}  T2D={R[2].mean():.2f}"
            f"  | P(IGT max)={frac_igt_max:.2f}  P(T2D>ND)={frac_t2d_gt_nd:.2f}")
