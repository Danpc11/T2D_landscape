"""Curva de potencia (suplemento): P(IGT maximo) vs fraccion de genes reorganizados,
con los n reales por tejido (15, 26, 11) y sin higado. Uso: python power_curve.py"""
import numpy as np, sys
sys.path.insert(0, __import__("os").path.dirname(__file__))
import sheaf_prototype_synthetic as sp
import stalk_comparison as sc
from power_real_n import make_state_n, E_tissue_emb
ns = [15, 26, 11]; sc.T = len(ns)
print("f_reorg  P(IGT max)  P(T2D>ND)   E_ND   E_IGT  E_T2D   (20 replicas)")
for f in [0.10, 0.15, 0.20, 0.25, 0.35, 0.50]:
    sp.F_REORG = f
    R = []
    for _ in range(20):
        m0 = np.repeat(np.arange(sp.K), sp.P // sp.K)
        R.append([E_tissue_emb(make_state_n(s, m0, ns, len(ns))) for s in ["ND", "IGT", "T2D"]])
    R = np.array(R)
    print(f"{f:6.2f}   {(R.argmax(1)==1).mean():8.2f}   {(R[:,2]>R[:,0]).mean():8.2f}   "
          f"{R[:,0].mean():.3f}  {R[:,1].mean():.3f}  {R[:,2].mean():.3f}")
