#!/usr/bin/env python3
"""Prueba rapida de que las funciones centrales hacen lo que el articulo afirma.
No necesita datos: construye casos con respuesta conocida. Uso: python tests/smoke_test.py
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../python/lib"))
import response as R
rng = np.random.default_rng(0); fails = []
def check(name, cond, detail=""):
    print(("  OK   " if cond else "  FALLA ") + name + (f"  [{detail}]" if detail else ""))
    if not cond: fails.append(name)

print("coherencia leave-one-out")
# sin direccion compartida debe dar ~0 a cualquier n; la version ingenua no
for n in (4, 10, 20):
    D = R.isotropic_null(rng.normal(size=(n, 5)) + np.r_[3, 0, 0, 0, 0], rng)
    loo = np.mean([R.coherence_loo(R.isotropic_null(D, rng)) for _ in range(200)])
    nai = np.mean([R.coherence_naive(R.isotropic_null(D, rng)) for _ in range(200)])
    check(f"insesgada a n={n}", abs(loo) < 0.05, f"loo={loo:+.3f} naive={nai:+.3f}")
    check(f"la version ingenua esta sesgada a n={n}", nai > loo, f"{nai:+.3f} > {loo:+.3f}")
# con direccion compartida debe ser alta
D = np.r_[3, 0, 0, 0, 0] + rng.normal(scale=0.3, size=(20, 5))
check("alta con direccion compartida", R.coherence_loo(D) > 0.9, f"{R.coherence_loo(D):.3f}")

print("alineamiento por persona")
A = np.r_[3, 0, 0, 0, 0] + rng.normal(scale=0.3, size=(15, 5))      # referencia concentrada
B = R.isotropic_null(rng.normal(size=(15, 5)), rng) * 3             # dispersa
D = np.vstack([A, B]); lb = np.r_[np.zeros(15), np.ones(15)].astype(int)
al = R.alignment_per_person(D, lb == 0)
check("separa grupos construidos", al[lb == 0].mean() - al[lb == 1].mean() > 0.5,
      f"{al[lb==0].mean():.2f} vs {al[lb==1].mean():.2f}")
obs, p, nullmean = R.perm_test_alignment(D, lb, rng, B=500)
check("la permutacion detecta la diferencia", p < 0.01, f"p={p:.4f}")
check("el nulo esta centrado en cero", abs(nullmean) < 0.08, f"media del nulo {nullmean:+.3f}")
# y NO debe detectar nada cuando los dos grupos vienen de la misma distribucion
D2 = np.r_[3, 0, 0, 0, 0] + rng.normal(scale=0.8, size=(30, 5)); lb2 = np.r_[np.zeros(15), np.ones(15)].astype(int)
_, p2, _ = R.perm_test_alignment(D2, lb2, rng, B=500)
check("no detecta diferencias inexistentes", p2 > 0.05, f"p={p2:.3f}")

print("nulo isotropo")
D = np.r_[5, 0, 0, 0, 0] + rng.normal(scale=0.2, size=(12, 5))
iso = R.isotropic_null(D, rng)
check("conserva las longitudes", np.allclose(np.linalg.norm(iso, axis=1), np.linalg.norm(D, axis=1)))
check("destruye la direccion compartida", R.coherence_loo(iso) < 0.4, f"{R.coherence_loo(iso):.3f}")

print("correccion de Benjamini-Hochberg")
p = np.r_[np.full(10, 1e-6), rng.uniform(size=990)]
check("monotona y acotada", np.all(np.diff(R.bh(p)[np.argsort(p)]) >= -1e-12) and R.bh(p).max() <= 1)
check("recupera los verdaderos positivos", (R.bh(p)[:10] < 0.05).all())

print("\n" + ("todo correcto" if not fails else f"FALLAN {len(fails)}: {fails}"))
sys.exit(1 if fails else 0)
