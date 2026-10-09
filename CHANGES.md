# Registro de cambios

## v9.1 — repositorio del artículo

- Raíz reducida a `README.md`, `THEORY.md` (marco conceptual + definiciones formales en un solo documento) y `CHANGES.md`.
- `docs/`: manuscrito, auditoría, tabla de replicación, exploración génica, hoja de ruta, borradores.
- README principal reescrito (el anterior había sido sobrescrito por error con `data/README.md`).

## v9 — reproducibilidad completa

- `python/lib/geo.py` (lectura de series matrix, anotación de plataformas, log-CPM) y `python/lib/response.py` (embedding, desplazamientos, coherencia leave-one-out, tests de permutación, interacción por gen, FDR).
- `python/analyses/01–08`: todos los análisis del manuscrito fuera del pipeline de descubrimiento, cada uno con su docstring y salidas en `results/`. `run_replication.sh` los encadena. `data/README.md` lista cada archivo de entrada y su origen.
- **Los números del manuscrito deben regenerarse desde estos scripts** (la coherencia es ahora leave-one-out en todas las tablas; los valores naive de borradores anteriores quedan sustituidos).
- Observación al regenerar: en GSE182120 la profundidad agregada respecto al estado sano muestra una asociación débil con el valor M (ρ ≈ −0.36) aunque ningún gen sobrevive FDR; el texto debe decir "señal difusa sin gen individual" en lugar de "sin información".

## v8 — auditoría

**Metodológicos**
- `landscape.py` [FIX-design]: la densidad se estimaba sobre los pacientes tal cual, de modo que la profundidad de cada valle reflejaba cuántos pacientes se reclutaron por estadio (32/15/36…), no la frecuencia de los estados. Ahora KDE, score, persistencia, umbral y barreras ponderan cada estadio a masa igual (`stage_weights`); `--no_balance` como sensibilidad (añadida a `run_all.sh`).
- `sheaf_coherence.py` [FIX-null-estimator]: el observado era una media de `reps` submuestras pero cada réplica del nulo era una sola submuestra → varianza del nulo inflada, test conservador. El nulo usa ahora el mismo estimador (media de `reps_null=2`). En sintéticos p(int max) pasa de 0.33 a 0.19 con la misma señal.
- `03`/`03b`: el test de permutación es válido (conserva tamaños de grupo) pero `obs_diff_*` sigue confundido por n; usar `*_sub` de `02` para el tamaño del efecto. Documentado.

**Rendimiento** (sin cambio numérico; verificado a 1e-9)
- `02`, `03`: Ḡ por eigen simétrico de D⁻¹ᐟ²WD⁻¹ᐟ² en vez de `expm()` (4× más rápido); R̄ con `only.values=TRUE` sin reconstruir L⁺ (3×, y O(p²) de memoria menos). Afecta a ~1 500 llamadas en 03.
- `sheaf_coherence.py`: matrices de expresión precomputadas por tejido (antes `expr.loc[genes]` en cada una de ~1 500 llamadas); `scipy.linalg.eigh(subset_by_index)` solo para los r+1 autovectores inferiores (~3×).

**Sin cambio, revisado y correcto**: bicor/β por tejido, submuestreo a n igual, nulos de configuración, Procrustes y referencia común, Hodge del puente, expansión de covariables GEO, `run_all.sh`.

## v7 — formalización y robustez

- `THEORY_FORMAL.md`: definiciones, 8 proposiciones y 1 teorema con demostraciones (haz ⇔ distancia de Grassmann; sesgo 1/n por Davis–Kahan y estimador extrapolado; correspondencia de Morse; conteo garantizado por el teorema de estabilidad; consistencia de la decisión de dos atractores; Kramers; Hodge del puente).
- `sheaf_coherence.py`: `E_half_n`, `E_n_extrapolated`, `sampling_bias_hat` (extrapolación de Richardson; en sintéticos b/n ≈ 0.08 frente a diferencias entre estadios ≈ 0.2).
- `landscape.py`: `persistence_threshold_guaranteed`, `stability_eps`, `n_basins_guaranteed` (bootstrap de Fasy et al. sobre el 80 % más denso); Fisher con separación mínima en θ (`[FIX-edge]`); persistencia también en el plano de la GMM (`n_basins_gmmplane`); expansión de `characteristics_ch1.N` de GEO y alias ampliados de covariables; Fisher para todos los controles disponibles (`fisher_peaks`); `--control`, `--no_covar`.
- Validado en sintéticos: biestable → 2 atractores nominales y garantizados; continuo → 0.

## v6 — renombrado

El repositorio pasa de `constructal_transcriptomics` a **`t2d_landscape`**. Se elimina la teoría constructal como marco (ver historial de versiones abajo para la justificación); el marco vigente es el de atractores / paisaje cuasi-potencial + haz celular (`THEORY.md`). Renombres: `02_networks_modularity_metrics.R`, `03_bootstrap_nulls_reference.R`; `*_network_metrics.tsv`; `CEI` → `NOI` (Network Organization Index); `approx_constructal_optimum` → `configuration_reference`; `results/optimum` → `results/reference` (`*_deviation_from_reference.tsv`).

# Historial — revisión v2

Todos los scripts R pasan `parse()`. Las métricas se verificaron numéricamente contra
sus definiciones exactas (R̄ = media de R_ij por par; EG = media de 1/d_ij). Los
scripts no se ejecutaron sobre datos reales (sin acceso a GEO en la revisión):
correr `01` primero y comprobar los avisos de conteo.

## Errores corregidos que afectaban resultados

| Tag | Dónde | Problema | Corrección |
|---|---|---|---|
| FIX-order | 01 | `hv_genes` quedaba en orden alfabético; `02/03` tomaban `[1:800]` → la "subred de alta varianza" eran los genes A–C | Ordenado por varianza media descendente |
| FIX-TRI-exact | 04, 04b | La "actualización rank-1 (Sherman-Morrison)" de tr(L⁺) es inválida: ΔL no es rango 1 y la fórmula no es la expansión de primer orden (error ~30× verificado) | Recalculo exacto por gen (EG+Hb en 04; R̄+Hb en 04b) |
| FIX-Gbar-exact | 02b, 03b | Truncar el espectro de D⁻½WD⁻½ a k autovectores: error 30–50 % (exp no decae en [−1,1]) | `eigen()` completo sobre la subred (segundos) |
| FIX-Rbar-v2 | todos | El "FIX-Rbar" anterior estaba invertido: Kf = p·tr(L⁺) ⇒ media por par = 2tr(L⁺)/(p−1) | Restaurado /(p−1) |
| FIX-EG | todos | EG multiplicado por 2 | Eliminado |
| FIX-beta-dataset | 02, 02b | β distinto por estado ⇒ W cambia de escala; las métricas comparaban β, no biología | β único por dataset, estimado en el estado sano |
| FIX-beta-min | 02b | Se elegía el β de máximo R² (sesgo a β altos) | Menor β con R² ≥ 0.80 y pendiente negativa (WGCNA) |
| FIX-geneset | 03b | El top-2000 por varianza se recalculaba en cada réplica bootstrap y grupo permutado | Fijo por dataset |
| FIX-regex | 05 | Los patrones capturaban `all_*.tsv` → filas duplicadas | Anclados `^GSE\d+_` |
| FIX-contrast | 04c | `combined_score` usaba el contraste de mayor |logFC| por gen (mejor de varios tests) | Contraste final-vs-referencia prefijado |
| FIX-dorng | 04 | `%dorng%` incondicional → fallo sin doRNG | Operador condicional `%dop%` |
| FIX-enrich | 04, 04b | Fallo de KEGG (sin internet) perdía también el GO | tryCatch separados, universo = subred |
| FIX-pval | 03, 03b | p = 0 posible | (b+1)/(B+1); se reporta nº de permutaciones distintas |
| FIX-subset | 02, 03 | bicor sobre ~5000 genes y recorte posterior (~40× más lento) | Recorte previo |
| FIX-disk | 02b | W densa de 14–22k genes por estado (1.7–3.8 GB) | Se guarda subred top-3000 + fuerzas completas |

## Errores corregidos en la 2ª pasada (tras ejecutar el pipeline completo sobre datos sintéticos)

| Tag | Dónde | Problema | Corrección |
|---|---|---|---|
| FIX-contrast-names | 04c | limma nombra los contrastes `"T2D - ND"`, no `"T2D_vs_ND"`; el contraste prefijado caía al fallback (IGT vs sano) | Nombres reales |
| FIX-enrich-optional | 04, 04b | `library(clusterProfiler/org.Hs.eg.db/enrichplot)` incondicional → abortaba sin ellos | `HAS_ENRICH`; enriquecimiento opcional |
| FIX-threads | 02 | `allowWGCNAThreads()` aborta con 1 core | Guardado |
| FIX-coroptions | 02 | `corOptions` string → avisos NA en WGCNA recientes | Lista |
| FIX-consistency | 02, 03 | 03 recalculaba la subred por su cuenta | 02 guarda `fixed_subnetwork_genes.rds`; 03 lo lee |
| FIX-dplyr | 05 | `cur_data()` deprecado | `pick()` |
| — | 02b, 03b | Mensajes obsoletos sobre truncación | Actualizados |

## Banderas de línea de comandos

Todos los scripts aceptan `--nombre=valor` (prioridad sobre la variable de entorno homónima, que se conserva) y los accessions como argumentos posicionales. `--help` muestra el uso.

| Bandera | Scripts | Defecto | Env equivalente |
|---|---|---|---|
| `--workers` | 03, 04, 02b, 03b, 04b | 40 (48 en 04b) | `N_WORKERS` |
| `--p_sub` | 02 | 800 | `P_SUB` |
| `--n_sub_reps` | 02 | 20 | `N_SUB_REPS` |
| `--n_boot` | 03, 03b | 100 | `N_BOOT` |
| `--n_perm` | 03, 03b | 1000 | `N_PERM` |
| `--max_genes_node`, `--max_genes_tri`, `--n_top_enrich` | 04b | 3000, 1500, 200 | `MAX_GENES_*`, `N_TOP_ENRICH` |

Prueba rápida antes de la corrida completa:
```bash
Rscript R_scripts_sub_networks/02_networks_modularity_metrics.R --p_sub=300 --n_sub_reps=3 GSE27951
Rscript R_scripts_sub_networks/03_bootstrap_nulls_reference.R --workers=8 --n_boot=5 --n_perm=20 GSE27951
```

## v5 — reformulación teórica y arquitectura

- `THEORY.md`: marco de atractores / paisaje cuasi-potencial (U, J) + haz celular; predicciones P0–P5 prefijadas; el marco anterior pasa a descripción secundaria.
- `python/landscape.py`: potencial condicional (covariables centradas por estadio), puntos críticos por score, persistencia de subnivel H₀ con umbral ln 2 y estabilidad en escala, ΔBIC, asimetría de barreras (bootstrap), índice de transición crítica I_c, Fisher–Rao a lo largo de HbA1c/glucosa si existe, puente de Schrödinger (Sinkhorn) + descomposición de Hodge de la deriva (‖J‖/‖∇U‖). Regla de decisión conjunta `two_attractors`. Validado con sintéticos biestable (detecta) y continuo (no detecta).
- Hallazgo metodológico: regresar BMI crudo (colineal con el estadio) borra la señal de los valles; se regresa centrado por estadio. `--no_covar` como sensibilidad.
- `run_all.sh` con `--quick`. Banderas CLI en todos los scripts R.
- Enriquecimiento opcional (sin clusterProfiler el pipeline corre igual).

## Nuevo

- **Métricas a n igual** (`*_sub`, `*_sub_sd`): media sobre 20 submuestras al n mínimo del dataset. El NOI se calcula sobre ellas. Motivo: con n pequeño |cor| es mayor ⇒ red más densa ⇒ EG↑, R̄↓; los grupos IGT son los más pequeños.
- **Métricas relativas al nulo de configuración** (`*_rel`, `NOI_rel`): W₀ = k_i k_j / 2m preserva fuerzas y elimina estructura. Separa estructura de densidad.
- **Análisis de coherencia cross-tejido (haz celular)**: `python/sheaf_coherence.py`. Lee `export_sheaf/` (lo genera `01`). Ver docstring.
- **Simulaciones** en `python/simulation/` (validación sintética, robustez, potencia con n reales).

## Reinterpretaciones (sin cambio de código, sí de texto)

- `configuration_reference`: W* ∝ (k_i k_j)^{1/α} es, salvo exponente, el modelo de configuración. "Desviación del óptimo" = cantidad de estructura, no ineficiencia.
- `obs_diff_CEI` → `obs_diff_EGHb`: la suma EG+Hb de las permutaciones no es el NOI z-scoreado.
- PTI/IRI/TRI/KO son sensibilidad topológica de la red de correlación al nodo, no evidencia causal de rescate/knock-out.

## Limitaciones que quedan (documentar en el artículo)

- Sin nulo a nivel de gen para PTI/IRI (requiere reconstruir redes bajo permutación; coste alto).
- GSE15653 (n = 5/4/9): bicor con n = 4 es ruido; solo 126 permutaciones distintas. Excluir de la pregunta IGT.
- Diseño transversal: "critical slowing down" no es medible; usar el marco DNB / haz para "coherencia".
- La energía del haz también sube con pérdida de modularidad (ruido), no solo con incoherencia. Ver curva de potencia.

## Orden de ejecución

```bash
Rscript R_scripts_preprocessing/01_download_qc_preprocess.R        # genera export_sheaf/
Rscript R_scripts_sub_networks/02_networks_modularity_metrics.R
Rscript R_scripts_sub_networks/03_bootstrap_nulls_reference.R
Rscript R_scripts_sub_networks/04_gene_drivers_and_enrichment.R
Rscript R_scripts_full_networks/02b_full_networks_biology.R
Rscript R_scripts_full_networks/03b_full_bootstrap_and_nulls.R
Rscript R_scripts_full_networks/04b_full_gene_drivers_and_enrichment.R
Rscript R_scripts_full_networks/04c_expression_differential.R
Rscript R_scripts_sub_networks/05_cross_dataset_summary.R
# Coherencia cross-tejido (python >= 3.9, numpy, scipy, pandas)
python python/sheaf_coherence.py --export_dir export_sheaf --out results/sheaf \
    --tissues GSE76895 GSE18732 GSE27951 --n_genes 800 --reps 30 --B 500
# Sensibilidad con hígado (solo contraste T2D vs sano):
python python/sheaf_coherence.py --tissues GSE76895 GSE18732 GSE27951 GSE15653 --min_n 4 --out results/sheaf_with_liver --no_loto
# Probar el módulo sin datos reales:
python python/simulation/make_fake_export.py && python python/sheaf_coherence.py --export_dir export_sheaf_FAKE --out results/sheaf_FAKE --n_genes 300 --reps 6 --B 40 --no_loto --no_sens
```

## Antes de confiar en los resultados

1. `results/qc/dataset_summary_all.tsv`: conteos = 83 / 118 / 18 / 33 y sin `NA`. Si no, revisar `map_conditions()`.
2. Log de `01`: `Common genes across all datasets` ≥ 3000; si no, el CDF custom de GSE18732 no mapeó.
3. `results/sheaf/main_permutation_tests.tsv`: estadísticos `int_max`, `T2D_gt_healthy`, `delta_int_gt_T2D`.
