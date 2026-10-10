# Registro de cambios

## v9.18 — correcciones de la revisión editorial

- **Fig 1f con el nulo real**: `03b` exporta las 200 permutaciones (`sheaf_null_draws.tsv`) y la figura las grafica; z pasa de −26 (100 permutaciones) a **−29.7** (200).
- **Sin fuga en sangre→tejido**: residualización, selección, escalado y ambos PCA dentro de cada fold; R² 0.279 / 0.229 / 0.217 (antes 0.296 / 0.233 / 0.228).
- **Welch** real en `05` y `09`; permutaciones 2 000 / 3 000 según la leyenda.
- **Nulo isótropo** (`isotropic_null`) en la simulación de potencia, en lugar del cambio de signo, que conservaba los ejes.
- **TOST de magnitud**: ninguna comparación es equivalente con margen 0.5 log2, así que "magnitud conservada" pasa a "no se detectó reducción".
- **Ejercicio**: `exercise_group_tests.tsv` con Δcoherencia, cosenos y P reales (músculo en recuperación P = 0.70, coseno 0.84).
- **Fig 2f** sustituye los "none" por descripciones por cohorte; `docs/REVIEW_RESPONSE.md` lleva el estado de los 13 puntos.

## v9.17 — especificidad, confusión y potencia

- `11_exercise_specificity.py`: el músculo diabético responde al ejercicio agudo con la misma coherencia que el sano (0.53 vs 0.45, P = 0.69, coseno 0.84; GSE202295), y lo mismo en adiposo (GSE198922). La pérdida es específica de la insulina. Panel 5j.
- `12_power_and_confounding.py`: el estadístico leave-one-out es insesgado a cualquier n (el ingenuo da 0.41 con n = 4); el lote de hibridación genera por sí solo coherencias de 0.75–0.87, por lo que el análisis de mismo lote pasa a ser primario; potencia por tamaño y efecto. ED Fig. 3.
- GSE159984 en reposo (85 donantes) añadido a Fig 2g,h; su respuesta ex vivo en ED Fig. 2.
- Manuscrito actualizado: resumen con especificidad, análisis primario de mismo lote, limitaciones reescritas sobre ED Fig. 3.

## v9.16 — manuscrito completo

- `docs/MANUSCRIPT_FULL.md`: manuscrito en formato Nature (resumen, introducción, resultados, discusión y limitaciones, métodos escritos por propósito, referencias, leyendas y material suplementario), compilado a PDF con las figuras al final.
- Material suplementario: nota sobre por qué el haz y qué añade frente a comparaciones más simples; nota sobre el órgano ausente con la propuesta de GSE159984 y GSE53949 para islote; cinco tablas suplementarias generadas desde los resultados.

## v9.15 — título

- Título definitivo: *Insulin resistance disrupts transcriptional coordination in skeletal muscle*. Aplicado en `docs/MANUSCRIPT.md` y `README.md`; se retira el provisional.

## v9.14 — manuscrito

- `docs/MANUSCRIPT.md`: borrador completo en formato Analysis de *Nature Metabolism* (resumen, introducción, seis secciones de resultados ancladas a las figuras, discusión, limitaciones, métodos resumidos, disponibilidad de datos y código), con el tono calibrado de las leyendas revisadas. El borrador anterior en formato Cell Press queda como `docs/DRAFT_CellMetab_format.md`.

## v9.13 — lienzo mas ancho

- Ancho de figura configurable con `T2D_FIG_WIDTH_MM` (250 mm por defecto para la redacción; 180 mm es el doble columna de Nature y basta poner `T2D_FIG_WIDTH_MM=180` para la versión de envío). Los tipos de letra son absolutos, así que ensanchar el lienzo da aire a los paneles sin reducir la legibilidad.

## v9.12 — las redes en la figura

- `python/analyses/10_network_panel.py`: red de coexpresión de ejemplo (disposición por fuerzas y módulos) y copia de las métricas por estadio a n igual; encadenado en `run_replication.sh`.
- Fig 1a–c: la red sobre la que se construye el haz y las métricas de organización ($E_G$, $\bar{R}$) por estadio y órgano. Paneles renumerados a–l y leyendas actualizadas.

## v9.11 — repositorio al día

- `.gitignore` (datos descargados, `results/`, `figures/`, caches) y mapa del repo actualizado en `README.md`: `docs/` desglosado, pasos 01–09, generación de figuras y variables de entorno (`T2D_RAW`, `T2D_EXPORT`, `T2D_RES`, `T2D_OUT`, `T2D_FIG`).
- `run_replication.sh`: un log por paso con el nombre del script (antes todos se llamaban `python.log`), exporta `T2D_EXPORT` y encadena la generación de figuras al final.
- `02_run_replication_landscape.sh`: ruta de `landscape.py` corregida y respeta `T2D_EXPORT`; `--quick` ya no duplica argumentos.
- `09_classical_de.py` falla con mensaje claro si no hay cohortes exportadas en lugar de escribir tablas vacías; `make_figures.py` tolera tablas ausentes o vacías.
- Título del repositorio alineado con el del manuscrito y nota sobre el reenfoque a Nature Metabolism.

## v9.10 — expresión diferencial canónica y revisión de maquetación

- `python/analyses/09_classical_de.py`: DE clásica (T2D/obeso vs control, en reposo) y enriquecimiento por conjuntos curados en cuatro cohortes; encadenado en `run_replication.sh`.
- Fig 2g,h: recuento de genes a FDR < 0.1 por cohorte y mapa de calor de enriquecimiento (islote y adiposo recuperan las firmas esperadas; músculo da 2 genes). Fig 1i: genes del eje sistémico compartido, para comparar drivers canónicos y propios.
- Revisión visual de las seis figuras: leyendas fuera de los datos, letras de panel separadas, etiquetas de genes sin solapes.

## v9.9 — sin títulos de panel

- Los paneles ya no llevan título ni encabezados de fila dentro de la figura (convención Nature Portfolio): el contenido se describe en `docs/FIGURE_LEGENDS.md`, con un pie por figura y una frase por panel. Los `set_title` del script quedan como documentación y se dibujan solo con `T2D_PANEL_TITLES=1` para revisión interna. Alturas ajustadas al espacio liberado; etiquetas de grupo (IS/IR/T2D) dentro de los paneles polares y aluviales.

## v9.8 — estilo tipográfico

- Jerarquía como en figuras de referencia: 3 paneles por fila, letras de panel 13 pt, títulos que afirman el resultado, ejes 8–9.5 pt, líneas 2 pt, marcadores grandes, etiquetas directas, paleta de dos colores + negro. Fig 1 y Fig 4 reducidas a 2×3; lo secundario a Extended Data (ED_Fig1 robustez y sangre; ED_Fig4 artefacto de biopsia; ED_Fig5 reloj, células y destino del programa).

## v9.7 — estética de los paneles

- Rosas de los vientos para la dirección de respuesta por persona (Fig 5d), relojes polares para la salida del reloj por grupo (Fig 5f), diagramas aluviales para el destino del programa (Fig 5i), coordenadas paralelas por donante frente a barajado (Fig 1d); helpers en `python/lib/response.py` (`rose`, `chord`, `alluvial`). `03b` escribe `donor_scores.tsv`.

## v9.6 — Fig 2 «dónde está la enfermedad»

- Nueva figura entre el estado sistémico y el músculo: búsqueda órgano por órgano en reposo y bajo insulina, con el mapa de evidencia que justifica el músculo. Seis figuras principales.

## v9.5 — figuras de la historia final

- Cinco figuras (estado sistémico y continuo → basal ciego → respuesta sana → fragmentación → adiposo e intervenciones); fuera de figura lo que no está en la historia.
- Puntos individuales con media ± s.d., corchetes de P, nulos como violines/bandas, leyendas fuera de los datos, letras de panel colocadas por geometría tras el layout, títulos de hasta tres líneas alineados a la izquierda.

## v9.4 — estilo Nature Portfolio

- Figuras a 180 mm de ancho (doble columna), Arial/Helvetica (Liberation Sans como sustituto), texto 5.5–6.5 pt, letras de panel 8 pt en negrita, líneas ≥ 0.5 pt, `constrained_layout` para evitar solapes; fuentes embebidas (Type 42) en PDF.

## v9.3 — figuras multipanel

- Seis figuras siguiendo `docs/FIGURE_PLAN.md` (órganos → continuo → basal → respuesta sana → fragmentación → adiposo), cada una con esquema, paneles de prueba y comprobaciones; nulos siempre como distribuciones; conjuntos de genes curados (SETS en `make_figures.py`) con nulo de conjuntos aleatorios en lugar de GSEA externo.

## v9.2 — figuras

- `python/figures/make_figures.py`: Fig 1–5 generadas solo desde `results/` (más `data/raw` para los vectores de GSE22309). `figures/` no se versiona; se regenera.
- Título provisional del manuscrito: *Insulin resistance fragments a coordinated, clock-coupled transcriptional response in humans*.

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
