# Plan de referencias

El manuscrito tiene 19. El objetivo habitual en *Nature Metabolism* es ~50. Faltan sobre todo las
publicaciones primarias de los datos reanalizados, que son obligadas en un trabajo de reanálisis, y
las fuentes metodológicas. Abajo, tres bloques.

**Cómo leer el estado:** *confirmada* = verificada en esta sesión contra el depósito o el artículo.
*por verificar* = autor y año provienen del inventario del repositorio y la cita completa debe
comprobarse en PubMed antes del envío. No se incluye ninguna cita cuyos datos no consten.

---

## Bloque A — Publicaciones primarias de los datos (21 pendientes)

| # | Accession | Cita propuesta | Estado |
|---|---|---|---|
| A1 | GSE231509 | Makhnovskii P. et al. Dysregulation of early gene response to a mixed meal in skeletal muscle in obesity and type 2 diabetes. *Physiol. Genomics* (2023). doi:10.1152/physiolgenomics.00046.2023 | **confirmada** |
| A2 | GTEx v11 | GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science* 369, 1318–1330 (2020) | por verificar (edición v11 frente a v8) |
| A3 | GSE76895 | Solimena M. et al. (2018), islotes humanos, IMIDIA | por verificar |
| A4 | GSE164416 | Wigger L. et al. (2021), islotes, multiómica de T2D | por verificar |
| A5 | GSE50244 | Fadista J. et al. (2014), islotes, *PNAS* | por verificar |
| A6 | GSE50398 | Fadista J. et al. (2014), mismo estudio, plataforma de microarray | por verificar |
| A7 | GSE159984 | Marselli L. et al. (2020), islotes, perturbación ex vivo | por verificar |
| A8 | GSE18732 | Gallagher I. J. et al. (2010), músculo, *Genome Med.* | por verificar |
| A9 | GSE25462 | Jin W. et al. (2011), músculo | por verificar |
| A10 | GSE7146 | Rome S. et al. (2007), músculo, infusión de insulina | por verificar |
| A11 | GSE224310 | (2023), ejercicio agudo y entrenamiento, músculo y adiposo | **falta identificar** |
| A12 | GSE27951 | (2010), adiposo, delgados frente a obesos | **falta identificar** |
| A13 | GSE64567 | (2015), adiposo, por glucosa en ayunas | **falta identificar** |
| A14 | GSE135134 | METSIM, adiposo subcutáneo, n = 434 | por verificar (¿Civelek et al.?) |
| A15 | GSE20950 | Hardy O. T. et al. (2011), adiposo, depósitos pareados | por verificar |
| A16 | GSE26637 | Soronen J. et al. (2012), adiposo, clamp | por verificar |
| A17 | GSE59034 | Mitchell A. et al. (2015), adiposo, cirugía bariátrica | por verificar |
| A18 | GSE67297 | Hanssen M. J. W. et al. (2015), aclimatación al frío, *Nat. Med.* | por verificar |
| A19 | GSE66306 | (2015), monocitos, cirugía bariátrica | **falta identificar** |
| A20 | GSE156993 | (2020), PBMC | **falta identificar** |
| A21 | GSE21321 | (2010), sangre | **falta identificar** |

Ya citadas en el manuscrito y correctas: GSE22309 (ref 6), GSE9105 (7), GSE106800 (11),
GSE182117/182120 (12), GSE130646 (14), GSE81965/GSE63887 (15), GSE202295 (16), GSE198922 (17),
Rydén 2016 (18), GSE157988 (19).

---

## Bloque B — Métodos y estadística (16 propuestas, todas verificables)

| # | Cita | Dónde se usa |
|---|---|---|
| B1 | Benjamini Y. & Hochberg Y. *J. R. Stat. Soc. B* 57, 289–300 (1995) | corrección por comparaciones múltiples |
| B2 | Mardia K. V. & Jupp P. E. *Directional Statistics* (Wiley, 2000) | distribución von Mises–Fisher |
| B3 | Fisher N. I. *Statistical Analysis of Circular Data* (CUP, 1993) | alternativa a B2 para el tratamiento direccional |
| B4 | Banerjee A. et al. Clustering on the unit hypersphere using von Mises–Fisher distributions. *JMLR* 6, 1345–1382 (2005) | estimación de κ por máxima verosimilitud |
| B5 | Schönemann P. H. A generalized solution of the orthogonal Procrustes problem. *Psychometrika* 31, 1–10 (1966) | alineamiento ortogonal de subespacios |
| B6 | Langfelder P. & Horvath S. WGCNA. *BMC Bioinformatics* 9, 559 (2008) | redes de coexpresión |
| B7 | von Luxburg U. A tutorial on spectral clustering. *Stat. Comput.* 17, 395–416 (2007) | laplaciano normalizado y encaje espectral |
| B8 | Ritchie M. E. et al. limma. *Nucleic Acids Res.* 43, e47 (2015) | modelos lineales y eliminación de lote |
| B9 | Robinson M. D. et al. edgeR. *Bioinformatics* 26, 139–140 (2010) | filtrado y normalización de conteos |
| B10 | Love M. I. et al. DESeq2. *Genome Biol.* 15, 550 (2014) | transformación de conteos |
| B11 | Wolf F. A. et al. SCANPY. *Genome Biol.* 19, 15 (2018) | procesamiento unicelular (paso 19a) |
| B12 | Traag V. A. et al. From Louvain to Leiden. *Sci. Rep.* 9, 5233 (2019) | agrupamiento de células |
| B13 | Aitchison J. The statistical analysis of compositional data. *J. R. Stat. Soc. B* 44, 139–177 (1982) | log-ratio centrado de proporciones |
| B14 | Chikina M. et al. CellCODE. *Bioinformatics* 31, 1584–1591 (2015) | deconvolución de proporciones celulares |
| B15 | Schuirmann D. J. A comparison of the two one-sided tests procedure. *J. Pharmacokinet. Biopharm.* 15, 657–680 (1987) | TOST para equivalencia de magnitud |
| B16 | Lawson C. L. & Hanson R. J. *Solving Least Squares Problems* (1974) | mínimos cuadrados no negativos |

---

## Bloque C — Contexto que refuerza la discusión (12 propuestas)

| # | Cita | Para qué |
|---|---|---|
| C1 | DeFronzo R. A. et al. Glucose clamp technique. *Am. J. Physiol.* 237, E214–E223 (1979) | definición del valor M |
| C2 | Batista T. M. et al. A cell-autonomous signature of dysregulated protein phosphorylation underlies muscle insulin resistance in type 2 diabetes. *Cell Metab.* 32, 844–859 (2020) | firmas individuales; antecedente del argumento |
| C3 | Barrett E. J. et al. The vascular actions of insulin control its delivery to muscle. *Diabetologia* 52, 752–764 (2009) | entrega de insulina como paso limitante |
| C4 | Clerk L. H. et al. Microvascular recruitment in insulin resistance | perfusión microvascular |
| C5 | Williams I. M. et al. Transendothelial insulin transport is impaired in obesity | transporte transendotelial |
| C6 | Stenvers D. J. et al. Diurnal rhythms in the white adipose tissue transcriptome are disturbed in obese individuals with type 2 diabetes. *Diabetologia* 62, 704–716 (2019) | ritmo diurno alterado en T2D |
| C7 | Dyar K. A. et al. Transcriptional programming of lipid and amino acid metabolism by the skeletal muscle circadian clock. *PLoS Biol.* 16, e2005886 (2018) | salida del reloj en músculo |
| C8 | Perrin L. et al. Transcriptomic analyses reveal rhythmic and CLOCK-driven pathways in human skeletal muscle | reloj en músculo humano |
| C9 | Scheer F. A. J. L. et al. Adverse metabolic and cardiovascular consequences of circadian misalignment. *PNAS* 106, 4453–4458 (2009) | desalineamiento y metabolismo |
| C10 | Sylow L. et al. Exercise-stimulated glucose uptake: regulation and implications for glycaemic control. *Nat. Rev. Endocrinol.* 13, 133–148 (2017) | por qué el ejercicio es el contraste correcto |
| C11 | Smith G. I. et al. Muscle protein and metabolic responses in obesity | contexto de obesidad y músculo |
| C12 | Sarwar R. et al. / revisión reciente sobre heterogeneidad individual en T2D | apoyo al argumento de variabilidad entre personas |

C4, C5, C8, C11 y C12 están señaladas por tema y deben completarse con la referencia exacta.

---

## Recuento

| Bloque | Propuestas |
|---|---|
| Actuales | 19 |
| A, fuentes de datos | 21 |
| B, métodos | 16 |
| C, contexto | 12 |
| **Total si se aceptan todas** | **68** |

Para quedar en ~50: el bloque A es obligatorio completo (21), el B puede recortarse a los diez
realmente usados, y del C bastan seis u ocho. Eso deja 19 + 21 + 10 + 7 ≈ 57, y afinando el solape
del bloque C con las citas actuales, alrededor de 50.
