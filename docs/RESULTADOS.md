# Resultados consolidados

Estado a 11 de octubre de 2026. Cada afirmación lleva su cohorte, su tamaño y su control.
Se separa lo verificado en personas vivas de lo observado solo en GTEx, donde la circunstancia
de la muerte influye de forma demostrada.

---

## 1. El estado en reposo no informa sobre la sensibilidad a la insulina

En 49 personas con clamp en todas (GSE182120), correlacionar la expresión basal con el valor M
sin ajustar devuelve **739 genes a FDR < 0.1**, con 85 de 211 genes mitocondriales entre ellos,
reproduciendo el enriquecimiento publicado para esa cohorte. La asociación **desaparece por
completo al ajustar por edad o por índice de masa corporal, cualquiera de los dos por separado**.
El IMC por sí solo no tiene ningún gen asociado, de modo que la señal procede de la varianza
compartida entre adiposidad, edad y sensibilidad, no de una firma de adiposidad.

Con correlación parcial y los dos lados residualizados: **0 genes de 21 595 a FDR < 0.1**, y 231
a P < 0.01 donde se esperan 215. Un diseño transversal en reposo no distingue un correlato de la
sensibilidad de un correlato de sus determinantes.

---

## 2. La respuesta a la insulina tiene un núcleo universal y una periferia individual

Dos criterios distintos, que antes mezclé y ahora separo. **Unanimidad estricta** (signo igual en
más del 90 % de las personas y |t| > 6): en GSE22309 (55 personas, >8 000 genes) la cumplen
**exactamente dos genes, RRAD y DDIT3**. **Respuesta altamente consistente** (signo igual en más del
80 %): la cumplen RRAD, TXNIP y PDK4, que son los que replican entre cohortes. En GSE13070
(68 personas, otra plataforma):

| Gen | GSE22309 | GSE13070 |
|---|---|---|
| RRAD | t = 9.6, 93 % | t = 10.4, 85 % |
| TXNIP | — | t = −13.3, 91 % |
| PDK4 | — | t = −11.0, 82 % |

El resto de la respuesta difiere entre personas. La heterogeneidad no es una propiedad de la
resistencia a la insulina: es la propiedad normal de la respuesta **en las condiciones estudiadas**,
también en personas sanas. Queda por evaluar si una predicción multigénica fuera de muestra
recupera señal donde el análisis gen a gen no la encuentra.

### Compartido entre tejidos, con componente propio

En adiposo bajo clamp (GSE26637, 10 personas) TXNIP, PDK4, HK2, IRS2 y ANGPTL4 responden **10/10**,
mientras **RRAD no se mueve** (5/10). El núcleo es compartido; RRAD es específico de músculo.

### Dos capas temporales

Con comida mixta a 1 h (GSE231509, 21 personas): ANGPTL4 responde 21/21 (t = 9.5), HK2 81 %,
IRS2 71 %; **TXNIP, PDK4 y RRAD apenas se mueven** (24 %, 33 %, 14 %) pese a ser los más fuertes
a 3–5 h. Capa rápida y capa lenta.

### Salida del reloj: cambio agudo sin cambio detectable en el oscilador

A 1 h de la comida: NR1D2 21/21 (t = −9.5), DBP y CIART 90 % (t = −8.3 y −8.8), PER1/PER3/TEF/HLF
con |t| 3.5–5.6. **No se detecta cambio en ARNTL (BMAL1) ni BHLHE40** (t = −0.71 y −0.08). El mismo
patrón aparece en adiposo bajo clamp (DBP, NR1D2, PER2 todos 10/10).

Dos cautelas. Una hora y dos genes sin cambio **no excluyen efectos de tiempo ni desplazamientos de
fase**: haría falta muestreo a varios tiempos con control horario. Y la ausencia de cambio agudo en
BMAL1 no demuestra un oscilador inalterado — de hecho el experimento de hidrocortisona de la sección
3bis **sí mueve ARNTL** (+0.78, P = 0.007). Lo que estos datos sostienen es que la insulina aguda a
1 h modifica la salida sin que se detecte cambio en el oscilador, no que el oscilador sea insensible.
Hay además literatura previa de que la insulina modifica PER1 y PER2 y puede desplazar la fase del
reloj adiposo.

---

## 3. Un eje individual sistémico de tono catabólico

### Qué es

En GTEx (647 donantes con músculo y adiposo), el eje que comparte una persona entre sus órganos
está dominado por la firma canónica de glucocorticoides: FKBP5 +0.84/+0.76, GLUL +0.79/+0.75,
MT1X +0.71/+0.65, ZBTB16 +0.67/+0.65. Frente a un fondo de 0.093:

| Conjunto | n | \|r\| concordante | P |
|---|---|---|---|
| Glucocorticoides (literatura) | 16 | 0.478 | 0.0002 |
| Núcleo de insulina | 6 | 0.340 | 0.0002 |
| Salida del reloj | 9 | 0.239 | 0.0012 |
| Oscilador central | 4 | 0.133 | 0.22 |
| Oxidativo / PGC-1α | 10 | 0.042 | 0.93 |

**Confirmado en personas vivas**: en 28 mujeres con músculo y adiposo (GSE41168), con un fondo de
+0.024, glucocorticoides 0.304 (P = 0.0002), núcleo de insulina 0.252 (P = 0.004), salida del
reloj 0.238 (P = 0.034), **oscilador −0.090** (P = 0.86), oxidativo 0.102 (n.s.). Se reproduce la
jerarquía completa, incluida la separación entre salida y oscilador.

### Glucocorticoides e insulina son el mismo eje; el reloj no

Correlación entre puntuaciones por donante: glucocorticoide–insulina **0.75**, metalotioneínas
0.66 y 0.59, **salida del reloj −0.14 y −0.24**, genético −0.22, **oxidativo −0.04 y 0.01**.

### Qué programa es

Enriquecimiento formal (MSigDB v2026.1.Hs, con conjuntos aleatorios del mismo tamaño como control):

- **Con el eje**: proteasas de ubiquitina, degradación mediada por ubiquitina, respuesta celular
  al estrés químico, p53, dianas de MYC. A mano: atrofia muscular (FBXO32, TRIM63) r = 0.447;
  autofagia/FOXO 0.416; metalotioneínas 0.441; lipólisis 0.211.
- **Contra el eje**: organización de matriz extracelular y colágeno, **transporte electrónico
  respiratorio**, contracción muscular, interferón gamma, integrinas.
- **Fuera**: inflamación (P = 0.25), estrés de retículo, angiogénesis, respuesta inmediata.

Es un programa **catabólico, antimiogénico y antimatriz**. El conjunto oxidativo **no es
homogéneo**: PPARGC1A (+0.216) y ESRRA (+0.118) son positivos mientras los componentes de la cadena
respiratoria son negativos (ATP5F1A −0.179, NDUFA4 −0.139, UQCRC1 −0.105). «Ortogonal» describe la
media del conjunto (−0.020), no el comportamiento de cada gen. Los conjuntos de inflamación y los
mitocondriales muestran relaciones distintas con este eje — nula la primera, inversa la segunda —
lo que **no excluye esos mecanismos de la enfermedad**, solo los sitúa fuera de este componente.

### El regulador

DoRothEA sobre 193 regulones da **FOXO3 como el más asociado** (r = 0.66 con el eje de
glucocorticoides y 0.67 con el de insulina), seguido de receptor de andrógenos (0.62), CEBPB
(0.55) y BCL6 en negativo (−0.53). Las dianas de factores de transcripción de MSigDB, construidas
desde sitios de unión y por tanto independientes de la expresión, apuntan en la misma dirección
pero **con soporte más débil del que informé primero** (ver «Errores corregidos»): FOXO3_01
z = +2.61, P = 0.0067, **FDR = 0.056**; CEBP_C z = +2.33, P = 0.0147, **FDR = 0.095**. Sobre 853
conjuntos, ninguno de los dos supera el 5 % tras corrección múltiple. Sí lo supera **MEF2 en
negativo** (z = −4.85, FDR = 0.015), y el motivo del elemento de respuesta a glucocorticoides no
aporta nada (GRE_C z = +0.86, FDR = 0.62), lo que era esperable porque esos conjuntos de motivo
son ruidosos.

La evidencia sobre FOXO3 descansa por tanto en la convergencia de tres fuentes parciales —DoRothEA
(r = 0.66), las dianas por motivo (P nominal 0.0067) y la literatura— y no en un único resultado
significativo tras corrección.

La literatura respalda el mecanismo pieza por pieza: FOXO3 media la comunicación entre
IRS-1/PI3K/Akt y la señal de glucocorticoides y coordina ubiquitina y atrogina-1; el promotor
humano de PDK4 requiere unión de FOXO3 y FOXO1 para su inducción máxima por glucocorticoides;
KLF15 induce atrogina-1, MuRF1, FOXO1 y FOXO3. **Lo nuevo no es el circuito, es que sea un componente
prominente de la variación entre personas medida en dos órganos a la vez.** No se ha cuantificado
qué fracción de la covariación intertisular explica frente a otros componentes bajo el mismo
procedimiento, de modo que «el mayor» no está demostrado: puede dominar la covariación entre órganos
sin dominar la variación total de cada órgano.

La actividad inferida de FOXO3 y el enriquecimiento de sus dianas **no demuestran causalidad**. El
enriquecimiento de anotaciones de unión sobre una lista seleccionada por expresión es evidencia
externa, no una validación experimental independiente.

---

## 3bis. El eje responde a una perturbación glucocorticoide controlada

**GSE148640**: 10 pacientes con insuficiencia suprarrenal, adiposo subcutáneo, **hidrocortisona
frente a placebo, diseño intrasujeto, personas vivas**. Se aplica la puntuación **congelada** (el
conjunto de glucocorticoides procede de la literatura, no de estos datos ni de GTEx) sin reajustar
nada, con conjuntos aleatorios del mismo tamaño como control:

| Conjunto | n | Δ (cortisol − placebo) | P vs azar |
|---|---|---|---|
| Glucocorticoides (congelada) | 9 | **+0.582** | 0.0002 |
| Núcleo de insulina | 6 | **+0.482** | 0.0002 |
| Oscilador | 5 | +0.166 | 0.012 |
| Atrofia / FOXO | 4 | +0.157 | 0.029 |
| Salida del reloj | 9 | +0.026 | 0.56 |
| Oxidativo | 5 | −0.041 | 0.48 |

Gen a gen: ZBTB16 +2.05 (P = 4·10⁻⁸, el gen más inducido del experimento), PER1 +1.50 (8·10⁻⁸),
TXNIP +1.07 (2.8·10⁻⁶), FKBP5 +1.04 (2.5·10⁻³), PDK4 +0.99 (1.8·10⁻⁴), GLUL +0.65, CEBPD +0.51.
También aparece CYP4B1, que era de los primeros genes del eje en GTEx.

Esto convierte «firma compatible con glucocorticoides» en **«programa que responde a una
perturbación glucocorticoide en humanos vivos»**. Tres matices que el propio experimento impone:

- **ARNTL sube** (+0.78, P = 0.007). El cortisol sí mueve el oscilador.
- **La atrofia apenas se mueve**: FBXO32 +0.03, TRIM63 +0.13, ninguno significativo. El programa de
  atrofia se asocia al eje en la variación entre personas pero no responde a una dosis aguda; la
  atrofia es un proceso de días. La etiqueta «catabólico» describe la asociación, no una respuesta
  demostrada a corto plazo.
- **La salida del reloj como conjunto no se mueve** (+0.03) aunque PER1 se dispare: el cortisol toca
  PER1 específicamente, no el conjunto en bloque.

No establece que el cortisol explique toda la variación basal entre personas, solo que el programa
es sensible a la perturbación.

**GSE17992** (dexametasona 4 días, músculo y adiposo) se evaluó y **se descartó por problemas de
procesamiento del depósito**, no por su resultado: identificadores de plataforma que no casan con
los archivos crudos, valores centinela de ±10⁵ mezclados con los datos, y orientación de intercambio
de colorante no reconstruible desde los metadatos. En la matriz depositada FKBP5 no se mueve tras
4 días de dexametasona, lo que es implausible y delata el problema.

---

## 4. El eje a través de los tejidos

En 20 tejidos de GTEx, el eje existe con coherencia interna 0.27–0.58 y concordancia 0.6–0.95 con
la referencia. **Dos tejidos quedan fuera, y son los esperables**: la suprarrenal, que produce el
cortisol, no correlaciona con el eje (r = 0.078); el hígado, con la expresión de HSD11B1 más alta
con diferencia (8.06 frente a 1.3–4.4), queda desacoplado (coherencia 0.190). La correlación entre
coherencia del eje y expresión del receptor NR3C1 es +0.34.

El **nervio tibial encabeza la lista** (0.581) con la mayor expresión del receptor.

El desacoplamiento hepático es **una hipótesis compatible con la expresión observada de HSD11B1**,
no una relación demostrada: la correlación entre coherencia del eje y HSD11B1 a través de tejidos es
de solo −0.16, con 20 puntos.

---

## 5. Una biopsia de grasa sirve de proxy sistémico; la sangre no

R² validado fuera de muestra, GTEx:

| Desde → hacia | n | R² |
|---|---|---|
| Adiposo → nervio tibial | 542 | 0.82 (ver abajo: inflado) |
| Adiposo → músculo | 647 | 0.64 |
| Músculo → adiposo | 647 | 0.64 |
| Adiposo → hígado | 201 | 0.36 |
| **Sangre → músculo** | 710 | **0.14** |
| Sangre → cualquier tejido sólido | — | 0.14–0.17 |

**Validación sin fuga de información.** Todo el aprendizaje ocurre dentro de los pliegues de
entrenamiento — ajuste de covariables, estandarización y modelo — y los dos tejidos de un donante
van siempre juntos. Predictor y resultado usan **genes disjuntos**: 9 de la firma congelada en
adiposo frente a 8 genes distintos en músculo.

| Configuración | R² validado |
|---|---|
| Adiposo (firma congelada) → músculo (genes distintos) | **0.601** |
| Ídem, ajustando además circunstancia de muerte | **0.505** |
| Solo covariables (edad, sexo, Hardy) | 0.237 |
| Ídem sin ajuste técnico (RIN, isquemia) | 0.651 |

El estado agónico aporta alrededor de un sexto de la predicción; las dos terceras partes restantes
son independientes de él. **El 0.82 de adiposo a nervio de la tabla anterior queda desautorizado**:
se obtuvo con la misma firma a ambos lados y sin estos controles. La cifra defendible es **R² = 0.50
tras ajustar circunstancia de muerte, con genes disjuntos**, referida a la covariación dentro de GTEx.

**Lo que esto permite afirmar y lo que no.** La formulación que sostienen los datos es que *el
adiposo predice el componente compartido en otros tejidos de GTEx*. «Una biopsia de grasa sirve de
proxy sistémico» exigiría demostrar la transferencia en personas vivas, donde la única medida
disponible (28 mujeres, GSE41168) da una concordancia de conjunto de 0.30, bastante menor que el
0.77 de GTEx. Pendientes además: validación con la selección de genes **dentro del entrenamiento**,
comparación frente a un modelo de covariables, predicción tras ajustar la circunstancia de la
muerte, y R² con genes distintos en predictor y resultado.

---

## 6. Demografía y circunstancia de la muerte

- **Edad: sin efecto** (ρ = 0.016; 0.004 ajustando sexo y Hardy).
- **Sexo**: mujeres +0.124 frente a hombres −0.060 (P = 0.001), coherente con que el conjunto
  Hallmark de respuesta a andrógenos sea el más positivo.
- **Circunstancia de la muerte**: gradiente monótono de más de una desviación típica, de −0.72
  en muerte violenta rápida a **+0.46 en muerte lenta**, con los casos en ventilador (la mitad de
  la cohorte, y receptores habituales de corticoides) en +0.14. **El estrés agónico amplifica el
  eje en GTEx**, aunque el eje existe en vivos (sección 3).

---

## 7. Tono catabólico y diabetes

En GSE22309 el tono basal separa a los diabéticos: IS −0.14, IR −0.26, **T2D +0.54**
(Kruskal P = 0.001). Es un salto asociado al diagnóstico, no un gradiente con la sensibilidad.

El tono **no** explica la coordinación de la respuesta (ρ = −0.027 con el alineamiento individual;
conjuntos aleatorios dan |ρ| medio 0.204, es decir, más).

**Esta asociación con T2D necesita replicación fuera de GSE22309**, cohorte en la que el lote y el
grupo están casi confundidos (sección «Lo que no se sostuvo»). La puntuación basal requiere su
propio control de lote, comparaciones dentro de un mismo lote, ajuste por edad, sexo, IMC y
tratamiento, y efectos con intervalos de confianza además de P. Si no replica, la diabetes debe
quedar como contexto y resultado exploratorio, no como eje clínico del trabajo.

---

## 8. Implicación terapéutica

El eje ya tiene farmacología: los inhibidores de 11β-HSD1 reducen el cortisol generado dentro del
tejido. INCB-13739 mejoró la sensibilidad en fase II en pacientes que fallaron con metformina;
el carbenoxolone redujo la producción hepática de glucosa; en músculo la enzima está sobreexpresada
en miotubos de pacientes con resistencia y diabetes.

Precisiones necesarias sobre esa evidencia: INCB-13739 mejoró HbA1c, glucosa y HOMA-IR, lo que **no
equivale** a mejora de sensibilidad medida directamente por clamp. El ensayo de AZD4017 con
prednisolona no alcanzó su desenlace de disposición de glucosa. El estudio de cicatrización fue un
ensayo distinto, con la actividad cutánea de 11β-HSD1 como desenlace primario. **No deben
presentarse como una secuencia demostrada de reorientación de un programa.**

La propuesta defendible a partir de estos datos no es que los ensayos midieran mal — nuestros datos
no permiten explicar por qué fallaron ni demostrar un beneficio clínico oculto — sino que **una
lectura transcripcional tisular podría complementar los desenlaces metabólicos** para evaluar
respuesta farmacodinámica, heterogeneidad entre pacientes y compromiso efectivo del programa
biológico. El eje es accesible por biopsia de grasa subcutánea y responde a una perturbación
glucocorticoide controlada (sección 3bis).

---

## Lo que no se sostuvo

Por orden en que se intentaron, con el control que los descartó:

- **Pérdida de coordinación de la respuesta a la insulina como hallazgo central.** Dependía de
  GSE22309: quitándola, b₁ pasa de −0.58 a −0.215 y la interacción de P = 0.0023 a 0.15. GSE13070
  (68 pares, clamp de 5 h) no replica (P = 0.81) y GSE24215 (resistencia inducida por reposo en
  las mismas personas) va en sentido contrario (P = 0.55).
- **Eje de empalme CLK1/SRSF**: set-level P = 0.39 y 0.077.
- **Transporte transendotelial por caveolas**: no replica en GSE202295 y la correlación no supera
  a conjuntos aleatorios (P = 0.32).
- **Actividad de factores de transcripción como marcador de grupo**: 39 regulones a FDR < 0.05 con
  todas las parejas, 0 restringiendo a mismo lote.
- **Eje de glucocorticoides en el destino del nervio diabético**: P = 0.38 (GSE148059, n = 77).
- **Edad como determinante del tono catabólico**: nulo.

---

## Errores corregidos en auditoría

- **Benjamini-Hochberg mal implementado** en el script auxiliar de enriquecimiento de esta sesión:
  tomaba el mínimo acumulado hacia delante en lugar de hacia atrás, lo que propaga el FDR del P más
  pequeño a todas las filas siguientes. Efecto: en dianas de factores de transcripción, 842
  conjuntos parecían pasar FDR < 0.05 cuando pasan **91**; en Reactome, 1302 frente a **288**.
  Los P nominales no cambian. **La función `bh()` del repositorio es correcta**; el error estaba
  solo en el código auxiliar, de modo que el pipeline del artículo no está afectado. Corregido y
  los tres enriquecimientos rehechos.
- **Dos cantidades distintas que conviene no confundir**: la correlación de un gen con la
  *puntuación* de glucocorticoides (FKBP5 = +0.817) y la *concordancia entre tejidos* de ese gen
  (FKBP5 = +0.838 en músculo / +0.762 en adiposo). El documento usa ambas y ahora las distingue.
- **Calibración del nulo verificada**: 300 conjuntos aleatorios de 50 genes dan 5.3 % de P < 0.05,
  el valor esperado. La dispersión de un conjunto real (0.197) es igual a la de uno aleatorio
  (0.198), de modo que la correlación entre genes de una vía no infla estos P como yo temía.
- **Cifras cotejadas contra los TSV almacenados**: modelo jerárquico b1 = −0.716 y b2 = +1.027 con
  208 participantes, ajuste por composición 0.463 sin ajustar y 0.082 ajustado. Coinciden.

## Controles de circularidad: hechos

El conjunto de glucocorticoides procede de la literatura, no de estos datos, de modo que su
enriquecimiento no es circular por construcción. Además:

- **Sin solapamiento** entre el conjunto de glucocorticoides (9 genes) y el núcleo de insulina
  (6 genes): su correlación de 0.76 no viene de compartir miembros.
- **No depende de genes dominantes.** Retirando FKBP5, la puntuación correlaciona 0.995 con la
  completa y su relación con el núcleo de insulina pasa de 0.765 a 0.759. Retirando GLUL, 0.738;
  retirando ZBTB16, 0.762.
- **El programa catabólico se conserva con genes ajenos a las puntuaciones.** Ocho genes de atrofia
  y autofagia, ninguno usado para definir nada, dan r = 0.52 con el eje frente a −0.01 del nulo.
  El nulo **equiparado en nivel de expresión y variabilidad** da el mismo resultado que el nulo por
  tamaño (P = 0.002 en ambos), que era la objeción específica.

Queda pendiente recalcular con decisiones congeladas el eje definido desde los 800 genes más
variables de GTEx, que sí es circular.

## Controles pendientes
2. **Fisiológico frente a agónico.** Presentar tres estimaciones: concordancia total en GTEx,
   concordancia dentro de estratos de circunstancia de muerte y tras ajustes técnicos, y
   concordancia en personas vivas. Distinguir circunstancia de muerte registrada de exposición
   documentada a corticoides: la ventilación no demuestra tratamiento en un donante concreto.
3. **Replicación de la asociación con T2D** fuera de GSE22309, con firma fijada de antemano.
4. **Cuantificar la fracción de covariación explicada** frente a otros componentes, para decidir
   entre «componente prominente» y cualquier afirmación más fuerte.
