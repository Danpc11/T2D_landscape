import re
#!/usr/bin/env python3
"""Inserta citas numeradas en el manuscrito y renumera por orden de aparicion (estilo Nature)."""
import re, pandas as pd
from docx import Document

def sent_end(t, i):
    """Fin de frase a partir de i: un punto que no sea decimal ni de abreviatura."""
    for m in re.finditer(r"\.(?=\s|$)", t[i:]):
        j = i + m.start()
        if j > 0 and t[j-1].isdigit() and j+1 < len(t) and t[j+1:j+2].isdigit(): continue
        return j
    return len(t)

EX = [l.strip() for l in open("existing_refs.txt") if l.strip()]
ex = {}                                   # clave -> cita, a partir de las 19 existentes
keys19 = ["defronzo09","richter25","james21","mootha03","patti03","wu07","coletta08","tuvia21",
          "crosby19","dyar14","wefers18","gabriel21","raju26","rubenstein20","varemo17","pillon22",
          "dollet24","ryden16","yoshino21"]
for k, l in zip(keys19, EX): ex[k] = re.sub(r"^\d+\.\s*", "", l)

V = pd.read_csv("/mnt/user-data/uploads/referencias_verificadas.csv")
def cite(bloque, idx):
    r = V[(V.bloque == bloque) & (V.id == f"{bloque}{idx}")]
    return str(r.cita_verificada.iloc[0]).strip() if len(r) else None
new = {
 # bloque A, fuentes primarias
 "makhnov23": cite("A",1), "gtex20": cite("A",2), "solimena18": cite("A",3), "wigger21": cite("A",4),
 "fadista14": cite("A",5), "marselli20": cite("A",7), "gallagher10": cite("A",8), "jin11": cite("A",9),
 "gse7146": cite("A",10), "dreher23": cite("A",11), "keller11": cite("A",12), "winnier15": cite("A",13),
 "raulerson19": cite("A",14), "hardy11": cite("A",15), "soronen12": cite("A",16), "arner_gse59034": cite("A",17),
 "hanssen15": cite("A",18), "poitou15": cite("A",19), "karolina11": cite("A",21),
 # bloque B, metodos
 "bh95": cite("B",1), "mardia00": cite("B",2), "banerjee05": cite("B",4), "schonemann66": cite("B",5),
 "wgcna08": cite("B",6), "luxburg07": cite("B",7), "limma15": cite("B",8), "scanpy18": cite("B",11),
 "leiden19": cite("B",12), "aitchison82": cite("B",13), "cellcode15": cite("B",14), "tost87": cite("B",15),
 # bloque C, contexto
 "clamp79": cite("C",1), "batista20": cite("C",2), "barrett09": cite("C",3), "clerk06": cite("C",4),
 "williams20": cite("C",5), "stenvers19": cite("C",6), "dyar18": cite("C",7), "scheer09": cite("C",9),
 "sylow17": cite("C",10),
}
new = {k: v for k, v in new.items() if v}
ALL = dict(ex); ALL.update(new)

# anclas: fragmento unico del parrafo -> claves a citar en ese punto
ANCH = [
 ("contributes to type 2 diabetes", ["defronzo09","richter25"]),
 ("lower expression of oxidative metabolism gene", ["mootha03","patti03","james21"]),
 ("can influence peripheral circadian clocks", ["wu07","coletta08","tuvia21","crosby19"]),
 ("the intrinsic muscle clock regulates aspects of g", ["dyar14","wefers18","gabriel21"]),
 ("point in different directions", ["raju26"]),
 ("their relative contributions can affect the measured transcriptome", ["rubenstein20"]),
 ("Cultured m", ["varemo17"]),
 ("spectral embeddings of their normalized graph Laplacians", ["luxburg07","wgcna08"]),
 ("Among 253 GTEx donors", ["gtex20"]),
 ("the islet datasets yielded 1,284, 45 and zero genes", ["wigger21","fadista14","marselli20"]),
 ("the METSIM adipose com", ["raulerson19"]),
 ("the independent muscle dataset yielded zero", ["jin11"]),
 ("49 participants studied in two laboratories", ["gabriel21"]),
 ("hyperinsulinemic–eugly", ["clamp79"]),
 ("Benjamini–Hochberg correction", ["bh95"]),
 ("0.20 at 30 min, 0.24 at 1 h after a mixed meal", ["coletta08","makhnov23"]),
 ("an independent prediabetic cohort", ["yoshino21"]),
 ("3 h into recovery from acute exercise", ["pillon22"]),
 ("An adip", ["dollet24"]),
 ("Other perturbations showed a broader range of coherence", ["dreher23","hanssen15"]),
 ("Two clamp comparisons retained similar mean response directions", ["ryden16","soronen12"]),
 ("published fiber-type marker panels and a single-cell reference", ["rubenstein20"]),
 ("Primary myotubes from 24 donors", ["varemo17"]),
 ("whereas adjustment for a fibe", ["cellcode15"]),
 ("Age, adiposity and study design influence both expression", ["batista20"]),
 ("Coherence describes agreement within a group", ["raju26","mardia00"]),
 ("altered coupling between insulin stimulation and clock output", ["dyar18","stenvers19","scheer09"]),
 ("another stimulus", ["sylow17"]),
 ("vascular, stromal and immune cells", ["barrett09","clerk06","williams20"]),
 # metodos
 ("absolute biweight midcorrelation raised to a soft-thresholding power", ["wgcna08"]),
 ("compared by orthogonal", ["schonemann66"]),
 ("von Mises–Fisher distribution by maximum likelihood", ["mardia00","banerjee05"]),
 ("The normalizing constant was evaluated", ["banerjee05"]),
 ("Table 1 of Rubenstein et al", ["rubenstein20"]),
 ("signature matrix derived from GSE130646", ["rubenstein20","scanpy18","leiden19"]),
 ("Before centered log-ratio transformation", ["aitchison82"]),
 ("two-one-sided tests", ["tost87"]),
 ("Welch t-tests", ["limma15"]),
 # fuentes de datos que solo aparecen en el inventario: se citan al presentarlo
 ("Supplementary Table 1 l", ["solimena18","gallagher10","gse7146","keller11","winnier15","hardy11",
                              "arner_gse59034","poitou15","karolina11"]),
]
doc = Document("/mnt/user-data/uploads/MANUSCRIPT_revised.docx")
order, assigned = [], 0
def num(k):
    if k not in order: order.append(k)
    return order.index(k) + 1
SUP = str.maketrans("0123456789,", "⁰¹²³⁴⁵⁶⁷⁸⁹𝄒")
def supfmt(ns):
    out = ",".join(str(n) for n in ns)
    return out.translate(str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")).replace(",", "ʼ")
paras = [p for p in doc.paragraphs]
stop = next((i for i, p in enumerate(paras) if p.text.strip() == "References"), len(paras))
for p in paras[:stop]:
    t = p.text
    if not t.strip(): continue
    # agrupar las anclas que caen en la misma frase, para no fusionar superindices
    hits = {}
    for anch, keys in ANCH:
        if anch not in t: continue
        i = t.find(anch) + len(anch); j = sent_end(t, i)
        hits.setdefault(j, []).extend(k for k in keys if k in ALL)
    for j in sorted(hits):            # primero asignar numeros en orden de lectura
        for k in hits[j]:
            if k in ALL: num(k)
    for j in sorted(hits, reverse=True):   # luego insertar de atras hacia delante
        seen, ks = set(), []
        for k in hits[j]:
            if k not in seen: seen.add(k); ks.append(k)
        ns = sorted(num(k) for k in ks)
        t = t[:j] + supfmt(ns) + t[j:]
        assigned += 1
    if t != p.text:
        for r in p.runs[1:]: r.text = ""
        p.runs[0].text = t
# reemplazar la seccion de referencias
for k in ALL:                      # las no citadas van al final, en bloque
    num(k)
ref_par = paras[stop]
idx = paras.index(ref_par)
# borrar las referencias antiguas
for p in paras[idx+1:]:
    if p.text.strip().startswith(("Figure legends", "Fig. 1")): break
    p.text = ""
body = ref_par._element.getparent()
anchor = ref_par._element
from copy import deepcopy
for n, k in enumerate(order, 1):
    new_p = deepcopy(ref_par._element); body.insert(list(body).index(anchor) + n, new_p)
    from docx.text.paragraph import Paragraph
    Paragraph(new_p, ref_par._parent).text = f"{n}. {ALL[k]}"
doc.save("/mnt/user-data/outputs/MANUSCRIPT_with_refs.docx")
import json
json.dump({"order": order, "refs": ALL}, open("ref_order.json", "w"), ensure_ascii=False, indent=1)
print(f"anclas aplicadas: {assigned} | referencias totales: {len(order)}")
print("primeras 12 por orden de aparicion:")
for n, k in enumerate(order[:12], 1): print(f"  {n}. {ALL[k][:88]}")
