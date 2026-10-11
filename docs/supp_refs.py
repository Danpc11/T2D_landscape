import re
#!/usr/bin/env python3
"""Corrige el inventario de la Informacion Suplementaria y le inserta las citas numeradas,
usando la misma numeracion que el manuscrito (ref_order.json)."""
import json, re
from docx import Document

def sent_end(t, i):
    """Fin de frase a partir de i: un punto que no sea decimal ni de abreviatura."""
    for m in re.finditer(r"\.(?=\s|$)", t[i:]):
        j = i + m.start()
        if j > 0 and t[j-1].isdigit() and j+1 < len(t) and t[j+1:j+2].isdigit(): continue
        return j
    return len(t)
D = json.load(open("ref_order.json")); ORDER = D["order"]; REFS = D["refs"]
N = {k: i + 1 for i, k in enumerate(ORDER)}
def sup(ns):
    s = "ʼ".join("".join("⁰¹²³⁴⁵⁶⁷⁸⁹"[int(c)] for c in str(n)) for n in sorted(ns))
    return s
# accession -> clave de referencia
ACC = {"GSE76895":"solimena18","GSE164416":"wigger21","GSE50244":"fadista14","GSE50398":"fadista14",
 "GSE159984":"marselli20","GSE18732":"gallagher10","GSE25462":"jin11","GSE182120":"gabriel21",
 "GSE182117":"gabriel21","GSE81965":"varemo17","GSE63887":"varemo17","GSE130646":"rubenstein20",
 "GSE22309":"wu07","GSE9105":"coletta08","GSE7146":"gse7146","GSE231509":"makhnov23",
 "GSE157988":"yoshino21","GSE202295":"pillon22","GSE224310":"dreher23","GSE106800":"wefers18",
 "GSE27951":"keller11","GSE64567":"winnier15","GSE135134":"raulerson19","GSE20950":"hardy11",
 "GSE26637":"soronen12","GSE59034":"arner_gse59034","Ryden2016":"ryden16","GSE67297":"hanssen15",
 "GSE66306":"poitou15","GSE156993":None,"GSE21321":"karolina11","GSE198922":"dollet24","GTEx":"gtex20"}
# correcciones de atribucion en los encabezados del inventario
FIX = {"GSE7146 Rome 2007":"GSE7146 Parikh, Groop and Mootha (depositors)",
 "GSE27951 2010":"GSE27951 Keller 2011","GSE64567 2015":"GSE64567 Winnier 2015",
 "GSE135134 METSIM":"GSE135134 Raulerson 2019 (METSIM)","GSE59034 Mitchell 2015":"GSE59034 Arner, Rydén and Mejhert (depositors)",
 "GSE66306 2015":"GSE66306 Poitou 2015","GSE21321 2010":"GSE21321 Karolina 2011",
 "GSE156993 2020":"GSE156993 no associated publication","GSE231509 2023":"GSE231509 Makhnovskii 2023",
 "GSE224310 2023":"GSE224310 Dreher 2023"}
BODY = {"cross-sectional (overlapping donors with GSE50244, array; two platforms in one series)":
        "SuperSeries grouping GSE50244 (RNA-seq) and GSE50397 (expression array); the array accession is GSE50397"}
# anclas en las notas suplementarias
ANCH = [("allowing for nonunique eigenbases", ["schonemann66","luxburg07"]),
        ("The single-cell CRISPR preprint reports correlations", ["raju26"]),
        ("alignment with a chosen reference describe different aspects", ["mardia00","banerjee05"]),
        ("the reference supplies 400 signature genes", ["rubenstein20","scanpy18","leiden19"]),
        ("relative signature contributions", ["aitchison82","cellcode15"]),
        ("ex vivo islet datasets examined donor preparations", ["marselli20"]),
        ("diurnal sampling and cold acclimation", ["wefers18","hanssen15","dreher23"])]
doc = Document("/mnt/user-data/uploads/SUPPLEMENTARY_T2D.docx")
nfix = nanch = nacc = 0
for p in doc.paragraphs:
    t = p.text
    if not t.strip(): continue
    orig = t
    for a, b in FIX.items():
        if t.startswith(a): t = t.replace(a, b, 1); nfix += 1
    for a, b in BODY.items():
        if a in t: t = t.replace(a, b); nfix += 1
    # cita tras el encabezado de cada accession del inventario
    m = re.match(r"^(GSE\d+|Ryden2016|GTEx[^\s]*)\s", t)
    if m and len(t) < 110:
        k = ACC.get(m.group(1))
        if k and k in N and sup([N[k]]) not in t:
            t = t.rstrip() + sup([N[k]]); nacc += 1
    hits = {}
    for anch, keys in ANCH:
        if anch in t:
            i = t.find(anch) + len(anch); j = sent_end(t, i)
            hits.setdefault(j, []).extend(k for k in keys if k in N)
    for j in sorted(hits, reverse=True):
        seen, ks = set(), []
        for k in hits[j]:
            if k not in seen: seen.add(k); ks.append(k)
        t = t[:j] + sup([N[k] for k in ks]) + t[j:]; nanch += 1
    if t != orig:
        for r in p.runs[1:]: r.text = ""
        if p.runs: p.runs[0].text = t
        else: p.text = t
# nota final con la correspondencia de numeracion
p = doc.add_paragraph()
p.text = ("Reference numbers follow the main text. Accessions without an associated publication are "
          "cited by accession only. GSE50398 is the SuperSeries grouping GSE50244 and GSE50397; "
          "GSE63887 supplies the baseline samples of GSE81965 and adds no donors.")
doc.save("/mnt/user-data/outputs/SUPPLEMENTARY_with_refs.docx")
print(f"correcciones de atribucion: {nfix} | citas en el inventario: {nacc} | citas en las notas: {nanch}")
