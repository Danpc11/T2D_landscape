"""Readers for GEO series matrices, platform annotations and common preprocessing."""
import gzip, re, os
import numpy as np, pandas as pd

RAW = os.environ.get("T2D_RAW", "data/raw")   # carpeta con los archivos descargados de GEO/GTEx

def path(name):
    p = os.path.join(RAW, name)
    if not os.path.exists(p): raise FileNotFoundError(f"{p} no existe: descargar segun data/README.md")
    return p

def read_series_matrix(name):
    """Devuelve (pheno, expr). pheno: gsm, title, y una columna por 'clave: valor' de characteristics.
    expr: tabla de la series matrix (None si no trae datos)."""
    rows, tab, intab = {}, [], False
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path(name), "rt", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"): intab = True; continue
            if line.startswith("!series_matrix_table_end"): intab = False; continue
            if intab: tab.append(line.rstrip("\n").split("\t")); continue
            if line.startswith("!Sample_geo_accession"): gsm = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_title"): title = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
            if line.startswith("!Sample_characteristics_ch1"):
                vals = [x.strip('"') for x in line.rstrip("\n").split("\t")[1:]]
                key = next((v.split(":")[0].strip() for v in vals if ":" in v), None)
                if key: rows[key] = [v.split(":", 1)[1].strip() if ":" in v else "" for v in vals]
    p = pd.DataFrame(rows); p.insert(0, "title", title); p.insert(0, "gsm", gsm)
    e = None
    if len(tab) > 1:
        e = pd.DataFrame(tab[1:], columns=[c.strip('"') for c in tab[0]]).set_index(tab[0][0].strip('"'))
        e.index = [i.strip('"') for i in e.index]; e = e.apply(pd.to_numeric, errors="coerce").dropna()
    return p, e

def read_gpl_annot(name, id_col=0):
    """Mapa probe -> simbolo desde GPLxxx.annot.gz o GPLxxx-NNNNN.txt (HTA)."""
    opener = gzip.open if name.endswith(".gz") else open
    with opener(path(name), "rt", errors="ignore") as fh:
        lines = [l for l in fh if not l.startswith(("#", "!", "^"))]
    hdr = lines[0].rstrip("\n").split("\t"); rows = [l.rstrip("\n").split("\t") for l in lines[1:]]
    df = pd.DataFrame(rows, columns=hdr[:len(rows[0])])
    if "gene_assignment" in df.columns:          # HTA 2.0
        def sym(s):
            if not isinstance(s, str) or s == "---": return None
            parts = [x.strip() for x in s.split("//")]; return parts[1] if len(parts) > 1 else None
        return {str(i): sym(s) for i, s in zip(df.iloc[:, id_col], df["gene_assignment"])}
    sym_col = [c for c in df.columns if "symbol" in c.lower()][0]
    return {str(i): (s if isinstance(s, str) and s and "///" not in s else None) for i, s in zip(df.iloc[:, id_col], df[sym_col])}

def probes_to_genes(e, amap):
    e = e.copy(); e.index = [amap.get(i) for i in e.index]
    e = e[[i is not None for i in e.index]]; return e.groupby(level=0).mean()

def ensembl_map_from_gtex(gct_name="gene_reads_adult_gtex_v11_liver.gct.gz"):
    g = pd.read_csv(path(gct_name), sep="\t", skiprows=2, usecols=["Name", "Description"])
    return dict(zip(g.Name.str.split(".").str[0], g.Description))

def log_cpm(counts, min_logcpm=1.0, frac=0.5):
    x = np.log2(counts / counts.sum(axis=0) * 1e6 + 1)
    return x[(x > min_logcpm).mean(axis=1) > frac]

def expand_characteristics(p):
    from landscape import expand_characteristics as _e
    return _e(p)
