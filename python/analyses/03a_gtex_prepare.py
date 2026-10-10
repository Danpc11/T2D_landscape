#!/usr/bin/env python3
"""Figura 5 / R1: GTEx v11 tejidos emparejados. 03a prepara matrices por donante; 03b: coherencia del haz, coherencia individual (nulo de barajado de donantes), sangre -> tejido."""
import pandas as pd, numpy as np, warnings, os; warnings.filterwarnings("ignore")
import os; U=os.environ.get("T2D_RAW","data/raw"); OUT=os.environ.get("T2D_OUT","results/gtex"); os.makedirs(OUT, exist_ok=True); os.makedirs("gtex", exist_ok=True)
s=pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt",sep="\t",low_memory=False)
s["donor"]=s.SAMPID.str.split("-").str[:2].str.join("-")
ph=pd.read_csv(f"{U}/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt",sep="\t").set_index("SUBJID")
# once tejidos: cinco metabolicos, sangre, y cinco no metabolicos que sirven de control de
# especificidad (la arquitectura compartida, es propia de los organos del metabolismo o de
# cualquier par de tejidos de la misma persona?)
tissues={"muscle":"muscle_skeletal","adipose":"adipose_subcutaneous","adipose_visceral":"adipose_visceral_omentum",
         "pancreas":"pancreas","liver":"liver","blood":"whole_blood",
         "adrenal":"adrenal_gland","kidney_cortex":"kidney_cortex","kidney_medulla":"kidney_medulla",
         "stomach":"stomach","ileum":"small_intestine_terminal_ileum"}
def load(t):
    g=pd.read_csv(f"{U}/gene_reads_adult_gtex_v11_{t}.gct.gz",sep="\t",skiprows=2)
    g=g[~g.Description.str.startswith(("MT-","RPL","RPS"))]
    g=g.set_index("Description").drop(columns=["Name","id"],errors="ignore"); g=g.groupby(level=0).sum()
    g=g.loc[:, g.columns.isin(s.SAMPID)]
    cpm=np.log2(g/g.sum(axis=0)*1e6+1); cpm=cpm[(cpm>1).mean(axis=1)>0.5]
    cpm.columns=[c.split("-")[0]+"-"+c.split("-")[1] for c in cpm.columns]
    cpm=cpm.loc[:, ~cpm.columns.duplicated()]
    return cpm
os.makedirs("gtex",exist_ok=True)
for k,t in tissues.items():
    e=load(t); e.to_pickle(f"{OUT}/{k}.pkl"); print(k, e.shape)
