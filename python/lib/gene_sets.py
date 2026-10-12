"""Registro UNICO de los conjuntos de genes del articulo.

Un solo lugar de definicion para evitar que el codigo y las tablas diverjan. Cada conjunto lleva
su fuente. Ningun conjunto se deriva de los datos que se ponen a prueba.

Los scripts y las tablas suplementarias deben importar de aqui; nunca redefinir in situ.
"""

SETS = {
    "Genetic": dict(
        genes=["GSTM1", "GSTT1", "ERAP2"],
        source="Genes silenced by common structural variants: whole-gene deletions of GSTM1 and "
               "GSTT1, and the haplotype silencing ERAP2 (Andres et al., PLoS Genet 2010, "
               "doi:10.1371/journal.pgen.1001157). Pre-specified as the primary genetic set."),
    "Genetic (HLA)": dict(
        genes=["HLA-DQA1", "HLA-DRB4"],
        source="HLA genes whose expression depends on haplotype-dependent gene presence. Reported "
               "separately because they behave differently from the deletion set: adding them to "
               "it raises between-cohort heterogeneity from I2 = 0% to 89%."),
    "Clock output": dict(
        genes=["DBP", "PER1", "PER2", "PER3", "NR1D1", "NR1D2", "CIART", "TEF", "HLF"],
        source="Output arm of the circadian oscillator; PAR-bZIP factors and Rev-erb targets."),
    "Core oscillator": dict(
        genes=["ARNTL", "CLOCK", "CRY1", "CRY2", "NPAS2"],
        source="Core transcriptional-translational feedback loop."),
    "Glucocorticoid": dict(
        genes=["FKBP5", "TSC22D3", "ZBTB16", "KLF15", "DDIT4", "SESN1", "GLUL", "SGK1", "CEBPD"],
        source="Canonical glucocorticoid-receptor target genes."),
    "Insulin core": dict(
        genes=["TXNIP", "PDK4", "HK2", "IRS2", "ANGPTL4", "RRAD"],
        source="Genes with the most consistent response to insulin in paired clamp biopsies."),
    "Oxidative": dict(
        genes=["PPARGC1A", "NDUFA4", "NDUFB8", "SDHB", "UQCRC1", "COX5A", "ATP5F1A", "CYCS",
               "ESRRA", "TFAM"],
        source="PGC-1alpha-associated oxidative phosphorylation programme (Mootha 2003, Patti 2003)."),
}

def genes(name):
    return list(SETS[name]["genes"])

def as_table():
    import pandas as pd
    return pd.DataFrame([dict(gene_set=k, n_genes=len(v["genes"]),
                              genes="; ".join(v["genes"]), source=v["source"])
                         for k, v in SETS.items()])
