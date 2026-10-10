# Data: what to download and where to put it

All inputs are public. Place them in `data/raw/` with exactly these names (GEO file names with
`_txt.gz` renamed to `.txt.gz` where shown). Nothing in `data/raw/` is tracked by git.

## Discovery (downloaded automatically by `R_scripts_preprocessing/01_*.R`)
GSE76895, GSE18732, GSE15653, GSE27951 (GEOquery; needs internet from the R session).

## Replication of the landscape (python/analyses/01_prepare_replication_cohorts.py)
| File | Source |
|---|---|
| GSE164416_series_matrix.txt.gz, GSE164416_DP_htseq_counts.txt.gz | GEO GSE164416 |
| GSE50244_series_matrix.txt.gz, GSE50244_Genes_counts_TMM_NormLength_atLeastMAF5_expressed.txt.gz | GEO GSE50244 |
| GSE50398-GPL6244_series_matrix.txt.gz | GEO GSE50398 |
| GSE25462_series_matrix.txt.gz | GEO GSE25462 |
| GSE64567_series_matrix.txt.gz | GEO GSE64567 |
| GSE59034_series_matrix.txt.gz | GEO GSE59034 |
| GSE135134_series_matrix.txt.gz, GSE135134_METSIM_subcutaneousAdipose_RNAseq_TPMs_n434.txt | GEO GSE135134 (supplementary file) |

## GTEx paired tissues (python/analyses/03*)
gene_reads_adult_gtex_v11_{muscle_skeletal,adipose_subcutaneous,pancreas,whole_blood,liver}.gct.gz (GTEx portal, open),
GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt, GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt (open).
The liver gct is also used as the Ensembl→symbol map.

## Insulin response (python/analyses/04_insulin_response.py)
| File | Source |
|---|---|
| GSE22309_series_matrix.txt.gz, GPL91.annot.gz | GEO GSE22309; platform GPL91 "Download full table" |
| GSE9105_series_matrix.txt.gz, GPL96.annot.gz | GEO GSE9105; platform GPL96 |
| GSE7146-GPL96_series_matrix.txt.gz | GEO GSE7146 |
| GSE231509_Read_counts.txt.gz | GEO GSE231509 |
| GSE157988_series_matrix.txt.gz, GSE157988_NMN_all_gene_counts.xlsx | GEO GSE157988 |
| GSE26637_series_matrix.txt.gz | GEO GSE26637 |
| ExpressionTable.txt, Cohort.txt | Rydén et al. 2016 Cell Reports export (export.uppmax.uu.se/b2013047/CellReportsTables/; mirror in this repo's release assets if the server is down) |

## Resting muscle vs insulin sensitivity (05) and myotubes (06)
GSE182120_series_matrix.txt.gz, GPL17586-45144.txt (HTA 2.0 annotation, GEO platform page); GSE182117_counts.tsv.gz.

## Supplementary (07)
GSE66306_PM_processed_counts.txt.gz, GSE156993_series_matrix.txt.gz, GSE21321-GPL6883_series_matrix.txt.gz, GSE129843_RESTRICT.txt.gz.

Set `T2D_RAW` to point elsewhere if needed (default `data/raw`).

## Inputs added for the response, specificity and model analyses

All from GEO unless stated. For RNA-seq series, use the NCBI-generated count matrix
(`GSE*_raw_counts_GRCh38_p13_NCBI.tsv.gz`) together with the shared annotation file, because several
series deposit no expression table in the series matrix.

| File | Used by | Analysis |
|---|---|---|
| `Human_GRCh38_p13_annot.tsv.gz` | NCBI, shared by all count matrices | gene symbol mapping |
| `GSE202295_raw_counts_GRCh38_p13_NCBI.tsv.gz` + series matrix | 11, 16, 18 | acute exercise in muscle, NGT vs T2D (specificity) |
| `GSE198922_raw_counts_GRCh38_p13_NCBI.tsv.gz` + series matrix | 11, 18 | the same protocol sampled in adipose tissue |
| `GSE224310_series_matrix.txt.gz` | Supplementary Note 3 | acute exercise versus weeks of training in the same participants |
| `GSE159984_raw_counts_GRCh38_p13_NCBI.tsv.gz` + both series matrices (GPL9115, GPL16791) | 09, ED Fig. 2 | islet donors at rest, and ex vivo perturbation with paired controls |
| `GSE106800_series_matrix.txt.gz` | Supplementary Note 3 | circadian misalignment, within person |
| `GSE67297_series_matrix.txt.gz` | Supplementary Note 3 | ten days of cold acclimation in T2D |
| `GSE182117_raw_counts_GRCh38_p13_NCBI.tsv.gz` + series matrix | 06 | myotube time series; preferred over the authors' counts for consistency |
| `gene_reads_adult_gtex_v11_{adipose_visceral_omentum,adrenal_gland,kidney_cortex,kidney_medulla,stomach,small_intestine_terminal_ileum}.gct.gz` | 03a, 14 | the non-metabolic tissues used as the specificity control for the shared individual position |

Note on a dataset that looks usable and is not: GSE182121 is the superseries of GSE182120 and
GSE182117. Its microarray samples are numbered oddly (001, 003, 007 …) and all are labelled Basal.
The clamp in that study was used to phenotype participants, not as a paired perturbation, so no
insulin-stimulated biopsies exist to request.
