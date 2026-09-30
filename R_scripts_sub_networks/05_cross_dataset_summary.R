#!/usr/bin/env Rscript
# =============================================================================
# 05_cross_dataset_summary.R
# Consolidate metrics and gene-driver scores into cross-dataset summary tables.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# BRANCH: Shared
#
# CORRECTIONS:
#   [C18] trend_summary: physiological ordering via phys_order before
#         first()/last(). Prevents NOI_change inversion due to alphabetical order.
#   [C19] Datasets and state_orders updated to the microarray quartet.
#   [NEW] Warn if the join generates NAs in phys_order (unrecognized conditions).
#   [FIX-regex] patrones anclados ^GSE: antes capturaban all_*.tsv -> filas duplicadas.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# [CLI] Banderas de linea de comandos: --nombre=valor (tienen prioridad sobre la
# variable de entorno homonima, que se conserva por compatibilidad).
#   Rscript script.R [--workers=32] [--n_boot=100] [...] [GSE... GSE...]
# Los argumentos sin "--" son los accessions a procesar.
# -----------------------------------------------------------------------------
.cli_raw   <- commandArgs(trailingOnly = TRUE)
.cli_flags <- grep("^--", .cli_raw, value = TRUE)
cli_args   <- grep("^--", .cli_raw, value = TRUE, invert = TRUE)
get_flag <- function(name, env = NULL, default) {
  hit <- grep(paste0("^--", name, "="), .cli_flags, value = TRUE)
  if (length(hit) > 0) return(sub(paste0("^--", name, "="), "", hit[1]))
  if (!is.null(env)) { v <- Sys.getenv(env, unset = NA); if (!is.na(v) && nzchar(v)) return(v) }
  as.character(default)
}
if (any(.cli_flags == "--help")) {
  cat("Uso: Rscript", basename(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE))),
      "[--flag=valor ...] [GSE...]\nBanderas: ver cabecera del script (get_flag).\n"); quit(status = 0)
}

options(stringsAsFactors = FALSE)

dir.create("results/summary", recursive = TRUE, showWarnings = FALSE)

targets <- if (length(cli_args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else cli_args

# -----------------------------------------------------------------------------
# MAPEO BIOLÓGICO 🔬
# -----------------------------------------------------------------------------
dataset_info <- data.frame(
  accession = c("GSE76895", "GSE18732", "GSE15653", "GSE27951"),
  organ     = c(
    "Pancreatic islets",   # GSE76895: Taneera et al. 2012
    "Skeletal muscle",     # GSE18732: Gallagher et al. 2010
    "Liver",               # GSE15653: Pihlajamaki et al. 2009 (verificar en GEO)
    "Adipose tissue"       # GSE27951: Keller et al. 2011 (verificar en GEO)
  ),
  stringsAsFactors = FALSE
)

state_orders <- list(
  GSE76895 = c("ND",  "IGT", "T2D"),
  GSE18732 = c("ND",  "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

message("=== Consolidación con contexto biológico ===")

# -----------------------------------------------------------------------------
# ORDEN FISIOLÓGICO
# -----------------------------------------------------------------------------
order_df <- bind_rows(lapply(targets, function(acc) {
  ord <- state_orders[[acc]]
  data.frame(
    accession  = acc,
    condition  = ord,
    phys_order = seq_along(ord),
    stringsAsFactors = FALSE
  )
}))

# -----------------------------------------------------------------------------
# LECTOR GENÉRICO
# -----------------------------------------------------------------------------
read_and_label <- function(dir_path, pattern, branch_name) {
  files <- list.files(dir_path, pattern = pattern, full.names = TRUE)
  if (length(files) == 0L) return(NULL)

  bind_rows(lapply(files, function(f) {
    d <- fread(f)
    acc_name <- sub("_.*", "", basename(f))
    if (!"accession" %in% colnames(d)) d$accession <- acc_name
    d$branch <- branch_name
    d
  }))
}

# -----------------------------------------------------------------------------
# 1. MÉTRICAS
# -----------------------------------------------------------------------------
message("1. Métricas...")

metrics <- bind_rows(
  read_and_label("results/metrics", "^GSE[0-9]+_network_metrics\\.tsv$", "shared"),
  read_and_label("results/full_metrics", "^GSE[0-9]+_full_network_metrics\\.tsv$", "full")
)

if (nrow(metrics) > 0) {

  metrics_ordered <- metrics |>
    left_join(order_df, by = c("accession", "condition")) |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, phys_order)

  # [C18] Verificar que el join no produjo NAs en phys_order
  # (indicaria condiciones no reconocidas en state_orders)
  na_rows <- metrics_ordered[is.na(metrics_ordered$phys_order), ]
  if (nrow(na_rows) > 0L) {
    warning(nrow(na_rows), " filas con phys_order=NA (condicion no reconocida): ",
            paste(unique(paste0(na_rows$accession, "/", na_rows$condition)),
                  collapse = ", "),
            "\nEstas filas se EXCLUYEN de trend_summary para evitar first()/last() incorrectos.")
    metrics_ordered <- metrics_ordered[!is.na(metrics_ordered$phys_order), ]
  }

  trend_summary <- metrics_ordered |>
    group_by(branch, accession, organ) |>
    summarise(
      healthiest  = first(condition),
      final_state = last(condition),
      EG_change   = last(EG)   - first(EG),
      Hb_change   = last(Hb)   - first(Hb),
      Rbar_change = last(Rbar) - first(Rbar),
      Gbar_change = last(Gbar) - first(Gbar),
      # [NEW] cambios a n igual y relativos al nulo (rama shared; NA en full)
      EG_sub_change  = if ("EG_sub" %in% names(pick(everything()))) last(EG_sub) - first(EG_sub) else NA_real_,
      EG_rel_change  = if ("EG_rel" %in% names(pick(everything()))) last(EG_rel) - first(EG_rel) else NA_real_,
      NOI_change     = if ("NOI" %in% names(pick(everything()))) last(NOI) - first(NOI) else NA_real_,
      .groups = "drop"
    )

  fwrite(metrics_ordered |> select(-phys_order),
         "results/summary/all_metrics_combined.tsv", sep = "\t")

  fwrite(trend_summary,
         "results/summary/cross_dataset_trends.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 2. DRIVERS
# -----------------------------------------------------------------------------
message("2. Drivers...")

driver_df <- bind_rows(
  read_and_label("results/gene_drivers", "_gene_driver_scores\\.tsv$", "shared"),
  read_and_label("results/full_gene_drivers", "_full_gene_driver_scores\\.tsv$", "full")
)

if (nrow(driver_df) > 0) {

  driver_df <- driver_df |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")"))

  fwrite(driver_df,
         "results/summary/all_gene_driver_scores_combined.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 3. OPTIMUM
# -----------------------------------------------------------------------------
message("3. Configuration reference...")

opt_df <- bind_rows(
  read_and_label("results/reference", "_deviation_from_reference\\.tsv$", "shared"),
  # [CORRECCIÓN] Carpeta cambiada a "results/full_reference"
  read_and_label("results/full_reference", "_full_deviation_from_reference\\.tsv$", "full") 
)

if (nrow(opt_df) > 0) {

  opt_df <- opt_df |>
    left_join(order_df, by = c("accession", "condition")) |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, phys_order) |>
    select(-phys_order)

  fwrite(opt_df,
         "results/summary/all_deviation_from_reference.tsv", sep = "\t")
}

# -----------------------------------------------------------------------------
# 4. DE INTEGRADO
# -----------------------------------------------------------------------------
message("4. DE integrado...")

de_int_df <- bind_rows(
  # [CORRECCIÓN] Regex anclado al inicio para evitar que cargue los archivos "full"
  read_and_label("results/differential_expression", "^[A-Z0-9]+_integrated_drivers_DE\\.tsv$", "shared"),
  read_and_label("results/differential_expression", "_full_integrated_drivers_DE\\.tsv$", "full")
)

if (nrow(de_int_df) > 0) {

  de_int_df <- de_int_df |>
    left_join(dataset_info, by = "accession") |>
    mutate(dataset_label = paste0(accession, " (", organ, ")")) |>
    arrange(branch, organ, accession, desc(combined_score))

  fwrite(de_int_df,
         "results/summary/all_integrated_drivers_DE.tsv", sep = "\t")
}

message("=== DONE ===")
