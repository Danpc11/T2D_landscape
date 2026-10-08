#!/usr/bin/env Rscript
# =============================================================================
# 04b_full_gene_drivers_and_enrichment.R
# Análisis a nivel de gen (drivers) para redes COMPLETAS (Rama Full).
#
# [FIX-TRI-exact] TRI/KO por recalculo exacto de Rbar (ver funciones).
# [FIX-enrich] GO/KEGG independientes.
# NOTA: 02b guarda ahora W como subred top-3000 por estado; reduce_networks_consistent
#       toma la interseccion entre estados (tipicamente 2000-2800 genes).
#
# LÓGICA DE PARALELIZACIÓN (Estilo 02b):
#   - Se usa mclapply (Forking) en lugar de PSOCK/foreach.
#   - Copy-On-Write (COW): Los 40 workers leen las matrices masivas (W_dis, V)
#     de la memoria compartida sin duplicarlas, ahorrando decenas de GB de RAM.
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(parallel)
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
# [FIX-threads] los forks de mclapply heredan BLAS; forzar 1 hilo por proceso
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)

# -----------------------------------------------------------------------------
# HPC / Paralelismo
# -----------------------------------------------------------------------------
n_workers_total <- as.integer(get_flag("workers", "N_WORKERS", 48))

# Evitar oversubscription BLAS
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

message("Motor HPC: dataset-level mclapply | workers totales = ", n_workers_total, " | BLAS=1")

# -----------------------------------------------------------------------------
# Directorios
# -----------------------------------------------------------------------------
dir.create("results/full_gene_drivers", recursive = TRUE, showWarnings = FALSE)
dir.create("results/full_enrichment",   recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Argumentos
# -----------------------------------------------------------------------------
targets <- if (length(cli_args) == 0) {
  c("GSE76895", "GSE18732", "GSE15653", "GSE27951")
} else {
  cli_args
}

state_orders <- list(
  GSE76895 = c("ND", "IGT", "T2D"),
  GSE18732 = c("ND", "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

# -----------------------------------------------------------------------------
# Parámetros de reducción
# -----------------------------------------------------------------------------
HAS_ENRICH <- requireNamespace("clusterProfiler", quietly = TRUE) &&
              requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (HAS_ENRICH) suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
if (!HAS_ENRICH) message("clusterProfiler/org.Hs.eg.db no disponibles: se omite el enriquecimiento")
MAX_GENES_NODE <- as.integer(get_flag("max_genes_node", "MAX_GENES_NODE", 3000))
MAX_GENES_TRI  <- as.integer(get_flag("max_genes_tri", "MAX_GENES_TRI", 1500))
MIN_COR        <- as.numeric(Sys.getenv("MIN_COR", unset = "0.2"))   # legacy, no usado en node_metrics
N_TOP_ENRICH   <- as.integer(get_flag("n_top_enrich", "N_TOP_ENRICH", 200))
# N_TOP_ENRICH=200: biologicamente interpretable y rapido.
# enrichGO sobre 3000 genes produce p-valores inflados (todo significativo).
# El universo correcto sigue siendo todos los genes de la red reducida.

# -----------------------------------------------------------------------------
# Utilidades
# -----------------------------------------------------------------------------
zscore <- function(x) {
  s <- scale(x)
  s[is.na(s)] <- 0
  as.numeric(s)
}

# nth_largest: umbral proporcional compatible con todas las versiones de R
nth_largest <- function(x, n) {
  n <- min(n, length(x))
  if (n <= 0L) return(max(x, na.rm = TRUE))
  as.numeric(quantile(x, probs = 1 - n / length(x), type = 1, na.rm = TRUE))
}

hb_from_W <- function(W) {
  k <- rowSums(W)
  prob <- k / sum(k)
  -sum(prob * log(prob + 1e-12))
}

# -----------------------------------------------------------------------------
# Selección consistente de genes entre estados
# Usa conectividad promedio entre redes para mantener comparabilidad
# -----------------------------------------------------------------------------
reduce_networks_consistent <- function(nets, max_genes = 3000L) {
  common_genes <- Reduce(intersect, lapply(nets, rownames))
  if (length(common_genes) == 0L) return(NULL)

  avg_strength <- Reduce("+", lapply(nets, function(W) {
    rowSums(W[common_genes, common_genes, drop = FALSE])
  })) / length(nets)

  k_use <- min(max_genes, length(avg_strength))
  top_genes <- names(sort(avg_strength, decreasing = TRUE))[seq_len(k_use)]

  lapply(nets, function(W) W[top_genes, top_genes, drop = FALSE])
}

# Métricas de nodo rápidas
# - strength: sobre la red completa reducida (O(p^2), exacto)
# - eigen centrality y participation: sobre red esparsificada top-1% aristas
#   El umbral min_cor^beta es INCORRECTO para beta alto (0.2^20 ~ 1e-14,
#   elimina casi nada). Usamos umbral proporcional identico a 02b.
# SIN betweenness (cuello de botella O(V*E))
node_metrics_fast <- function(W, beta = NULL, min_cor = NULL) {
  p      <- nrow(W)
  k_full <- rowSums(W)

  # Esparsificacion por top-1% de aristas (escala con p, independiente de beta)
  total_ar <- as.numeric(p) * (p - 1L) / 2
  n_keep   <- min(ceiling(total_ar * 0.01), 500000L)
  w_upper  <- W[upper.tri(W)]
  thr      <- nth_largest(w_upper, n_keep)
  rm(w_upper)

  W_sp <- W
  W_sp[W_sp < thr] <- 0
  diag(W_sp) <- 0

  g_sp <- graph_from_adjacency_matrix(
    W_sp, mode = "undirected", weighted = TRUE, diag = FALSE
  )

  ev <- tryCatch(
    eigen_centrality(g_sp, weights = E(g_sp)$weight)$vector,
    error = function(e) rep(0, p)
  )
  names(ev) <- rownames(W)

  mem <- tryCatch({
    cl <- cluster_fast_greedy(g_sp)
    membership(cl)
  }, error = function(e) {
    setNames(seq_len(p), rownames(W))
  })

  Pc <- sapply(seq_along(k_full), function(i) {
    ki <- k_full[i]
    if (ki == 0) return(0)
    parts <- tapply(W[i, ], mem, sum)
    1 - sum((parts / ki)^2, na.rm = TRUE)
  })

  data.frame(
    gene          = names(k_full),
    strength      = k_full,
    eigen         = unname(ev[names(k_full)]),
    participation = Pc
  )
}

# -----------------------------------------------------------------------------
# [FIX-TRI-exact] TRI y KO por RECALCULO EXACTO de Rbar (y Hb).
# La actualizacion "rank-1" anterior no era valida (dL = diag(dk) - dW no es
# rango 1 y la formula no era la expansion de primer orden; error ~30x en
# simulacion). Para p = 1500, eigen(L, only.values=TRUE) ~0.3-0.5 s/gen, en
# paralelo con mclapply (COW) es viable.
# Definiciones (rama full; en la rama shared con p<=800 se usa EG exacto):
#   TRI_i = [Rbar(W_dis) - Rbar(W_dis rescatada en i)] / Rbar(W_dis) + [Hb_resc - Hb_dis]
#   KO_i  = [Rbar(W_dis KO i) - Rbar(W_dis)] / Rbar(W_dis)   (>0: el gen sostiene la red)
# -----------------------------------------------------------------------------
rbar_exact <- function(W) {
  p   <- nrow(W)
  lam <- eigen(diag(rowSums(W), p) - W, symmetric = TRUE, only.values = TRUE)$values
  tol <- max(abs(lam)) * p * 1e-10
  2 * sum(ifelse(lam > tol, 1 / lam, 0)) / (p - 1)
}

reduce_common_tri <- function(W_dis, common) {
  if (length(common) > MAX_GENES_TRI) {
    k_tri  <- rowSums(W_dis[common, common, drop = FALSE])
    common <- names(sort(k_tri, decreasing = TRUE))[seq_len(MAX_GENES_TRI)]
  }
  common
}

compute_TRI_exact_hpc <- function(W_dis, W_ref, common, n_workers_inner = 1L) {
  common <- reduce_common_tri(W_dis, common)
  W_dis  <- W_dis[common, common, drop = FALSE]
  W_ref  <- W_ref[common, common, drop = FALSE]
  R0 <- rbar_exact(W_dis); H0 <- hb_from_W(W_dis)
  vals <- mclapply(seq_along(common), function(i) {
    W <- W_dis; W[i, ] <- W_ref[i, ]; W[, i] <- W_ref[, i]; W[i, i] <- 0
    (R0 - rbar_exact(W)) / R0 + (hb_from_W(W) - H0)
  }, mc.cores = n_workers_inner)
  setNames(unlist(vals), common)
}

compute_KO_exact_hpc <- function(W_dis, common, n_workers_inner = 1L) {
  common <- reduce_common_tri(W_dis, common)
  W_dis  <- W_dis[common, common, drop = FALSE]
  R0 <- rbar_exact(W_dis)
  vals <- mclapply(seq_along(common), function(i) {
    W <- W_dis; W[i, ] <- 0; W[, i] <- 0
    # el nodo aislado se excluye del calculo (su componente tendria R infinita)
    (rbar_exact(W[-i, -i]) - R0) / R0
  }, mc.cores = n_workers_inner)
  setNames(unlist(vals), common)
}

# -----------------------------------------------------------------------------
# Enriquecimiento
# -----------------------------------------------------------------------------
safe_enrich_robust <- function(genes, universe_genes, prefix) {
  if (!HAS_ENRICH) return(invisible(NULL))
  tryCatch({
    suppressWarnings({
      eg <- suppressMessages(
        bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )
      univ <- suppressMessages(
        bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
      )
    })

    if (is.null(eg) || nrow(eg) < 10L) return(FALSE)
    if (is.null(univ) || nrow(univ) < 10L) return(FALSE)

    ego <- tryCatch(suppressMessages(
      enrichGO(
        gene = eg$ENTREZID,
        universe = univ$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        readable = TRUE
      )
    ), error = function(e) { message("  GO error: ", conditionMessage(e)); NULL })

    # [FIX-enrich] KEGG en su propio tryCatch (requiere internet)
    ekk <- tryCatch(suppressMessages(
      enrichKEGG(
        gene = eg$ENTREZID,
        universe = univ$ENTREZID,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05
      )
    ), error = function(e) { message("  KEGG error: ", conditionMessage(e)); NULL })

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      fwrite(as.data.frame(ego),
             file.path("results/full_enrichment", paste0(prefix, "_GO_3k.tsv")),
             sep = "\t")
    }

    if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
      fwrite(as.data.frame(ekk),
             file.path("results/full_enrichment", paste0(prefix, "_KEGG_3k.tsv")),
             sep = "\t")
    }

    TRUE
  }, error = function(e) {
    message("  [Warning] Enrichment failed (", prefix, "): ", conditionMessage(e))
    FALSE
  })
}

# -----------------------------------------------------------------------------
# Procesamiento por dataset
# -----------------------------------------------------------------------------
process_dataset <- function(acc, workers_per_dataset = 1L) {
  message("=====================================")
  message("=== FULL drivers: ", acc, " ===")

  ord <- state_orders[[acc]]

  nets <- list()
  betas <- list()

  for (st in ord) {
    fp <- file.path("results/full_networks", paste0(acc, "_", st, "_full_network.rds"))
    if (file.exists(fp)) {
      obj <- readRDS(fp)
      nets[[st]] <- if (is.list(obj) && !is.null(obj$W)) obj$W else obj
      betas[[st]] <- if (is.list(obj) && !is.null(obj$beta)) obj$beta else 6
    }
  }

  if (length(nets) < 2L) {
    message("  Saltando ", acc, ": menos de 2 estados con red.")
    return(NULL)
  }

  # Reducción consistente
  nets <- reduce_networks_consistent(nets, max_genes = MAX_GENES_NODE)
  if (is.null(nets)) {
    message("  Saltando ", acc, ": no hay genes comunes.")
    return(NULL)
  }

  message("  Reducción consistente: ", nrow(nets[[1]]), " genes")

  # Métricas de nodo
  message("  Métricas de nodo rápidas...")
  nm <- lapply(names(nets), function(st) {
    node_metrics_fast(nets[[st]])
  })
  names(nm) <- names(nets)

  base <- Reduce(
    function(x, y) full_join(x, y, by = "gene", suffix = c("", "")),
    Map(function(df, nm_) rename_with(df, ~ paste0(., "_", nm_), -gene), nm, names(nm))
  )

  first <- names(nets)[1]
  last  <- names(nets)[length(nets)]
  mid   <- if (length(nets) >= 3L) names(nets)[2] else last

  df <- base |>
    mutate(
      PTI =
        zscore(abs(.data[[paste0("strength_", mid)]]      - .data[[paste0("strength_", first)]])) +
        zscore(abs(.data[[paste0("eigen_", mid)]]         - .data[[paste0("eigen_", first)]])) +
        zscore(abs(.data[[paste0("participation_", mid)]] - .data[[paste0("participation_", first)]])),
      IRI =
        zscore(abs(.data[[paste0("strength_", last)]]      - .data[[paste0("strength_", first)]])) +
        zscore(abs(.data[[paste0("eigen_", last)]]         - .data[[paste0("eigen_", first)]])) +
        zscore(abs(.data[[paste0("participation_", last)]] - .data[[paste0("participation_", first)]]))
    )

  W_ref <- nets[[first]]
  W_dis <- nets[[last]]
  common <- intersect(rownames(W_ref), rownames(W_dis))

  message("  TRI y KO-support (", length(common), " genes comunes; límite TRI=", MAX_GENES_TRI, ")...")

  tri_vec <- compute_TRI_exact_hpc(W_dis, W_ref, common, n_workers_inner = workers_per_dataset)
  ko_vec  <- compute_KO_exact_hpc(W_dis, common, n_workers_inner = workers_per_dataset)

  out <- df |>
    inner_join(data.frame(gene = names(tri_vec), TRI = zscore(tri_vec)), by = "gene") |>
    inner_join(data.frame(gene = names(ko_vec),  KO_support = ko_vec), by = "gene") |>
    mutate(
      class = case_when(
        PTI >= quantile(PTI, 0.90, na.rm = TRUE) &
          TRI >= quantile(TRI, 0.90, na.rm = TRUE) ~ "transition_rescue",

        IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
          TRI >= quantile(TRI, 0.90, na.rm = TRUE) ~ "irreversibility_rescue",

        IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
          KO_support >= quantile(KO_support, 0.90, na.rm = TRUE) ~ "irreversibility_lock",

        PTI >= quantile(PTI, 0.90, na.rm = TRUE) ~ "transition_driver",
        TRUE ~ "other"
      )
    ) |>
    arrange(desc(PTI + IRI + TRI))

  fwrite(
    out,
    file.path("results/full_gene_drivers", paste0(acc, "_full_gene_driver_scores.tsv")),
    sep = "\t"
  )

  # Tops
  top_trans_100 <- slice_head(arrange(out, desc(PTI)), n = 100)
  top_irrev_100 <- slice_head(arrange(out, desc(IRI)), n = 100)
  top_resc_100  <- slice_head(arrange(out, desc(TRI)), n = 100)

  fwrite(top_trans_100,
         file.path("results/full_gene_drivers", paste0(acc, "_full_top_transition.tsv")),
         sep = "\t")
  fwrite(top_irrev_100,
         file.path("results/full_gene_drivers", paste0(acc, "_full_top_irreversibility.tsv")),
         sep = "\t")
  fwrite(top_resc_100,
         file.path("results/full_gene_drivers", paste0(acc, "_full_top_rescue.tsv")),
         sep = "\t")

  # Enriquecimiento
  top_trans_enrich <- slice_head(arrange(out, desc(PTI)), n = min(N_TOP_ENRICH, nrow(out)))
  top_irrev_enrich <- slice_head(arrange(out, desc(IRI)), n = min(N_TOP_ENRICH, nrow(out)))
  top_resc_enrich  <- slice_head(arrange(out, desc(TRI)), n = min(N_TOP_ENRICH, nrow(out)))

  universe_genes <- out$gene

  message("  Enriquecimiento...")
  safe_enrich_robust(top_trans_enrich$gene, universe_genes, paste0(acc, "_full_transition"))
  safe_enrich_robust(top_irrev_enrich$gene, universe_genes, paste0(acc, "_full_irreversibility"))
  safe_enrich_robust(top_resc_enrich$gene,  universe_genes, paste0(acc, "_full_rescue"))

  # Recurrentes globales
  bind_rows(
    mutate(top_trans_100, score_type = "PTI", accession = acc),
    mutate(top_irrev_100, score_type = "IRI", accession = acc),
    mutate(top_resc_100,  score_type = "TRI", accession = acc)
  )
}

# -----------------------------------------------------------------------------
# Ejecución paralela por dataset
# -----------------------------------------------------------------------------
n_datasets <- length(targets)
mc_cores_outer <- min(n_datasets, n_workers_total)
workers_per_dataset <- max(1L, floor(n_workers_total / mc_cores_outer))

message("Datasets en paralelo: ", mc_cores_outer,
        " | workers internos por dataset: ", workers_per_dataset)

all_top_list <- mclapply(
  targets,
  function(acc) {
    tryCatch(
      process_dataset(acc, workers_per_dataset = workers_per_dataset),
      error = function(e) {
        message("ERROR en dataset ", acc, ": ", conditionMessage(e))
        NULL
      }
    )
  },
  mc.cores = mc_cores_outer,
  mc.preschedule = FALSE
)

all_top <- all_top_list[!vapply(all_top_list, is.null, logical(1))]

if (length(all_top) > 0L) {
  recurrent <- bind_rows(all_top) |>
    distinct(accession, gene, score_type) |>
    count(gene, score_type, name = "n_datasets") |>
    arrange(desc(n_datasets), gene)

  fwrite(
    recurrent,
    "results/full_gene_drivers/recurrent_full_driver_genes_across_datasets.tsv",
    sep = "\t"
  )
}

message("Done.")
