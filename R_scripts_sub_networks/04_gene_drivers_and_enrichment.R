#!/usr/bin/env Rscript
# =============================================================================
# 04_gene_drivers_and_enrichment.R
# Compute gene-driver scores (PTI, IRI, TRI, KO-support) and perform
# functional enrichment (GO:BP and KEGG) in the SHARED gene space.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# BRANCH: Shared
#
# CORRECTIONS:
#   [C13] IRI includes the "participation" dimension (as PTI does).
#   [C14] TRI documented as a computational estimate of rescue potential.
#   [C15] compute_metrics identical to scripts 02 and 03.
#   [C16] Datasets updated: GSE76895, GSE18732, GSE15653, GSE27951.
#   [FIX-TRI-exact] TRI y KO por recalculo exacto (la "actualizacion rank-1"
#         anterior era matematicamente invalida; ver comentario en la funcion).
#   [FIX-dorng] paralelismo funciona con o sin doRNG.
#   [FIX-enrich] GO/KEGG independientes, con universo = subred.
#   NOTA: PTI/IRI/TRI/KO son medidas de SENSIBILIDAD topologica de la red de
#         correlacion al nodo, no evidencia causal de "rescate" o "knock-out".
#   [FIX-PTI] Two-condition datasets: mid = last -> PTI != IRI by
#         construction. The intermediate state is used as mid when available,
#         otherwise a warning is issued.
#   [FIX-EG] EG sin el *2 espurio.
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(igraph)
  library(parallel)
  library(foreach)
  library(doParallel)
})

# doRNG opcional — reproducibilidad en paralelo
HAS_DORNG <- requireNamespace("doRNG", quietly = TRUE)
if (HAS_DORNG) {
  library(doRNG)
  message("doRNG disponible: resultados reproducibles.")
} else {
  message("doRNG no disponible. Instalar con: install.packages('doRNG')")
}

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Parallelization (CRITICAL)
# -----------------------------------------------------------------------------

n_workers <- as.integer(Sys.getenv("N_WORKERS", 40))

# [CORRECTED] Prevents oversubscription in matrix multiplications (%*%)
Sys.setenv(
  OMP_NUM_THREADS = 1,
  OPENBLAS_NUM_THREADS = 1,
  MKL_NUM_THREADS = 1
)

RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

# [CORRECTED] Explicit cluster creation (PSOCK)
cl <- makeCluster(n_workers)
registerDoParallel(cl)
if (HAS_DORNG) registerDoRNG(1234)
# [NEW] Enriquecimiento opcional: si faltan clusterProfiler/org.Hs.eg.db se omite
# (los scores de drivers se calculan igual).
HAS_ENRICH <- requireNamespace("clusterProfiler", quietly = TRUE) &&
              requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (HAS_ENRICH) suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
if (!HAS_ENRICH) message("clusterProfiler/org.Hs.eg.db no disponibles: se omite el enriquecimiento")
# [FIX-dorng] %dorng% solo existe si doRNG esta instalado; antes el script fallaba sin el.
`%dop%` <- if (HAS_DORNG) doRNG::`%dorng%` else foreach::`%dopar%`

message("Parallel engine: ", n_workers, " workers | BLAS=1 | doRNG=", HAS_DORNG)

# -----------------------------------------------------------------------------

dir.create("results/gene_drivers", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment",   recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

state_orders <- list(
  GSE76895 = c("ND", "IGT", "T2D"),
  GSE18732 = c("ND", "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

zscore <- function(x) {
  s <- scale(x); s[is.na(s)] <- 0; as.numeric(s)
}

# -----------------------------------------------------------------------------
# Network metrics—identical to scripts 02 and 03
# [FIX-Rbar] Rbar = 2*tr(L+)/p
# -----------------------------------------------------------------------------

compute_metrics <- function(W) {
  p    <- nrow(W)
  k    <- rowSums(W)
  g    <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
  D    <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / D; diag(invD) <- NA_real_
  EG   <- mean(invD[upper.tri(invD)], na.rm = TRUE)
  prob <- k / sum(k)
  Hb   <- -sum(prob * log(prob + 1e-12))
  c(EG = EG, Hb = Hb)
}

# -----------------------------------------------------------------------------
# Node metrics
# -----------------------------------------------------------------------------

node_metrics <- function(W) {
  g   <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                     diag = FALSE)
  k   <- strength(g)
  ev  <- eigen_centrality(g, weights = E(g)$weight)$vector
  btw <- betweenness(g, weights = 1 / (E(g)$weight + 1e-6), normalized = TRUE)
  cl  <- cluster_fast_greedy(g)
  mem <- membership(cl)
  Pc  <- sapply(seq_along(k), function(i) {
    ki <- k[i]
    if (ki == 0) return(0)
    kim <- tapply(W[i, ], mem, sum)
    1 - sum((kim / ki)^2, na.rm = TRUE)
  })
  data.frame(gene = names(k), strength = k, eigen = ev,
             betweenness = btw, participation = Pc)
}

# -----------------------------------------------------------------------------
# [FIX-TRI-exact] TRI y KO por RECALCULO EXACTO.
# La version anterior usaba una "actualizacion rank-1 (Sherman-Morrison)" de tr(L+).
# Eso no es valido: cambiar la fila/columna i de W cambia L en
#   dL = diag(dk) - dW, que NO es rango 1, y la expresion usada tampoco era la
# expansion de primer orden -tr(L+ dL L+) (verificado numericamente: error ~30x).
# Con p <= 800 el recalculo exacto de EG (distancias) y Hb cuesta ~0.2 s/gen.
# TRI_i = [EG(W_dis con fila/col i de W_ref) - EG(W_dis)] + [Hb(...) - Hb(W_dis)]
# KO_i  = EG(W_dis) - EG(W_dis con fila/col i = 0)
# -----------------------------------------------------------------------------
rescue_row <- function(W_dis, W_ref, i) {
  W <- W_dis
  W[i, ] <- W_ref[i, ]; W[, i] <- W_ref[, i]; W[i, i] <- 0
  W
}
knockout_row <- function(W_dis, i) {
  W <- W_dis
  W[i, ] <- 0; W[, i] <- 0
  W
}

compute_TRI_exact <- function(W_dis, W_ref, common) {
  W_dis <- W_dis[common, common, drop = FALSE]
  W_ref <- W_ref[common, common, drop = FALSE]
  m0    <- compute_metrics(W_dis)
  vals  <- foreach(i = seq_along(common), .combine = "c",
                   .packages = "igraph",
                   .export = c("compute_metrics", "rescue_row")) %dop% {
    m <- compute_metrics(rescue_row(W_dis, W_ref, i))
    (m["EG"] - m0["EG"]) + (m["Hb"] - m0["Hb"])
  }
  names(vals) <- common
  vals
}

compute_KO_exact <- function(W_dis, common) {
  W_dis <- W_dis[common, common, drop = FALSE]
  m0    <- compute_metrics(W_dis)
  vals  <- foreach(i = seq_along(common), .combine = "c",
                   .packages = "igraph",
                   .export = c("compute_metrics", "knockout_row")) %dop% {
    m <- compute_metrics(knockout_row(W_dis, i))
    unname(m0["EG"] - m["EG"])
  }
  names(vals) <- common
  vals
}

# -----------------------------------------------------------------------------
# Enriquecimiento funcional
# -----------------------------------------------------------------------------

# [FIX-enrich] GO y KEGG en tryCatch SEPARADOS: si KEGG falla (sin internet en HPC)
# antes se perdia tambien el GO ya calculado. Universo = genes de la subred.
enrich_gene_set <- function(genes, universe, prefix) {
  if (!HAS_ENRICH) return(invisible(NULL))
  eg <- tryCatch(suppressMessages(bitr(genes, fromType = "SYMBOL",
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)),
                 error = function(e) NULL)
  un <- tryCatch(suppressMessages(bitr(universe, fromType = "SYMBOL",
                                       toType = "ENTREZID", OrgDb = org.Hs.eg.db)),
                 error = function(e) NULL)
  if (is.null(eg) || nrow(eg) < 10L) return(invisible(NULL))
  tryCatch({
    ego <- suppressMessages(enrichGO(eg$ENTREZID, universe = un$ENTREZID,
                                     OrgDb = org.Hs.eg.db, ont = "BP",
                                     pAdjustMethod = "BH", readable = TRUE))
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0)
      fwrite(as.data.frame(ego),
             file.path("results/enrichment", paste0(prefix, "_GO.tsv")), sep = "\t")
  }, error = function(e) message("  GO error (", prefix, "): ", e$message))
  tryCatch({
    ekk <- suppressMessages(enrichKEGG(eg$ENTREZID, universe = un$ENTREZID,
                                       organism = "hsa", pAdjustMethod = "BH"))
    if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0)
      fwrite(as.data.frame(ekk),
             file.path("results/enrichment", paste0(prefix, "_KEGG.tsv")), sep = "\t")
  }, error = function(e) message("  KEGG error (", prefix, "): ", e$message,
                                  " (requiere acceso a rest.kegg.jp)"))
}

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

all_top <- list()

for (acc in targets) {
  message("=== ", acc, " ===")
  ord <- state_orders[[acc]]

  nets <- lapply(setNames(ord, ord), function(st) {
    fp <- file.path("results/networks", paste0(acc, "_", st, "_network.rds"))
    if (!file.exists(fp)) return(NULL)
    readRDS(fp)$W
  })
  nets <- nets[!vapply(nets, is.null, logical(1))]

  if (length(nets) < 2L) {
    message("  Menos de 2 redes disponibles; saltando.")
    next
  }

  message("  Calculando metricas de nodo...")
  nm   <- lapply(nets, node_metrics)
  base <- Reduce(function(x, y) full_join(x, y, by = "gene", suffix = c("", "")),
                 Map(function(df, nm_) rename_with(df, ~ paste0(., "_", nm_), -gene),
                     nm, names(nm)))

  first <- names(nets)[1]
  last  <- names(nets)[length(nets)]

  # [FIX-PTI] mid = estado intermedio si existe, de lo contrario = last
  if (length(nets) >= 3L) {
    mid <- names(nets)[2]
  } else {
    mid <- last
    message("  AVISO: solo 2 estados; PTI = IRI (no hay estado intermedio)")
  }

  df <- base |>
    dplyr::mutate(
      PTI = zscore(abs(.data[[paste0("strength_",     mid)]] - .data[[paste0("strength_",     first)]])) +
            zscore(abs(.data[[paste0("eigen_",        mid)]] - .data[[paste0("eigen_",        first)]])) +
            zscore(abs(.data[[paste0("betweenness_",  mid)]] - .data[[paste0("betweenness_",  first)]])) +
            zscore(abs(.data[[paste0("participation_", mid)]] - .data[[paste0("participation_", first)]])),
      IRI = zscore(abs(.data[[paste0("strength_",     last)]] - .data[[paste0("strength_",     first)]])) +
            zscore(abs(.data[[paste0("eigen_",        last)]] - .data[[paste0("eigen_",        first)]])) +
            zscore(abs(.data[[paste0("betweenness_",  last)]] - .data[[paste0("betweenness_",  first)]])) +
            zscore(abs(.data[[paste0("participation_", last)]] - .data[[paste0("participation_", first)]]))
    )

  W_ref  <- nets[[first]]
  W_dis  <- nets[[last]]
  common <- intersect(rownames(W_ref), rownames(W_dis))

  message("  TRI (recalculo exacto, ", length(common), " genes)...")
  tri_vec <- compute_TRI_exact(W_dis, W_ref, common)

  message("  KO-support (recalculo exacto, ", length(common), " genes)...")
  ko_vec  <- compute_KO_exact(W_dis, common)

  out <- df |>
    inner_join(data.frame(gene = names(tri_vec), TRI = zscore(tri_vec)),
               by = "gene") |>
    inner_join(data.frame(gene = names(ko_vec),  KO_support = ko_vec),
               by = "gene") |>
    dplyr::mutate(class = dplyr::case_when(
      PTI >= quantile(PTI, 0.90, na.rm = TRUE) &
        TRI >= quantile(TRI, 0.90, na.rm = TRUE)         ~ "transition_rescue",
      IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
        TRI >= quantile(TRI, 0.90, na.rm = TRUE)         ~ "irreversibility_rescue",
      IRI >= quantile(IRI, 0.90, na.rm = TRUE) &
        KO_support >= quantile(KO_support, 0.90, na.rm = TRUE) ~ "irreversibility_lock",
      PTI >= quantile(PTI, 0.90, na.rm = TRUE)            ~ "transition_driver",
      TRUE ~ "other"
    )) |>
    dplyr::arrange(dplyr::desc(PTI + IRI + TRI))

  fwrite(out,
         file.path("results/gene_drivers",
                   paste0(acc, "_gene_driver_scores.tsv")), sep = "\t")

  top_trans <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(PTI)), n = 50)
  top_irrev <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(IRI)), n = 50)
  top_resc  <- dplyr::slice_head(dplyr::arrange(out, dplyr::desc(TRI)), n = 50)

  fwrite(top_trans,
         file.path("results/gene_drivers", paste0(acc, "_top_transition.tsv")),
         sep = "\t")
  fwrite(top_irrev,
         file.path("results/gene_drivers", paste0(acc, "_top_irreversibility.tsv")),
         sep = "\t")
  fwrite(top_resc,
         file.path("results/gene_drivers", paste0(acc, "_top_rescue.tsv")),
         sep = "\t")

  message("  Enriquecimiento funcional...")
  enrich_gene_set(top_trans$gene, out$gene, paste0(acc, "_transition"))
  enrich_gene_set(top_irrev$gene, out$gene, paste0(acc, "_irreversibility"))
  enrich_gene_set(top_resc$gene,  out$gene, paste0(acc, "_rescue"))

  all_top[[acc]] <- dplyr::bind_rows(
    dplyr::mutate(top_trans, score_type = "PTI", accession = acc),
    dplyr::mutate(top_irrev, score_type = "IRI", accession = acc),
    dplyr::mutate(top_resc,  score_type = "TRI", accession = acc)
  )
}

if (length(all_top) > 0L) {
  recurrent <- dplyr::bind_rows(all_top) |>
    dplyr::distinct(accession, gene, score_type) |>
    dplyr::count(gene, score_type, name = "n_datasets") |>
    dplyr::arrange(dplyr::desc(n_datasets), gene)
  fwrite(recurrent,
         "results/gene_drivers/recurrent_driver_genes_across_datasets.tsv",
         sep = "\t")
}

# [CORRECTED] Cierre de clúster explícito
stopCluster(cl)
message("Done.")
