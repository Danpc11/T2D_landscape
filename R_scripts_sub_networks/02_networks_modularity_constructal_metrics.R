#!/usr/bin/env Rscript
# =============================================================================
# 02_networks_modularity_constructal_metrics.R
# Construye redes de co-expresion por estado sobre el espacio COMUN de genes
# (interseccion de los 4 datasets de microarray), calcula metricas constructales
# y detecta modulos.
#
# DATASETS: GSE76895, GSE18732, GSE15653, GSE27951
# RAMA: Compartida (genes fijos, metricas directamente comparables entre datasets)
#
# CORRECCIONES:
#   [C5]  Subred fija de genes comunes filtrada por varianza global (script 01).
#   [C6]  CEI calculado globalmente sobre todos los datasets/estados juntos.
#   [C7]  find_modules maximiza n_mod * Q (modularidad Newman) en lugar de solo
#         n_mod. Fallback a corte mediano si ningun corte produce modulos validos.
#   [C8]  Orden fisiologico en archivos de salida.
#   [FIX-Rbar-v2] Rbar = 2*tr(L+)/(p-1) (media de R_ij por par; Kirchhoff Kf = p*tr(L+)).
#   [FIX-EG]      EG = media de 1/d_ij sobre pares (se elimino un *2 espurio).
#   [FIX-beta-dataset] beta unico por dataset (estado sano).
#   [FIX-subset]  bicor solo sobre la subred fija (antes: 5000 genes y recorte).
#   [NEW] metricas a n igual (submuestreo) y relativas al nulo de configuracion.
#   [FIX-match] match(states, condition) en lugar del match invertido anterior.
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(igraph)
  library(expm)
  library(Matrix)
})

options(stringsAsFactors = FALSE)
allowWGCNAThreads()

dir.create("results/networks", recursive = TRUE, showWarnings = FALSE)
dir.create("results/metrics",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/modules",  recursive = TRUE, showWarnings = FALSE)

args    <- commandArgs(trailingOnly = TRUE)
targets <- if (length(args) == 0) c("GSE76895", "GSE18732", "GSE15653", "GSE27951") else args

# Canonical physiological order
state_orders <- list(
  GSE76895 = c("ND", "IGT", "T2D"),
  GSE18732 = c("ND", "IGT", "T2D"),
  GSE15653 = c("Lean", "Obese_noT2D", "Obese_T2D"),
  GSE27951 = c("NGT", "IGT", "T2D")
)

# -----------------------------------------------------------------------------
# Network functions
# -----------------------------------------------------------------------------

# [FIX-beta-dataset] beta se estima UNA vez por dataset (estado sano) y se
# reutiliza en todos los estados. Con beta distinto por estado, W cambia de
# escala y las metricas comparan beta, no biologia.
pick_beta <- function(expr, cor_method = "bicor") {
  expr_t <- t(expr)
  sft    <- pickSoftThreshold(expr_t, dataIsExpr = TRUE, corFnc = cor_method,
                              corOptions = "use = 'p', maxPOutliers = 0.1",
                              powerVector = 1:20, RsquaredCut = 0.80, verbose = 0)
  beta <- sft$powerEstimate
  if (is.na(beta)) {
    fit       <- sft$fitIndices
    candidate <- fit$Power[fit$SFT.R.sq >= 0.80]
    beta      <- if (length(candidate) > 0) min(candidate) else 6L
  }
  beta
}

compute_adjacency <- function(expr, beta, cor_method = "bicor") {
  expr_t  <- t(expr)
  cor_mat <- if (cor_method == "bicor") bicor(expr_t, maxPOutliers = 0.1)
             else cor(expr_t, method = "pearson")
  cor_mat[is.na(cor_mat)] <- 0
  adj           <- abs(cor_mat)^beta
  diag(adj)     <- 0
  rownames(adj) <- colnames(adj) <- rownames(expr)
  list(adj = adj, beta = beta)
}

# [NEW] Nulo de configuracion (Chung-Lu ponderado): misma secuencia de fuerzas,
# sin estructura. Las metricas se reportan tambien RELATIVAS a este nulo para
# separar estructura de densidad.
configuration_null <- function(W) {
  k  <- rowSums(W); m2 <- sum(k)
  W0 <- outer(k, k) / m2
  diag(W0) <- 0
  # reescalar para conservar exactamente la suma de pesos
  W0 * (sum(W) / sum(W0))
}

# [C7] Seleccion de corte maximizando n_mod * Q (no solo n_mod)
find_modules <- function(W, min_size = 20L) {
  d  <- as.dist(1 - W)
  hc <- hclust(d, method = "average")
  g  <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                    diag = FALSE)
  best <- NULL
  for (q in seq(0.55, 0.80, by = 0.05)) {
    h       <- quantile(hc$height, probs = q, na.rm = TRUE)
    cl      <- cutree(hc, h = h)
    tab     <- table(cl)
    keep_cl <- names(tab)[tab >= min_size]
    n_mod   <- length(keep_cl)
    if (n_mod < 1L) next
    module_tmp        <- paste0("M", cl)
    small             <- names(table(module_tmp))[table(module_tmp) < min_size]
    module_tmp[module_tmp %in% small] <- "grey"
    mem     <- as.integer(factor(module_tmp))
    Q_val   <- tryCatch(
      modularity(g, membership = mem, weights = E(g)$weight),
      error = function(e) 0
    )
    score <- n_mod * max(Q_val, 0)
    candidate <- list(h = h, cl = cl, n_mod = n_mod, Q = Q_val, score = score)
    if (is.null(best) || candidate$score > best$score) best <- candidate
  }
  if (is.null(best)) {
    # Fallback: corte mediano si ningun percentil produce modulos de tamano >= min_size
    h    <- median(hc$height)
    cl   <- cutree(hc, h = h)
    best <- list(cl = cl, n_mod = 1L, Q = 0, score = 0)
  }
  module        <- paste0("M", best$cl)
  tab           <- table(module)
  small         <- names(tab)[tab < min_size]
  module[module %in% small] <- "grey"
  names(module) <- rownames(W)
  list(module    = module,
       n_modules = length(setdiff(unique(module), "grey")),
       Q         = best$Q)
}

# -----------------------------------------------------------------------------
# Metricas constructales
# [FIX-Rbar-v2] Rbar = 2*tr(L+)/(p-1): media de la resistencia efectiva por par
# -----------------------------------------------------------------------------

compute_strength_entropy <- function(W) {
  k <- rowSums(W)
  p <- k / sum(k)
  -sum(p * log(p + 1e-12))
}

compute_global_efficiency <- function(W) {
  g    <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE,
                                      diag = FALSE)
  D    <- distances(g, weights = 1 / (E(g)$weight + 1e-6))
  invD <- 1 / D
  diag(invD) <- NA_real_
  mean(invD[upper.tri(invD)], na.rm = TRUE)
}

compute_communicability <- function(W) {
  k    <- rowSums(W)
  Dinv <- diag(1 / sqrt(k + 1e-12), nrow(W))
  G    <- expm(Dinv %*% W %*% Dinv)
  mean(G[row(G) != col(G)])
}

compute_avg_effective_resistance <- function(W) {
  p        <- nrow(W)
  k        <- rowSums(W)
  L        <- diag(k, p) - W
  eig      <- eigen(L, symmetric = TRUE)
  tol      <- max(abs(eig$values)) * p * 1e-10
  inv_vals <- ifelse(eig$values > tol, 1 / eig$values, 0)
  Lplus    <- eig$vectors %*% diag(inv_vals, p) %*% t(eig$vectors)
  as.numeric(2 * sum(diag(Lplus)) / (p - 1))   # [FIX-Rbar-v2] media por par
}

# Wrapper unico — mismo en scripts 03 y 04 para consistencia
compute_metrics <- function(W) {
  c(EG   = compute_global_efficiency(W),
    Gbar = compute_communicability(W),
    Rbar = compute_avg_effective_resistance(W),
    Hb   = compute_strength_entropy(W))
}

# -----------------------------------------------------------------------------
# [C5] Carga del conjunto fijo de genes comunes (generado por script 01)
# -----------------------------------------------------------------------------

hv_genes_path <- "data/processed/high_variance_genes.rds"
if (!file.exists(hv_genes_path)) {
  stop("'high_variance_genes.rds' no encontrado. Ejecuta primero el script 01.")
}
hv_genes    <- readRDS(hv_genes_path)
p_sub       <- min(800L, length(hv_genes))
fixed_genes <- hv_genes[seq_len(p_sub)]
message("Fixed subnetwork nodes: ", length(fixed_genes))

# -----------------------------------------------------------------------------
# Bucle principal
# -----------------------------------------------------------------------------

metrics_list <- list()
n_sub_reps   <- 20L   # [NEW] replicas de submuestreo a n igual entre estados

for (acc in targets) {
  obj   <- readRDS(file.path("data/processed", paste0(acc, "_processed.rds")))
  pheno <- obj$pheno
  # [FIX-subset] recortar a la subred fija ANTES de correlacionar
  genes_use <- intersect(fixed_genes, rownames(obj$expr))
  expr      <- obj$expr[genes_use, , drop = FALSE]
  message("=== ", acc, ": ", nrow(expr), " genes en subred fija ===")

  ord    <- state_orders[[acc]]
  states <- ord[ord %in% as.character(unique(pheno$condition))]
  n_by_state <- sapply(states, function(st) sum(pheno$condition == st))
  states <- states[n_by_state[states] >= 4L]
  n_min  <- min(n_by_state[states])
  message("  n por estado: ", paste(states, n_by_state[states], sep = "=", collapse = " | "),
          " -> submuestreo a n_min=", n_min)

  # beta fijo del estado sano (primer estado en orden fisiologico)
  healthy_samples <- pheno$.sample_id[pheno$condition == states[1]]
  beta_acc <- pick_beta(expr[, healthy_samples, drop = FALSE])
  message("  beta (estado sano, fijo para el dataset) = ", beta_acc)

  acc_metrics <- list()
  for (st in states) {
    samples  <- pheno$.sample_id[pheno$condition == st]
    expr_sub <- expr[, samples, drop = FALSE]

    # Red con TODAS las muestras del estado (se guarda para scripts 03/04)
    Wsub <- compute_adjacency(expr_sub, beta_acc)$adj
    mods <- find_modules(Wsub)
    fwrite(data.frame(gene = names(mods$module), module = unname(mods$module),
                      accession = acc, condition = st),
           file.path("results/modules", paste0(acc, "_", st, "_modules.tsv")), sep = "\t")
    m    <- compute_metrics(Wsub)
    m0   <- compute_metrics(configuration_null(Wsub))

    # [NEW] Metricas a n igual: media/sd sobre submuestras de tamano n_min
    set.seed(1234)
    sub_m <- t(replicate(n_sub_reps, {
      sel <- sample(samples, n_min)
      compute_metrics(compute_adjacency(expr_sub[, sel, drop = FALSE], beta_acc)$adj)
    }))

    acc_metrics[[st]] <- data.frame(
      accession = acc, condition = st,
      n_samples = ncol(expr_sub), n_sub = n_min, n_subgenes = nrow(Wsub), beta = beta_acc,
      EG = m["EG"], Gbar = m["Gbar"], Rbar = m["Rbar"], Hb = m["Hb"],
      EG_rel = m["EG"] / m0["EG"], Gbar_rel = m["Gbar"] / m0["Gbar"],
      Rbar_rel = m["Rbar"] / m0["Rbar"], Hb_rel = m["Hb"] / m0["Hb"],
      EG_sub = mean(sub_m[, "EG"]),   EG_sub_sd = sd(sub_m[, "EG"]),
      Gbar_sub = mean(sub_m[, "Gbar"]), Gbar_sub_sd = sd(sub_m[, "Gbar"]),
      Rbar_sub = mean(sub_m[, "Rbar"]), Rbar_sub_sd = sd(sub_m[, "Rbar"]),
      Hb_sub = mean(sub_m[, "Hb"]),   Hb_sub_sd = sd(sub_m[, "Hb"]),
      mean_k = mean(rowSums(Wsub)), n_modules = mods$n_modules, Q_modularity = mods$Q
    )
    saveRDS(list(W = Wsub, genes = rownames(Wsub), beta = beta_acc),
            file.path("results/networks", paste0(acc, "_", st, "_network.rds")))
  }
  metrics_df <- dplyr::bind_rows(acc_metrics)
  metrics_df <- metrics_df[match(states, metrics_df$condition), ]
  metrics_list[[acc]] <- metrics_df
}

# [C6] CEI calculado GLOBALMENTE sobre todos los datasets/estados
all_metrics <- dplyr::bind_rows(metrics_list)

zscore_global <- function(x) as.numeric(scale(x))

all_metrics <- all_metrics |>
  dplyr::mutate(
    # CEI sobre metricas a n IGUAL (evita confundir densidad por ruido de muestreo)
    CEI = zscore_global(EG_sub) + zscore_global(Gbar_sub) -
          zscore_global(Rbar_sub) + zscore_global(Hb_sub),
    # CEI relativo al nulo de configuracion (estructura, no densidad)
    # (Hb_rel == 1 por construccion: el nulo preserva las fuerzas; se excluye)
    CEI_rel = zscore_global(EG_rel) + zscore_global(Gbar_rel) -
              zscore_global(Rbar_rel)
  )

# Reescribir archivos por dataset con CEI global incluido
for (acc in targets) {
  df_acc <- all_metrics[all_metrics$accession == acc, ]
  fwrite(df_acc,
         file.path("results/metrics",
                   paste0(acc, "_constructal_metrics.tsv")), sep = "\t")
}

fwrite(all_metrics, "results/metrics/all_constructal_metrics.tsv", sep = "\t")
message("Done.")
