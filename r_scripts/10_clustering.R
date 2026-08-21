# Computes a neighbor graph/UMAP and Leiden-clusters one tissue's
# harmony-integrated data across a range of resolutions, scoring each with
# approximate silhouette width, neighborhood purity, and graph modularity,
# for the user to pick a resolution to carry forward. Runs as a SLURM job
# array (see jobs/10_clustering.sh), one task per tissue, since all 3
# tissues are fully independent (same restructuring as 04_doubletfinder.R,
# 07_norm_pca.R, and 08_sketch_pca.R). 11_clustering.R (FindAllMarkers())
# only needs cluster labels plus the expression data already saved by
# 09_*_integration.R, which this script doesn't touch -- so metadata.rds
# (with one "resX_clusters" column per tested resolution) is the only
# required output. The harmony UMAP is also saved for convenience/
# plotting, but isn't needed by 11.
#
# Cluster-quality scoring uses centroid-approximated silhouette width and
# a kNN-graph-based neighborhood purity/modularity, rather than
# cluster::silhouette() on a full cell x cell distance matrix -- that
# O(cells^2) matrix OOM-killed this job on the larger tissues. Both
# approximations only need cell-to-centroid or cell-to-neighbor
# comparisons (O(cells x clusters) or reusing the already-computed
# neighbor graph), never a full pairwise distance matrix.

# Load libraries
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(ggbeeswarm)
  library(ggalluvial)
  library(igraph)
  library(reticulate)
  library(BPCells)
})

options(reticulate.verbose = TRUE)
# py_module_available(module = "leidenalg")

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Cluster-quality helpers ----------------------------------------------------

# Approximates cluster::silhouette() without ever building a cell x cell
# distance matrix. For each cell, computes distance to its own cluster's
# centroid (a) vs. distance to the nearest other cluster's centroid (b),
# and returns (b - a) / max(a, b) -- same range/interpretation as classic
# silhouette width. Memory scales with cells x clusters, not cells^2.
approx_silhouette <- function(embedding, clusters){
  clusters <- as.character(clusters)
  cluster_levels <- sort(unique(clusters))
  k <- length(cluster_levels)

  centroids <- rowsum(embedding, clusters)
  centroids <- sweep(centroids, 1, as.vector(table(clusters)[cluster_levels]), "/")

  dist_to_centroids <- matrix(NA_real_, nrow = nrow(embedding), ncol = k,
                              dimnames = list(NULL, cluster_levels))
  for (j in seq_len(k)){
    dist_to_centroids[, j] <- sqrt(rowSums(sweep(embedding, 2, centroids[j, ]) ^ 2))
  }

  own_idx <- match(clusters, cluster_levels)
  a <- dist_to_centroids[cbind(seq_len(nrow(embedding)), own_idx)]

  dist_to_others <- dist_to_centroids
  dist_to_others[cbind(seq_len(nrow(embedding)), own_idx)] <- Inf
  nearest_other <- apply(dist_to_others, 1, which.min)
  b <- dist_to_others[cbind(seq_len(nrow(embedding)), nearest_other)]

  data.frame(cluster = clusters,
             neighbor = cluster_levels[nearest_other],
             sil_width = (b - a) / pmax(a, b))
}

# For each cell, the fraction of its k nearest neighbors that share its
# cluster label. Reuses the Neighbor object from the first FindNeighbors()
# call below (an N x k.param matrix of neighbor indices) instead of
# computing any new distances.
neighbor_purity <- function(obj, clusters, neighbor_name = "RNA.nn"){
  idx <- Indices(obj[[neighbor_name]])
  neighbor_clusters <- matrix(clusters[idx], nrow = nrow(idx))
  rowMeans(neighbor_clusters == clusters)
}

# The same objective Leiden already optimized when it clustered. Reuses
# the SNN graph from the second FindNeighbors() call below -- it's already
# a sparse Graph/dgCMatrix, so igraph builds from it without densifying.
graph_modularity <- function(obj, clusters, graph_name = "RNA_snn"){
  snn <- as(obj[[graph_name]], "dgCMatrix")
  g <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::modularity(g, membership = as.integer(factor(clusters)))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/10_clustering.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/10_clustering.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

data_in_dir <- paste0("data/09_integration/", tissue_file, "/")

plots_dir <- paste0("plots/10_clustering/", tissue_file, "/")
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- paste0("data/10_clustering/", tissue_file, "/")
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

tab_data_out_dir <- paste0("tab_data/10_clustering/", tissue_file, "/")
dir.create(tab_data_out_dir, showWarnings = F,
           recursive = T)

# Load the harmony-integrated object from 09 ---------------------------------

message2("Reading in integrated object")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
harmony <- readRDS(paste0(data_in_dir, "harmony.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony

# Compute the NN graph and UMAP ------------------------------------------
# This is closer to how Scanpy runs these computations, although not
# exactly -- this should give a "nicer" UMAP and clusters than Seurat's
# default would, i.e., clusters will mostly be together with fewer random
# cells in what look like the wrong places.

message2("Computing neighbor graph and UMAP")

obj <- obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:20,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:20,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

# Cluster cells at several resolutions ------------------------------------

res_tests <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
               0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 4, 6, 8, 10, 12)

for (res in res_tests){
  message2(paste0("Clustering at resolution = ", res))

  obj <- FindClusters(obj,
                      resolution = res,
                      algorithm = 4,
                      graph.name = "RNA_snn",
                      cluster.name = paste0("res", res, "_clusters"),
                      method = "igraph")
}

# DimPlots with clustering labels -----------------------------------------

message2("Making per-resolution DimPlots")

for (res in res_tests){
  p <- DimPlot_scCustom(obj,
                        group.by = paste0("res", res, "_clusters"),
                        label = F,
                        reduction = "harmony_umap")

  ggsave(p,
         filename = paste0(plots_dir, "res", res, ".png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
}

# Score each clustering resolution -----------------------------------------
# The embedding is computed once (on the harmony reduction) and reused
# across all resolutions below, rather than recomputed each time.

message2("Computing embedding for cluster scoring")

embedding <- Embeddings(object = obj[["harmony"]])[, 1:20]
cells <- rownames(embedding)

sil_list <- vector("list", length(res_tests))
purity_list <- vector("list", length(res_tests))
modularity_vec <- numeric(length(res_tests))

for (i in seq_along(res_tests)){
  res <- res_tests[i]
  message2(paste0("Scoring resolution = ", res))

  col <- paste0("res", res, "_clusters")
  clusters <- obj@meta.data[[col]]

  sil <- approx_silhouette(embedding, clusters)
  sil$cell <- cells
  sil$resolution <- res
  sil_list[[i]] <- sil

  purity <- neighbor_purity(obj, clusters)
  purity_list[[i]] <- data.frame(cell = cells, resolution = res, purity = purity)

  modularity_vec[i] <- graph_modularity(obj, clusters)

  # Make two plots that really don't provide unique information from each
  # other. This one is just a jittered dot plot of silhouette scores,
  # where each dot is a cell.

  pal <- DiscretePalette_scCustomize(num_colors = length(unique(sil$cluster)),
                                     palette = "varibow")

  p <- sil %>%
    mutate(cluster = factor(cluster)) %>%
    ggplot(aes(x = cluster,
               y = sil_width)) +
    geom_quasirandom(aes(fill = cluster),
                     shape = 21,
                     color = "black",
                     show.legend = F) +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_fill_manual(values = pal) +
    ggtitle(paste0("Resolution = ", res)) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = 0.5))

  ggsave(p,
         filename = paste0(plots_dir, "cluster_silhouette_res",
                           res, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 6)

  # Alluvial plot showing each cluster's nearest neighbor. This, paired
  # with looking at the silhouette scores and UMAP plots, can help us infer
  # if a cluster is scoring low because the resolution is bad or there are
  # other similar cells nearby that may result in low-confidence
  # clustering.

  p <- sil %>%
    dplyr::select(cell, cluster, neighbor) %>%
    pivot_longer(!cell) %>%
    mutate(value = factor(value)) %>%
    ggplot(aes(x = name,
               stratum = value,
               alluvium = cell,
               fill = value,
               label = value)) +
    geom_flow() +
    geom_stratum(alpha = .7) +
    geom_text(stat = "stratum", size = 3) +
    scale_fill_manual(values = pal) +
    scale_x_discrete(labels = c("Cluster", "Nearest neighbor")) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(color = "black", face = "bold"))
  ggsave(p,
         filename = paste0(plots_dir, "cluster_alluvial_res", res, ".png"),
         units = "in", dpi = 600, bg = "white",
         height = 12, width = 3)
}

message2("Saving cluster-quality tables")

write.csv(bind_rows(sil_list),
          file = paste0(tab_data_out_dir, "approx_silhouette.csv"),
          row.names = F)

write.csv(bind_rows(purity_list),
          file = paste0(tab_data_out_dir, "neighbor_purity.csv"),
          row.names = F)

write.csv(data.frame(resolution = res_tests, modularity = modularity_vec),
          file = paste0(tab_data_out_dir, "graph_modularity.csv"),
          row.names = F)

message2("Saving clustered metadata and UMAP")

# The only output 11_clustering.R (FindAllMarkers()) needs -- expression
# data is unchanged from data/09_integration/<tissue>/bpcells_data, so it
# isn't rewritten here.
saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

# Not required by 11, but saved for convenience/plotting since it isn't
# cheap to recompute.
saveRDS(obj[["harmony_umap"]],
        file = paste0(data_out_dir, "harmony_umap.rds"))
