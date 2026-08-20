# Computes a neighbor graph/UMAP and Leiden-clusters one tissue's
# Harmony-integrated data across a range of resolutions, scoring each with
# silhouette width, for the user to pick a resolution to carry forward.
# Runs as a SLURM job array (see jobs/10_clustering.sh), one task per
# tissue, since all 3 tissues are fully independent (same restructuring as
# 04_doubletfinder.R, 07_norm_pca.R, and 08_sketch_pca.R). 11_clustering.R
# (FindAllMarkers()) only needs cluster labels plus the expression data
# already saved by 09_*_integration.R, which this script doesn't touch --
# so metadata.rds (with one "resX_clusters" column per tested resolution)
# is the only required output. The Harmony UMAP is also saved for
# convenience/plotting, but isn't needed by 11.

# Load libraries
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(cluster)
  library(ggbeeswarm)
  library(factoextra)
  library(ggalluvial)
  library(reticulate)
  library(BPCells)
})

options(reticulate.verbose = TRUE)
py_module_available(module = "leidenalg")

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
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

plots_dir <- "plots/10_clustering/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- paste0("data/10_clustering/", tissue_file, "/")
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# Load the Harmony-integrated object from 09 ---------------------------------

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
         filename = paste0(plots_dir, tissue_file, "_res", res, ".png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
}

# Compute silhouette scores for each clustering resolution ----------------
# The distance matrix is computed once (on the Harmony embedding) and
# reused across all resolutions below, rather than recomputed each time --
# note this is an O(cells^2) matrix, which may be memory-heavy for the
# larger tissues.

message2("Computing distance matrix for silhouette scoring")

dist_mat <- dist(x = Embeddings(object = obj[["harmony"]])[, 1:20])

sil_df <- data.frame(cell = attr(dist_mat, "Labels"))

for (res in res_tests){
  message2(paste0("Scoring silhouette at resolution = ", res))

  col <- paste0("res", res, "_clusters")

  clusters <- obj@meta.data %>%
    dplyr::select(all_of(col)) %>%
    na.exclude()

  sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters[, 1])),
                             dist = dist_mat)

  sil[, 1] <- sil[, 1] - 1
  sil[, 2] <- sil[, 2] - 1

  # Make two silhouette plots that really don't provide unique information
  # from each other. This one is just a jittered dot plot of silhouette
  # scores, where each dot is a cell.

  pal <- DiscretePalette_scCustomize(num_colors = length(unique(clusters[, 1])),
                                     palette = "varibow")

  p <- sil %>%
    as.data.frame() %>%
    mutate(cluster = cluster - 1) %>%
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
         filename = paste0(plots_dir, tissue_file, "_cluster_silhouette_res",
                           res, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 6)

  # This is a silhouette widths plot from the Factoextra R package

  p <- fviz_silhouette(sil,
                  label = F) +
    scale_color_manual(values = pal) +
    scale_fill_manual(values = pal)

  ggsave(p,
         filename = paste0(plots_dir, tissue_file, "_cluster_sil_width_res",
                           res, ".png"),
         units = "in", dpi = 600,
         height = 4, width = 6)

  # Alluvial plot showing each cluster's nearest neighbor. This, paired
  # with looking at the silhouette scores and UMAP plots, can help us infer
  # if a cluster is scoring low because the resolution is bad or there are
  # other similar cells nearby that may result in low-confidence
  # clustering.
  sil <- as.data.frame(sil)

  p <- sil %>%
    mutate(cell = sil_df$cell) %>%
    pivot_longer(!cell) %>%
    filter(name != "sil_width") %>%
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
         filename = paste0(plots_dir, tissue_file,
                          "_cluster_alluvial_res", res, ".png"),
         units = "in", dpi = 600, bg = "white",
         height = 12, width = 3)
}

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
