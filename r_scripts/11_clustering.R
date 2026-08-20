options(reticulate.verbose = TRUE)
library(Seurat)
library(tidyverse)
library(scCustomize)
library(cluster)
library(ggbeeswarm)
library(factoextra)
library(ggalluvial)
library(stringr)
library(reticulate)

py_module_available(module = "leidenalg")
py_module_available(module = "umap")

# Directories, data -------------------------------------------------------

b1169 <- "/projects/b1169/boles/als_multitissue_scfrp/"

data_in_dir <- paste0(b1169, "data/10_integration")
files <- list.files(data_in_dir,
                    full.names = T)[str_detect(list.files(data_in_dir), 
                                               "cca")]

plots_dir <- paste0(b1169, "plots/11_clustering/")
dir.create(plots_dir, 
           showWarnings = F,
           recursive = T)

data_out_dir <- paste0(b1169, "data/11_clustering/")
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

obj_list <- list()
length(obj_list) <- 3

obj_list <- map(.x = files, .f = readRDS)

obj_list <- map(.x = obj_list, .f = JoinLayers)

# Function to compute NN graphs and run the UMAP on them ------------------
# This is closer to how Scanpy runs these computations, although not exactly
# This should give a "nicer" UMAP and clusters than Seurat's default would
# i.e., clusters will mostly be together with fewer random cells in what look like 
# the wrong places

process_seurat <- function(s){

  s <- s %>%
    FindNeighbors(reduction = "cca_pca",
                  dims = 1:20,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  return.neighbor = T) %>%
    FindNeighbors(reduction = "cca_pca",
                  dims = 1:20,
                  k.param = 15,
                  nn.method = "annoy",
                  annoy.metric = "euclidean",
                  compute.SNN = T) %>%
    RunUMAP(umap.method = "uwot",
            # reduction = "cca_pca",
            nn.name = "sketch.nn",
            metric = "euclidean",
            min.dist = 0.5,
            n_neighbors = 15L,
            reduction.name = "cca_umap",
            return.model = T)
  
  return(s)
}

obj_list <- map(obj_list, process_seurat)

# Cluster cells at several resolutions ------------------------------------

res_tests <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
               0.9, 1, 1.2, 1.4, 1.6, 1.8, 2, 4, 6, 8, 10, 12)

find_clusters <- function(obj){
  obj <- FindClusters(obj, 
                      resolution = res, 
                      algorithm = 4, 
                      graph.name = "sketch_snn",
                      cluster.name = paste0("res", res, "_clusters"),
                      method = "igraph")
}

for (i in seq_along(res_tests)){

  res <- res_tests[i]
  
  message(paste0("Resolution = ", res))
  
  obj_list <- map(obj_list, find_clusters)
}

# DimPlots with clustering labels -----------------------------------------

names(obj_list) <- c("brain", "muscle", "sc")

for (i in seq_along(res_tests)){
  res <- res_tests[i]
  
  for (j in seq_along(obj_list)){
    p <- DimPlot_scCustom(obj_list[[j]],
                          group.by = paste0("res", res, "_clusters"),
                          label = F,
                          reduction = "cca_umap")
    
    ggsave(p,
           filename = paste0(plots_dir, names(obj_list)[j], "_res", res, ".png"),
           units = "in", dpi = 600,
           height = 5, width = 6)
  }
}

# Compute silhouette scores for each clustering resolution ----------------

dist_mat_list <- list()
length(dist_mat_list) <- 3
names(dist_mat_list) <- names(obj_list)
for (i in seq_along(obj_list)){
  dist_mat_list[[i]] <- dist(x = Embeddings(object = obj_list[[i]][["cca_pca"]])[, 1:20])
}

sil_list <- list()
length(sil_list) <- 3
names(sil_list) <- names(obj_list)
for (i in seq_along(sil_list)){
  sil_list[[i]] <- data.frame(cell = attr(dist_mat_list[[i]], "Labels"))
}

for (i in seq_along(obj_list)){
  for (j in seq_along(res_tests)){
    col <- paste0("res", res_tests[j], "_clusters")

    clusters <- obj_list[[i]]@meta.data %>%
      dplyr::select(col) %>%
      na.exclude()

    sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters[, 1])),
                               dist = dist_mat_list[[i]])

    sil[, 1] <- sil[, 1] - 1
    sil[, 2] <- sil[, 2] - 1

    # Make two silhouette plots that really don't provide unique information from each other
    # This one is just a jittered dot plot of silhouette scores, where each dot is a cell

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
      ggtitle(paste0("Resolution = ", res_tests[j])) +
      theme_bw() +
      theme(axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5))

    ggsave(p,
           filename = paste0(plots_dir, names(obj_list)[i], "_cluster_silhouette_res",
                             res_tests[j], ".png"),
           units = "in", dpi = 600,
           height = 4, width = 6)

    # This is a silhouette widths plot from the Factoextra R package

    p <- fviz_silhouette(sil,
                    label = F) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal)

    ggsave(p,
           filename = paste0(plots_dir, names(obj_list)[i], "_cluster_sil_width_res",
                             res_tests[j], ".png"),
           units = "in", dpi = 600,
           height = 4, width = 6)

    # Alluvial plot showing each cluster's nearest neighbor
    # This, paired with looking at the silhouette scores and UMAP plots, can help us
    # infer if a cluster is scoring low because the resolution is bad or there are other
    # similar cells nearby that may result in low-confidence clustering
    sil <- as.data.frame(sil)

    p <- sil %>%
      mutate(cell = sil_list[[i]]$cell) %>%
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
           filename = paste0(plots_dir, names(obj_list)[i],
                            "_cluster_alluvial_res", res_tests[j], ".png"),
           units = "in", dpi = 600, bg = "white",
           height = 12, width = 3)

    # Bind to the large data frame with other clusters' scores

    colnames(sil) <- paste0(colnames(sil), res_tests[j])

    sil_list[[i]] <- cbind(sil_list[[i]], sil)
  }
}

message("Saving brain")
saveRDS(obj_list[[1]],
        file = paste0(data_out_dir, "brain.rds"))

message("Saving muscle")
saveRDS(obj_list[[2]],
        file = paste0(data_out_dir, "muscle.rds"))

message("Saving spinal cord")
saveRDS(obj_list[[3]],
        file = paste0(data_out_dir, "sc.rds"))