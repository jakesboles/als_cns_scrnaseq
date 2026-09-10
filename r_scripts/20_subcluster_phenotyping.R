suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(scater)
  library(igraph)
  library(scales)
  library(BPCells)
  library(dittoSeq)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

set.seed(256)

# Assemble object -----------

cell_type <- "microglia"

data_dir <- paste0("data/20_subcluster_phenotyping/", cell_type, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0("results/20_subcluster_phenotyping/", cell_type, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

in_dir <- paste0("data/19_subclustering3/", cell_type)

meta <- readRDS(paste0(in_dir, "/metadata.rds"))
harmony <- readRDS(paste0(in_dir, "/harmony.rds"))
umap <- readRDS(paste0(in_dir, "/harmony_umap.rds"))

mat <- open_matrix_dir(paste0(in_dir, "/bpcells_data"))

obj <- CreateSeuratObject(counts = mat,
                          meta.data = meta)
obj[["RNA"]]$data <- mat
obj[["harmony"]] <- harmony
obj[["umap"]] <- umap

obj$group <- factor(obj$group, levels = c("Control", "sALS", "C9orf72"))

# Find neighbors and cluster ----------------------------------------------

obj <- FindNeighbors(obj,
                     reduction = "harmony",
                     dims = 1:15,
                     k.param = 15,
                     nn.method = "annoy",
                     annoy.metric = "euclidean",
                     compute.SNN = T)

obj <- FindClusters(obj,
                    algorithm = 4,
                    method = "igraph",
                    resolution = 2)

DimPlot_scCustom(obj,
                 label = F) + 

dittoBarPlot(obj,
             var = c("group"), 
             group.by = "seurat_clusters")
