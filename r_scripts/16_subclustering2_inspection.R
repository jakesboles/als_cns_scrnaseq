suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(BPCells)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")


# Make object manifest ----------------------------------------------------

subclustering_targets <- list(
  list(tissue = "brain", cell_types = "EN"),
  list(tissue = "brain", cell_types = "IN"),
  list(tissue = "brain", cell_types = "Myeloid"),
  list(tissue = "brain", cell_types = "Lymphocyte"),
  list(tissue = "brain", cell_types = "Vascular"),
  list(tissue = "muscle", cell_types = "Muscle fiber"),
  list(tissue = "muscle", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell")),
  list(tissue = "muscle", cell_types = c("Macrophage", "Monocyte", "Neutrophil")),
  list(tissue = "muscle", cell_types = c("T-cell", "NK cell")),
  list(tissue = "sc", cell_types = "Myeloid"),
  list(tissue = "sc", cell_types = "Lymphocyte"),
  list(tissue = "sc", cell_types = c("IN", "EN", "SN", "MN")),
  list(tissue = "sc", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell"))
)

# Assemble object, load marker genes --------------------------------------

i <- 1

tissue <- subclustering_targets[[i]]$tissue
cell <- paste(subclustering_targets[[i]]$cell_types, collapse = "_")

data_in_dir <- paste0("data/15_subclustering2/", tissue, "/", cell, "/")
plots_dir <- paste0("plots/15_subclustering2/", tissue, "/", cell, "/")
tab_dir <- paste0("tab_data/15_subclustering2/", tissue, "/", cell, "/")

data_mat <- open_matrix_dir(paste0(data_in_dir, "/bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
umap <- readRDS(paste0(data_in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony_umap"]] <- umap

markers <- read.csv(paste0(tab_dir, "markers.csv"))

# Get clustering resolution -----------------------------------------------

graph_modularity <- read.csv(paste0(tab_dir, "graph_modularity.csv"))

res <- graph_modularity %>% 
  arrange(desc(modularity)) %>% 
  pull(resolution) %>% 
  .[1]

Idents(obj) <- paste0("res", res, "_clusters")

# Plotting ----------------------------------------------------------------

DimPlot_scCustom(obj)
