# Integrates the cervical spinal cord sketch data across samples using CCA,
# on top of the PCA computed in 08_sketch_pca.R. Writes the integrated
# ("cca_pca") reduction alongside the sketch data as an on-disk BPCells
# matrix, for 11_clustering.R to load and Leiden-cluster. Clustering itself
# (and any marker-gene/cell-type inspection) is deliberately left to
# 11/later steps -- this script's job ends at producing a corrected
# embedding, plus a few plots to sanity check that the correction worked.

# Load libraries
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_in_dir <- "data/08_sketch_pca/sc/"

plots_dir <- "plots/10_integration/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/10_integration/sc/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# Load the sketch object from 08 --------------------------------------------

message2("Reading in sketch object")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
pca <- readRDS(paste0(data_in_dir, "pca.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "sketch")
obj[["sketch"]]$data <- data_mat
obj[["pca"]] <- pca

# IntegrateLayers() defaults to using VariableFeatures(obj), which is empty
# on a freshly-built object -- without this, integration silently runs
# against every gene instead of the ~2000 already-selected variable ones,
# which is likely a major reason this has been so slow historically (the
# same failure mode debugged in 08_sketch_pca.R's SketchData() call).
VariableFeatures(obj) <- readRDS(paste0(data_in_dir, "variable_features.rds"))

# Split by sample right before IntegrateLayers(), which needs per-sample
# layers to integrate across -- 08_sketch_pca.R saves this data already
# joined into one matrix, so it isn't split until it's actually needed here.
obj[["sketch"]] <- split(obj[["sketch"]], f = obj$orig.ident)

message2("Integrating samples using CCA")

# reference restricts anchor-finding to these 2 samples vs. every other
# sample, rather than every possible pair (O(n) instead of O(n^2) anchor
# computations) -- carried over from the original script's choice of
# reference samples.
obj <- IntegrateLayers(obj,
                       method = "CCAIntegration",
                       orig.reduction = "pca",
                       new.reduction = "cca_pca",
                       k.anchor = 20,
                       reference = which(Layers(obj, search = "data") %in%
                                           c("data.GWF21-56_s", "data.GBB-18-13_s")),
                       dims = 1:20)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join the
# assay's existing layers, so this is still split by sample at this point.
obj[["sketch"]] <- JoinLayers(obj[["sketch"]])

message2("Computing UMAP for integration diagnostics")

# For sanity-checking the integration below only -- 11_clustering.R computes
# its own UMAP with different parameters as part of its Leiden clustering
# workflow, so this one isn't saved.
obj <- RunUMAP(obj,
               dims = 1:20,
               reduction = "cca_pca",
               reduction.name = "cca_umap",
               reduction.key = "cca_umap_")

message2("Making integration diagnostic plots")

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "tissue",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "sc_cca_umap_dimplot_tissue.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "batch",
                      colors_use = DiscretePalette_scCustomize(6,
                                                               palette = "ditto_seq"))
ggsave(p,
       filename = paste0(plots_dir, "sc_cca_umap_dimplot_batch.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "group",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "sc_cca_umap_dimplot_group.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      group.by = "id",
                      reduction = "cca_umap")
ggsave(p,
       filename = paste0(plots_dir, "sc_cca_umap_dimplot_id.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

message2("Saving CCA-integrated sketch data")

write_matrix_dir(mat = obj[["sketch"]]$data,
                 dir = paste0(data_out_dir, "bpcells_data"))

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

saveRDS(obj[["cca_pca"]],
        file = paste0(data_out_dir, "cca_pca.rds"))
