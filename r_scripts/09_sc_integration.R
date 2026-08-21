# Integrates the cervical spinal cord data across samples using Harmony, on
# top of the PCA computed in 07_norm_pca.R. Writes the integrated
# ("harmony") reduction alongside the normalized data as an on-disk BPCells
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

data_in_dir <- "data/07_norm_pca/sc/"

plots_dir <- "plots/09_integration/sc/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/09_integration/sc/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# Load the normalized object from 07 -----------------------------------------

message2("Reading in normalized object")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
pca <- readRDS(paste0(data_in_dir, "pca.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["pca"]] <- pca

# IntegrateLayers() defaults to using VariableFeatures(obj), which is empty
# on a freshly-built object -- without this, integration silently runs
# against every gene instead of the already-selected variable ones.
VariableFeatures(obj) <- readRDS(paste0(data_in_dir, "variable_features.rds"))

# Split by sample right before IntegrateLayers(), which needs per-sample
# layers to integrate across -- 07_norm_pca.R saves this data already
# joined into one matrix, so it isn't split until it's actually needed here.
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

message2("Integrating samples using Harmony")

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:20)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join the
# assay's existing layers, so this is still split by sample at this point.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

message2("Computing UMAP for integration diagnostics")

# For sanity-checking the integration below only -- 11_clustering.R computes
# its own UMAP with different parameters as part of its Leiden clustering
# workflow, so this one isn't saved.
obj <- RunUMAP(obj,
               dims = 1:20,
               reduction = "harmony",
               reduction.name = "harmony_umap",
               reduction.key = "harmonyumap_")

message2("Making integration diagnostic plots")

p <- DimPlot_scCustom(obj,
                      reduction = "harmony_umap",
                      group.by = "tissue",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "harmony_umap_dimplot_tissue.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "harmony_umap",
                      group.by = "batch",
                      colors_use = DiscretePalette_scCustomize(6,
                                                               palette = "ditto_seq"))
ggsave(p,
       filename = paste0(plots_dir, "harmony_umap_dimplot_batch.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "harmony_umap",
                      group.by = "group",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "harmony_umap_dimplot_group.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      group.by = "id",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, "harmony_umap_dimplot_id.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

message2("Saving Harmony-integrated data")

# Removed first if present, so a rerun (e.g. while tuning integration
# parameters) doesn't fail on "Path already exists" against a stale
# directory from a previous attempt.
bpcells_data_dir <- paste0(data_out_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

saveRDS(obj[["harmony"]],
        file = paste0(data_out_dir, "harmony.rds"))
