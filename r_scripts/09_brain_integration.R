# Integrates the brain (motor cortex) data across samples using CCA, on
# top of the PCA computed in 07_norm_pca.R. Writes the integrated
# ("cca") reduction alongside the normalized data as an on-disk BPCells
# matrix, for 10_clustering.R to load and Leiden-cluster. Clustering itself
# (and any marker-gene/cell-type inspection) is deliberately left to
# 10/later steps -- this script's job ends at producing a corrected
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

data_in_dir <- "data/07_norm_pca/brain/"

plots_dir <- "plots/09_integration/brain/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/09_integration/brain/"
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

message2("Integrating samples using CCA")

# reference restricts anchor-finding to these 2 samples vs. every other
# sample, rather than every possible pair (O(n) instead of O(n^2) anchor
# computations) -- carried over from the earlier (sketch-data) version of
# this script's choice of reference samples.
#
# reference must be an integer index into the per-sample "data" layers, in
# the order IntegrateLayers() sees them. The original approach here,
# `which(Layers(obj, search = "data") %in% c(...))`, is what crashed the
# earlier sketch-data CCA attempt: on this Seurat/SeuratObject version,
# Layers(obj, search = "data") warns "No layers found matching search
# pattern provided" and returns nothing, so `reference` silently becomes
# integer(0) instead of the intended 2 indices -- IntegrateLayers() then
# mishandles that empty reference internally (repeated "only the first
# layer is used" warnings, then a crash trying to subset barcodes that
# don't exist under a mangled "_<sample>_<barcode>" name).
#
# split() (above) sorts its resulting per-sample layers alphabetically by
# orig.ident (orig.ident is a plain character column throughout this
# pipeline, never a factor with a custom level order -- see 03_qc2.R/
# 04_doubletfinder.R using the same sort(unique(...)) assumption), so the
# reference indices can be computed directly from sample IDs instead,
# without going through Layers() at all.
#
# CCA needs dense access to expression data (unlike Harmony, which only
# needs the PCA embedding), so expect a "on-disk CCA Integration is not
# currently supported" warning -- Seurat falls back to converting the
# BPCells matrix to an in-memory dgCMatrix internally to run this.
obj <- IntegrateLayers(obj,
                       method = "CCAIntegration",
                       orig.reduction = "pca",
                       new.reduction = "cca",
                       k.anchor = 20,
                       reference = which(sort(unique(obj$orig.ident)) %in%
                                           c("GBB-23-11_b", "GWF19-47_b")),
                       dims = 1:20)

# IntegrateLayers() only adds the new reduction -- it doesn't touch/join the
# assay's existing layers, so this is still split by sample at this point.
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

message2("Computing UMAP for integration diagnostics")

# For sanity-checking the integration below only -- 10_clustering.R computes
# its own UMAP with different parameters as part of its Leiden clustering
# workflow, so this one isn't saved.
obj <- RunUMAP(obj,
               dims = 1:20,
               reduction = "cca",
               reduction.name = "cca_umap",
               reduction.key = "ccaumap_")

message2("Making integration diagnostic plots")

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "tissue",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "cca_umap_dimplot_tissue.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "batch",
                      colors_use = DiscretePalette_scCustomize(6,
                                                               palette = "ditto_seq"))
ggsave(p,
       filename = paste0(plots_dir, "cca_umap_dimplot_batch.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "group",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, "cca_umap_dimplot_group.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      group.by = "id",
                      reduction = "cca_umap")
ggsave(p,
       filename = paste0(plots_dir, "cca_umap_dimplot_id.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

message2("Saving CCA-integrated data")

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

saveRDS(obj[["cca"]],
        file = paste0(data_out_dir, "cca.rds"))
