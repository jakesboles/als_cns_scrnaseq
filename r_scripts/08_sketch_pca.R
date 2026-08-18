# Sketches each tissue's data and runs PCA on the sketch, following Seurat's
# recommended workflow for datasets this size. Splitting layers by
# orig.ident (deferred from 07_norm.R -- see that script's comments) happens
# here, right before SketchData(), since that's the only place it's needed.

# Load libraries
suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_out_dir <- "data/08_sketch_pca/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

message2("Reading in full counts matrix")

# Only opened, not read into memory -- each tissue's cells are subset from
# this lazily below, same pattern as elsewhere in this pipeline.
counts_all <- open_matrix_dir("data/06_obj_reassembly/bpcells")

for (i in seq_along(tissues$file)){
  tissue_file <- tissues$file[i]

  message2(paste0("Building object for ", tissues$title[i]))

  meta <- readRDS(paste0("data/07_norm/", tissue_file, "/metadata.rds"))
  data_mat <- open_matrix_dir(paste0("data/07_norm/", tissue_file, "/bpcells_data"))
  counts_mat <- counts_all[, rownames(meta)]

  obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta)
  obj[["RNA"]]$data <- data_mat

  # Deferred from 07_norm.R: splitting by sample here, right before
  # SketchData(), which needs per-sample layers to sample representative
  # cells across samples rather than across the tissue as a whole.
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

  message2(paste0("Sketching data for ", tissues$title[i]))

  obj <- SketchData(obj,
                    ncells = 2000,
                    method = "LeverageScore",
                    sketched.assay = "sketch")

  DefaultAssay(obj) <- "sketch"

  message2(paste0("Finding variable features in sketched ", tissues$title[i], " data"))

  obj <- FindVariableFeatures(obj)

  message2(paste0("Scaling sketched expression data for ", tissues$title[i]))

  obj <- ScaleData(obj)

  message2(paste0("Running PCA on sketched data for ", tissues$title[i]))

  obj <- RunPCA(obj, npcs = 100)

  message2(paste0("Saving sketched data, PCA, and variable features for ",
                  tissues$title[i]))

  tissue_out_dir <- paste0(data_out_dir, tissue_file, "/")
  dir.create(tissue_out_dir, showWarnings = F, recursive = T)

  write_matrix_dir(mat = obj[["sketch"]]$data,
                   dir = paste0(tissue_out_dir, "bpcells_data"))

  saveRDS(obj@meta.data,
          file = paste0(tissue_out_dir, "metadata.rds"))

  saveRDS(obj[["pca"]],
          file = paste0(tissue_out_dir, "pca.rds"))

  saveRDS(VariableFeatures(obj),
          file = paste0(tissue_out_dir, "variable_features.rds"))
}
