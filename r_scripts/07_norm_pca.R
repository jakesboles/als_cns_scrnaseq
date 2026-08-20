# Normalizes counts and finds variable features per tissue, ahead of the
# sketch-based dimensional reduction workflow Seurat recommends for datasets
# this size (see 08, not yet written). Splitting layers by orig.ident for
# SketchData() is deferred to 08 -- it's cheap structural bookkeeping, not
# real computation, so there's nothing to gain from persisting it here.

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

data_out_dir <- "data/07_norm_pca/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

message2("Reading in object")

counts <- open_matrix_dir("data/06_obj_reassembly/bpcells")
meta <- readRDS("data/06_obj_reassembly/metadata.rds")
full_obj <- CreateSeuratObject(counts = counts, meta.data = meta)

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

for (i in seq_along(tissues$title)){
  message2(paste0("Processing ", tissues$title[i]))

  obj <- subset(full_obj, 
                subset = tissue == tissues$title[i])

  message2("Normalizing data")

  obj <- NormalizeData(obj)
  
  message2("Finding variable features")
  
  obj <- FindVariableFeatures(obj)
  
  message2("Scaling data")
  
  obj <- ScaleData(obj)
  
  message2("Running PCA")
  
  obj <- RunPCA(obj, npcs = 100)
  
  message2("Running JackStraw")
  
  obj <- JackStraw(obj, num.replicate = 100, dims = 100)
  obj <- ScoreJackStraw(obj, dims = 1:100)
  
  # Saved as a plain data.frame (PC, Score) rather than the full
  # JackStrawData object attached to the reduction, so 09 can evaluate/plot
  # this without needing scale.data or the full object reloaded.
  jackstraw_scores <- JS(obj[["pca"]], slot = "overall.p.values") %>%
    as.data.frame()
  
  message2("Saving sketched data, PCA, variable features, and JackStraw scores")
  
  tissue_out_dir <- paste0(data_out_dir, tissue_file, "/")
  dir.create(tissue_out_dir, showWarnings = F, recursive = T)
  
  # Removed first if present, so a rerun (e.g. after a bug fix or a failed
  # job) doesn't fail on "Path already exists" against a stale/incomplete
  # directory from a previous attempt.
  bpcells_data_dir <- paste0(tissue_out_dir, "bpcells_data")
  if (dir.exists(bpcells_data_dir)){
    unlink(bpcells_data_dir, recursive = T)
  }
  
  write_matrix_dir(mat = obj[["RNA"]]$data,
                   dir = bpcells_data_dir)
  
  saveRDS(obj@meta.data,
          file = paste0(tissue_out_dir, "metadata.rds"))
  
  saveRDS(obj[["pca"]],
          file = paste0(tissue_out_dir, "pca.rds"))
  
  saveRDS(VariableFeatures(obj),
          file = paste0(tissue_out_dir, "variable_features.rds"))
  
  saveRDS(jackstraw_scores,
          file = paste0(tissue_out_dir, "jackstraw_scores.rds"))
  
}
