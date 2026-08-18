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

plots_dir <- "plots/07_norm/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/07_norm/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

message2("Reading in object")

counts <- open_matrix_dir("data/06_obj_reassembly/bpcells")
meta <- readRDS("data/06_obj_reassembly/metadata.rds")
obj <- CreateSeuratObject(counts = counts, meta.data = meta)

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

for (i in seq_along(tissues$title)){
  message2(paste0("Processing ", tissues$title[i]))

  tissue_obj <- subset(obj, subset = tissue == tissues$title[i])

  message2("Normalizing data")

  tissue_obj <- NormalizeData(tissue_obj)

  message2("Finding variable features")

  tissue_obj <- FindVariableFeatures(tissue_obj, verbose = T)

  tissue_out_dir <- paste0(data_out_dir, tissues$file[i], "/")
  dir.create(tissue_out_dir, showWarnings = F, recursive = T)

  message2("Saving normalized data as BPCells on-disk matrix")

  # No convert_matrix_type() here, unlike the counts-matrix saves elsewhere
  # in this pipeline -- log-normalized values are genuinely continuous, so
  # forcing an integer type would be wrong, not just a missed optimization.
  write_matrix_dir(mat = tissue_obj[["RNA"]]$data,
                   dir = paste0(tissue_out_dir, "bpcells_data"))

  message2("Saving metadata and variable features")

  # Raw counts aren't re-saved per tissue -- 08 should pull this tissue's
  # cells from data/06_obj_reassembly/bpcells directly (matching this
  # metadata's rownames) rather than duplicating the counts matrix here.
  saveRDS(tissue_obj@meta.data,
          file = paste0(tissue_out_dir, "metadata.rds"))

  saveRDS(VariableFeatures(tissue_obj),
          file = paste0(tissue_out_dir, "variable_features.rds"))
}
