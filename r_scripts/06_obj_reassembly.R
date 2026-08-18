# Reassembles the full 90-sample object from 04_doubletfinder.R's per-sample
# outputs (bpcells_persample/<sample>/ + metadata_persample/<sample>.rds),
# mirroring 01_obj_creation.R's own per-sample-load -> merge -> JoinLayers
# pattern, which is already proven to work end-to-end on this dataset.

# Load libraries
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(scCustomize)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_in_dir <- "data/04_doubletfinder/"

data_out_dir <- "data/06_obj_reassembly/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# Find completed samples ----------------------------------------------------

# Discovered from what 04_doubletfinder.R actually finished writing, rather
# than assumed to be a fixed count of 90, so this fails loudly and
# specifically if any array tasks are missing instead of silently working
# with a partial cohort.
message2("Finding completed samples")

meta_files <- list.files(paste0(data_in_dir, "metadata_persample/"),
                         pattern = "\\.rds$")
samples <- sort(str_remove(meta_files, "\\.rds$"))

expected <- sort(unique(readRDS("data/03_qc2/metadata.rds")$orig.ident))
missing <- setdiff(expected, samples)
if (length(missing) > 0){
  stop(paste0(length(missing), " sample(s) missing from ", data_in_dir,
              "metadata_persample/ -- 04_doubletfinder.R may not have ",
              "finished (or failed) for: ", paste(missing, collapse = ", ")))
}

message2(paste0("Found all ", length(samples), " expected samples"))

# Load and merge per-sample objects ------------------------------------------

message2("Loading per-sample objects")

obj_list <- list()

for (i in seq_along(samples)){
  message(paste0("Loading ", samples[i]))

  mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_persample/", samples[i]))
  meta <- readRDS(paste0(data_in_dir, "metadata_persample/", samples[i], ".rds"))

  obj_list[[i]] <- CreateSeuratObject(counts = mat, meta.data = meta)
}

message2("Merging Seurat objects")

# No add.cell.ids here, unlike 01_obj_creation.R's merge -- cell barcodes in
# 04_doubletfinder.R's per-sample output are already sample-prefixed (04
# applies that prefix itself before saving), unlike 01's raw per-sample
# matrices, which aren't. Prefixing again here would double it.
obj <- Merge_Seurat_List(obj_list)
obj

message2("Joining layers")

obj <- JoinLayers(obj)

# Sanity check metadata rather than assuming it's clean or blindly patching
# specific ids -- 02_qc1.R/03_qc2.R's metadata handling has already been
# fixed once (see PRs #14-#16), so any legacy patch for specific ids here
# would either be redundant or, worse, wrong if circumstances changed. If
# this prints anything, check what's actually happening with the object
# loaded before deciding whether/how to patch it.
message2("Checking metadata for unexpected missing values")

na_cols <- names(which(sapply(obj@meta.data, anyNA)))
if (length(na_cols) > 0){
  message(paste0("Columns with NA values: ", paste(na_cols, collapse = ", ")))
} else {
  message("No NA values found in metadata.")
}

# Save --------------------------------------------------------------------

message2("Saving counts matrix as BPCells on-disk matrix")

counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")

write_matrix_dir(mat = counts_out,
                 dir = paste0(data_out_dir, "bpcells"))

message2("Saving metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts <- open_matrix_dir(paste0(data_out_dir, "bpcells"))
#   meta <- readRDS(paste0(data_out_dir, "metadata.rds"))
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta)
