# Normalizes counts, finds variable features, scales, runs PCA, and runs
# JackStraw for one tissue, on the full (non-sketched) data. Runs as a
# SLURM job array (see jobs/07_norm_pca.sh), one task per tissue, since all
# 3 tissues are fully independent and there's no reason to process them one
# at a time in a single job (same restructuring as 04_doubletfinder.R and
# 08_sketch_pca.R). Saves the normalized data, PCA, variable features, and
# JackStraw scores as BPCells/RDS output for the next steps: appraising the
# PCA to pick integration parameters, then integration itself.

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

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/07_norm_pca.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/07_norm_pca.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

# Build the tissue's object --------------------------------------------------

message2("Reading in object")

counts <- open_matrix_dir("data/06_obj_reassembly/bpcells")
meta <- readRDS("data/06_obj_reassembly/metadata.rds")
full_obj <- CreateSeuratObject(counts = counts, meta.data = meta)

obj <- subset(full_obj,
              subset = tissue == tissue_title)

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

message2("Saving normalized data, PCA, variable features, and JackStraw scores")

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
