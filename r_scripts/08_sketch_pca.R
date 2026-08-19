# Sketches one tissue's data and runs PCA on the sketch, following Seurat's
# recommended workflow for datasets this size. Runs as a SLURM job array
# (see jobs/08_sketch_pca.sh), one task per tissue, since all 3 tissues are
# fully independent and there's no reason to process them one at a time in
# a single job. Splitting layers by orig.ident (deferred from 07_norm.R --
# see that script's comments) happens here, right before SketchData(),
# since that's the only place it's needed. JackStraw() also runs here
# rather than in 09 -- it needs the object's scale.data to repeat PCA on
# permuted data, which 09 doesn't reload.

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

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/08_sketch_pca.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/08_sketch_pca.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

# Build the tissue's object --------------------------------------------------

message2("Reading in full counts matrix")

# Only opened, not read into memory -- this tissue's cells are subset from
# this lazily below, same pattern as elsewhere in this pipeline.
counts_all <- open_matrix_dir("data/06_obj_reassembly/bpcells")

meta <- readRDS(paste0("data/07_norm/", tissue_file, "/metadata.rds"))
data_mat <- open_matrix_dir(paste0("data/07_norm/", tissue_file, "/bpcells_data"))
counts_mat <- counts_all[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts_mat, meta.data = meta)
obj[["RNA"]]$data <- data_mat

# Deferred from 07_norm.R: splitting by sample here, right before
# SketchData(), which needs per-sample layers to sample representative
# cells across samples rather than across the tissue as a whole.
obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

# SketchData()'s LeverageScore method needs variable features already set
# per (per-sample) layer -- without them it falls back to using every
# gene, forcing a dense conversion of the full sparse matrix and blowing
# through LeverageScore's internal time budget ("too slow"). This is
# separate from the FindVariableFeatures() call below on the sketch
# subset, which serves a different purpose (feeding ScaleData/RunPCA).
obj <- FindVariableFeatures(obj)

message2(paste0("Sketching data for ", tissue_title))

obj <- SketchData(obj,
                  ncells = 2000,
                  method = "LeverageScore",
                  sketched.assay = "sketch")

DefaultAssay(obj) <- "sketch"

message2(paste0("Finding variable features in sketched ", tissue_title, " data"))

obj <- FindVariableFeatures(obj)

message2(paste0("Scaling sketched expression data for ", tissue_title))

obj <- ScaleData(obj)

message2(paste0("Running PCA on sketched data for ", tissue_title))

obj <- RunPCA(obj, npcs = 100)

message2(paste0("Running JackStraw for ", tissue_title))

obj <- JackStraw(obj, num.replicate = 100, dims = 100)
obj <- ScoreJackStraw(obj, dims = 1:100)

# Saved as a plain data.frame (PC, Score) rather than the full
# JackStrawData object attached to the reduction, so 09 can evaluate/plot
# this without needing scale.data or the full object reloaded.
jackstraw_scores <- JS(obj[["pca"]], slot = "overall.p.values") %>%
  as.data.frame()

message2(paste0("Saving sketched data, PCA, variable features, and ",
                "JackStraw scores for ", tissue_title))

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

saveRDS(jackstraw_scores,
        file = paste0(tissue_out_dir, "jackstraw_scores.rds"))
