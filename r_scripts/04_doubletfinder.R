# This script runs as a SLURM job array (see jobs/04_doubletfinder.sh), one
# task per sample, instead of looping over all 90 samples in a single job.
# Each task independently loads its own sample's raw per-sample BPCells
# matrix from 01_obj_creation.R and the shared metadata saved by 03_qc2.R,
# so no task ever materializes the whole-cohort object.

# Load libraries
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(DoubletFinder)
  library(scCustomize)
  library(patchwork)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "plots/04_doubletfinder/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/04_doubletfinder/"
dir.create(paste0(data_out_dir, "bpcells_persample/"), showWarnings = F,
           recursive = T)
dir.create(paste0(data_out_dir, "metadata_persample/"), showWarnings = F,
           recursive = T)

# Figure out which sample this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/04_doubletfinder.sh), one task per sample, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

message2("Loading metadata")
meta_all <- readRDS("data/03_qc2/metadata.rds")

# Sorting guarantees the same sample <-> array index mapping every run.
samples <- sort(unique(meta_all$orig.ident))

if (task_id < 1 | task_id > length(samples)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(samples), " samples -- check the --array range in ",
              "jobs/04_doubletfinder.sh against length(unique(orig.ident))."))
}

sample_id <- samples[task_id]

message2(paste0("Processing ", sample_id, " (task ", task_id, "/",
                length(samples), ")"))

# Load and filter this sample's matrix ---------------------------------------

message2("Loading raw per-sample matrix")

# Cell barcodes in the per-sample matrices from 01_obj_creation.R are NOT
# sample-prefixed (that prefixing happens later, at the whole-cohort merge
# step), so it's reapplied here to match 03_qc2.R's metadata, whose cell
# names are "<sample_id>_<barcode>".
mat <- open_matrix_dir(paste0("data/01_obj_creation/bpcells_persample/",
                              sample_id))
colnames(mat) <- paste0(sample_id, "_", colnames(mat))

meta_sample <- meta_all %>%
  filter(orig.ident == sample_id)

# The doublet rate below uses this PRE-QC-filter cell count on purpose:
# doublets form during the original partitioning/loading step, so the rate
# should reflect how many cells were originally captured for this sample,
# not how many survived 03_qc2.R's QC filtering.
n_cells_preqc <- nrow(meta_sample)

meta_sample <- meta_sample %>%
  filter(discard == F) %>%
  column_to_rownames(var = "cell")

mat <- mat[, rownames(meta_sample)]

doublet_rate <- (n_cells_preqc / 10000) * 0.08

obj <- CreateSeuratObject(counts = mat, meta.data = meta_sample)

# Need to remove non-probe genes again
probes <- read.csv("tab_data/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]",
              i = 2) %>% 
  unique()

idx <- rownames(obj) %in% genes
# print(table(idx))
obj <- obj[idx, ]

# Define function to run DoubletFinder ---------------------------------------

run_doubletfinder <- function(s, doublet_rate) {

  # Standard normalization and scaling
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()

  # Default Seurat clustering and UMAP
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunUMAP(s, dims = 1:10)

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))

  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # Run DoubletFinder without homotypic adjustment
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)

  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"

  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  print(pANN)
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)

  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"

  # Plot unadjusted vs. adjusted doublets in UMAP coordinates
  p1 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.unadj",
                      pt.size = 1, shuffle = TRUE, alpha = 0.6) +
       ggtitle("Unadjusted")
  
  p2 <- DimPlot_scCustom(s, reduction = "umap", group.by = "DF.adj",
                      pt.size = 1, shuffle = TRUE, alpha = 0.6) +
       ggtitle("Adjusted for Homotypic Proportion")
  
  p <- p1 + p2 + 
    plot_layout(ncol = 2,
                nrow = 1,
                guides = "collect")

  ggsave(p,
         filename = paste0(plots_dir, unique(s$orig.ident), "_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 10)

  return(s)
}

message2("Running DoubletFinder")

obj <- run_doubletfinder(obj, doublet_rate)

# Save this sample's filtered, DoubletFinder-processed matrix + metadata ----
message2("Saving filtered, DoubletFinder-processed sample")

# Explicit type conversion before writing -- same fix as 01_obj_creation.R's
# write_matrix_dir() segfault (BPCells' compressed writer needs an explicitly
# integer-typed matrix). Cheap insurance in case type tagging doesn't survive
# the round-trip through CreateSeuratObject() and DoubletFinder's own layer
# manipulation.
counts_out <- convert_matrix_type(obj[["RNA"]]$counts, type = "uint32_t")

write_matrix_dir(mat = counts_out,
                 dir = paste0(data_out_dir, "bpcells_persample/", sample_id))

saveRDS(obj@meta.data,
        file = paste0(data_out_dir, "metadata_persample/", sample_id, ".rds"))
