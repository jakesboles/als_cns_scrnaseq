# This script is to assemble the object from CellBender-corrected count matrices for
# all 90 samples from the multi-system ALS project.
# The next scripts will assemble QC metrics and discard low-quality cells,
# but I expect the object will be so large that it will not be feasible to
# tinker with QC thresholds and other analysis parameters at the same time
# as creating the main object.
# With >1.3M cells prior to QC, the counts matrix is saved on-disk with BPCells
# instead of as a single in-memory RDS, so downstream scripts can load it without
# reading the full matrix into memory.

# Load libraries
suppressMessages({
  library(plyr)
  library(tidyverse)
  library(Seurat)
  library(stringr)
  library(scCustomize)
  library(BPCells)
})

# Function to print clear log progress updates
message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Filter operator
`%notin%` <- Negate(`%in%`)

# Create directories ------------------------------------------------------
message2("Creating directories")

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_in_dir <- "/projects/b1042/Gate_Lab/boles/als_multitissue/CellbenderOutput"

data_out_dir <- "data/01_obj_creation/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

meta <- read.csv("tab_data/metadata.csv")

probes <- read.csv("tab_data/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv",
                   skip = 5)
genes <- probes$probe_id %>%
  str_split_i(pattern = "[|]",
              i = 2) %>% 
  unique()

# Create Seurat objects from corrected counts -----------------------------
message2("Creating Seurat objects")

# Each sample's matrix is written to its own on-disk BPCells directory right
# after reading, so only one sample's matrix is fully in memory at a time
# instead of holding all 90 samples' matrices at once through the merge below.
bpcells_persample_dir <- paste0(data_out_dir, "bpcells_persample/")
dir.create(bpcells_persample_dir, showWarnings = F, recursive = T)

create_object <- function(file, id){
  counts <- Read_CellBender_h5_Mat(file)
  bp_dir <- paste0(bpcells_persample_dir, id)
  write_matrix_dir(mat = counts, dir = bp_dir)
  mat <- open_matrix_dir(dir = bp_dir)
  return(CreateSeuratObject(counts = mat,
                            project = id)) # These thresholds are new as of 12/26/2024
}

cellbender_dirs <- list.dirs(data_in_dir,
                             recursive = F)[2:7]
samples <- list()
files <- list()
for (i in seq_along(cellbender_dirs)){
  samples[[i]] <- list.dirs(cellbender_dirs[i],
                          recursive = F)
  samples[[i]] <- samples[[i]][str_detect(samples[[i]], "UWA",
                                          negate = T)]
  
  files[[i]] <- paste0(samples[[i]], "_filtered.h5")
}

list <- list()

samples <- list_c(samples)
files <- list_c(files)

table(file.exists(files))

ids <- str_split(samples, "/") %>%
  lapply(`[`, 9) %>%
  unlist()

# ids %>% str_split_i("_", i = 1) %>% unique()

ids <- str_replace_all(ids, "GWF-19-47", "GWF19-47")
ids <- str_replace_all(ids, "GWF-20-54", "GWF20-54")
ids <- str_replace_all(ids, "GWF-21-56", "GWF21-56")

for (i in seq_along(samples)){
  message(paste0("Creating object for ", ids[i]))
  list[i] <- create_object(files[i], ids[i])
  print(list[i])
}

for (i in seq_along(samples)){
  message(paste0("Removing non-probe list genes from ", ids[i]))
  idx <- rownames(list[[i]]) %in% genes
  print(table(idx))
  list[[i]] <- list[[i]][idx, ]
}

message2("Merging Seurat objects")

obj <- Merge_Seurat_List(list,
                         add.cell.ids = ids)
obj

message2("Joining layers")

obj <- JoinLayers(obj)

message2("Saving counts matrix as BPCells on-disk matrix")

write_matrix_dir(mat = obj[["RNA"]]$counts,
                 dir = paste0(data_out_dir, "bpcells"))

message2("Saving metadata as RDS")

saveRDS(obj@meta.data,
        file = paste0(data_out_dir,
                      "metadata.rds"))

# Downstream scripts should reconstruct the object from these on-disk pieces:
#   counts <- open_matrix_dir(paste0(data_out_dir, "bpcells"))
#   meta <- readRDS(paste0(data_out_dir, "metadata.rds"))
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta)
