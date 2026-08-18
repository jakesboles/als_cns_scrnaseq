suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(stringr)
  library(glmGamPoi)
  library(BPCells)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

b1169_dir <- "/projects/b1169/boles/als_multitissue_scfrp/"
b1042_dir <- "/projects/b1042/Gate_Lab/boles/als_multitissue/"

plots_dir <- paste0(b1169_dir, "plots/07_norm/")
dir.create(plots_dir,
           showWarnings = F, recursive = T)

data_dir <- paste0(b1169_dir, "data/07_norm/")
dir.create(data_dir,
           showWarnings = F, recursive = T)

message("Reading in object")
t0 <- Sys.time()

big_obj <- readRDS(paste0(b1169_dir, "data/06_obj_assembly/obj.rds"))

t1 <- Sys.time()
t1 - t0

tissues <- c("Motor cortex", "Skeletal muscle", "Cervical spinal cord")
files <- c("brain", "muscle", "sc")

for (i in seq_along(tissues)){
  obj <- subset(big_obj, subset = tissue == tissues[i])
  
  message("Normalizing data")
  t2 <- Sys.time()
  
  obj <- NormalizeData(obj)
  
  t3 <- Sys.time()
  t3 - t2
  
  message("Splitting layers")
  
  obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
  
  t4 <- Sys.time()
  t4 - t3
  
  message("Finding variable features")
  
  obj <- FindVariableFeatures(obj, 
                              verbose = T)
  
  t5 <- Sys.time()
  t5 - t4
  
  message("Saving data")
  
  saveRDS(obj,
          file = paste0(data_dir, "lognorm_", files[i], "_obj.rds"))
  
  t6 <- Sys.time()
  t6 - t5
}
