library(Seurat)
library(scCustomize)
library(tidyverse)

b1169 <- "/projects/b1169/boles/als_multitissue_scfrp/"
b1042 <- "/projects/b1042/Gate_Lab/boles/als_multitissue/"
p31535 <- "/projects/p31535/boles/als_multitissue_scfrp/"

plots_dir <- paste0(b1169, "plots/08_sketch_pca/")
dir.create(plots_dir,
           showWarnings = F, 
           recursive = T)

data_out_dir <- paste0(b1169, "data/08_sketch_pca/")
dir.create(data_out_dir, 
           showWarnings = F,
           recursive = T)

message("Reading data")
t0 <- Sys.time()

obj_files <- list.files(paste0(b1169, "data/07_norm/"),
                        full.names = T)

obj_list <- list()
obj_list <- map(obj_files, readRDS)
names(obj_list) <- c("brain", "muscle", "sc")

Sys.time() - t0

for (i in seq_along(obj_list)){
  
  tissue <- names(obj_list)[i]
  
  message(paste0("Sketching data for ", tissue))
  
  t0 <- Sys.time()
  
  obj_list[[i]] <- SketchData(obj_list[[i]],
                              ncells = 2000,
                              method = "LeverageScore",
                              sketched.assay = "sketch")
  
  print(Sys.time() - t0)
  
  DefaultAssay(obj_list[[i]]) <- "sketch"
  
  message(paste0("Finding variable features in sketched ", tissue, " dataset"))
  
  t0 <- Sys.time()
  
  obj_list[[i]] <- FindVariableFeatures(obj_list[[i]])
  
  print(Sys.time() - t0)
  
  message(paste0("Scaling sketched expression data for ", tissue))
  
  t0 <- Sys.time()
  
  obj_list[[i]] <- ScaleData(obj_list[[i]])
  
  print(Sys.time() - t0)
  
  message(paste0("Running PCA on sketched data for ", tissue))
  
  t0 <- Sys.time()
  
  obj_list[[i]] <- RunPCA(obj_list[[i]],
                          npcs = 100)
  
  print(Sys.time() - t0)
  
  message(paste0("Saving ", tissue, " object"))
  
  t0 <- Sys.time()
  
  saveRDS(obj_list[[i]],
          file = paste0(data_out_dir, tissue, "_obj.rds"))
  
  print(Sys.time() - t0)
}