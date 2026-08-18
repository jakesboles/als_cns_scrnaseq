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

# message("Making elbow plot")
# 
# p <- ElbowPlot(obj)
# ggsave(p,
#        filename = paste0(plots_dir, "elbow_plot.png"),
#        units = "in", dpi = 600,
#        height = 6, width = 6)
# 
# message("Making PC loading plots")
# 
# Iterate_PC_Loading_Plots(obj,
#                          file_path = plots_dir,
#                          file_name = "pca_loadings.pdf")
# 
# pca_plot2 <- function(dims){
#   DimPlot_scCustom(obj,
#                    reduction = "pca",
#                    label = F,
#                    dims = dims) +
#     # guides(color = guide_legend(position = "inside")) +
#     theme(axis.text = element_text(size = 8),
#           axis.title = element_text(size = 10),
#           legend.justification = "left")
# }
# 
# pca_grid <- function(group.by){
# 
#   Idents(obj) <- group.by
# 
#   design <- "
#   AA######
#   BBCC####
#   DDEEFF##
#   GGHHIIJJ
#   "
# 
#   pca_plot2(c(1,2)) +
#     pca_plot2(c(1,3)) + pca_plot2(c(2,3)) +
#     pca_plot2(c(1,4)) + pca_plot2(c(2,4)) + pca_plot2(c(3,4)) +
#     pca_plot2(c(1,5)) + pca_plot2(c(2,5)) + pca_plot2(c(3,5)) + pca_plot(c(4,5)) +
#     plot_layout(design = design)
# }
# 
# message("Making gridded PCA DimPlots, colored by variables of interest")
# 
# p <- pca_grid("tissue")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_tissue.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 10)
# 
# p <- pca_grid("orig.ident")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_sample.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 10)
# 
# p <- pca_grid("id")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_id.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 10)
# 
# p <- pca_grid("Batch")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_batch.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 10)
# 
# p <- pca_grid("Group")
# ggsave(p,
#        filename = paste0(plots_dir, "pca_dimplot_group.png"),
#        units = "in", dpi = 600,
#        height = 10, width = 10)