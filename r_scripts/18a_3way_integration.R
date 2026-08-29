options(future.globals.maxSize = 264 * 1024^3)
library(Seurat)
library(tidyverse)
library(scCustomize)

setwd("/projects/b1169/boles/als_multitissue_scfrp/")

plots_dir <- "plots/22_full_integration/"
dir.create(plots_dir)

data_dir <- "data/22_full_integration/"
dir.create(data_dir)

# files <- list.files("data/19_more_cleanup/",
#                     full.names = T)
# files <- files[c(1, 5, 6)]
# 
# objs <- map(files,
#             readRDS)
# 
# for (i in seq_along(objs)){
#   
#   message(paste0("Data set ", i, " out of 3"))
#   
#   objs[[i]][["sketch"]] <- NULL
#   
#   objs[[i]]@reductions <- list()
#   objs[[i]]@neighbors <- list() 
#   objs[[i]]@graphs <- list() 
#   
#   if (i == 1) {
#     new_meta <- read.csv(paste0("tab_data/19_more_cleanup/ex_neurons/brain_final_annotations.csv"))
#     
#     new_meta <- new_meta %>%
#       dplyr::select(c(bc, y)) %>%
#       dplyr::rename("final_label2" = "y") %>%
#       mutate(bc = if_else(bc %in% colnames(objs[[i]]), 
#                           bc,
#                           str_remove(bc, "_"))) %>%
#       column_to_rownames(var = "bc")
#     
#     objs[[i]]@meta.data <- objs[[i]]@meta.data %>%
#       dplyr::select(-c(final_label2))
#     
#     objs[[i]] <- AddMetaData(objs[[i]],
#                        new_meta)
#     
#     objs[[i]] <- objs[[i]] %>%
#       subset(final_label2 != "remove")
#   }
#   
# }
# 
# obj <- Merge_Seurat_List(objs)
# 
# obj <- NormalizeData(obj) %>%
#   FindVariableFeatures() %>%
#   ScaleData() %>%
#   RunPCA()
# 
# p <- ElbowPlot(obj,
#                ndims = 50)
# ggsave(p,
#        filename = paste0(plots_dir, "nosketch_pca_elbow.png"),
#        units = "in", dpi = 600,
#        bg = "white",
#        height = 6, width = 6)
# 
# Iterate_PC_Loading_Plots(obj,
#                          file_path = plots_dir,
#                          file_name = "nosketch_pca_loadings")
# 
# # Integration parameters
# integration_method <- "HarmonyIntegration"
num_pcs <- 1:15
# 
# DefaultAssay(obj) <- "RNA"
# 
# obj <- IntegrateLayers(obj,
#                        method = integration_method,
#                        orig.reduction = "pca",
#                        new.reduction = "integrated_pca",
#                        k.anchor = 20,
#                        reference = which(Layers(obj, search = "data") %in% 
#                                            c("data.GBB-23-11_b.1.1", "data.GWF19-47_b.1.1",
#                                              "data.GBB-23-11_m.2.1", "data.AU-073_m.2.1",
#                                              "data.GWF21-56_s.2", "data.GBB-18-13_s.2")),
#                        dims = num_pcs)

# if this has been run already and you're trying to tweak the UMAP
obj <- readRDS(paste0(data_dir, "harmonyintegration_nosketch.rds"))

message("Computing NN and UMAP")

obj <- obj %>%
  FindNeighbors(reduction = "integrated_pca",
                dims = num_pcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "integrated_pca",
                dims = num_pcs,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T)
  # RunUMAP(umap.method = "uwot",
  #         # reduction = "cca_pca",
  #         nn.name = "RNA.nn",
  #         metric = "euclidean",
  #         min.dist = 0.5,
  #         n_neighbors = 15L,
  #         reduction.name = "integrated_umap",
  #         return.model = F)

obj <- RunUMAP(obj,
               umap.method = "uwot",
               reduction = "integrated_pca",
               dims = num_pcs,
               # nn.name = "RNA.nn",
               metric = "euclidean",
               min.dist = 0.5,
               n_neighbors = 15L,
               reduction.name = "integrated_umap2",
               return.model = F)

for (group in c("final_label2", "Batch", "tissue", "orig.ident", "Group")){
  
  if (group == "final_label2"){
    file1 <- "celltype"
  } else { 
    file1 <- str_to_lower(group)
  }

  if (group %in% c("final_label2", "orig.ident")) { 
    w = 15 
  } else { 
    w = 11
  }
 
  # p <- DimPlot_scCustom(obj,
  #                       reduction = "integrated_umap",
  #                       group.by = group)
  # ggsave(p,
  #        filename = paste0(plots_dir, file1, "_harmonyintegration_nosketch_umap1.png"),
  #        units = "in", dpi = 600,
  #        height = 8, width = w)
  
  p <- DimPlot_scCustom(obj,
                        reduction = "integrated_umap2",
                        group.by = group)
  ggsave(p,
         filename = paste0(plots_dir, file1, "_harmonyintegration_nosketch_umap2.png"),
         units = "in", dpi = 600,
         height = 8, width = w)
   
}

saveRDS(obj,
        file = paste0(data_dir, "harmonyintegration_nosketch.rds"))

