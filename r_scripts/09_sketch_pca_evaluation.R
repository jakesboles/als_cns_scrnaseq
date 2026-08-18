library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- paste0(b1169, "plots/09_sketch_pca_evaluation/")
dir.create(plots_dir,
           showWarnings = F, 
           recursive = T)

data_in_dir <- paste0(b1169, "data/08_sketch_pca/")
dir.create(data_out_dir, 
           showWarnings = F,
           recursive = T)

files <- list.files(data_in_dir, 
                    full.names = T)
obj_list <- map(files, 
                readRDS)
names(obj_list) <- c("brain", "muscle", "sc")

for (i in seq_along(obj_list)){
  print(names(obj_list)[i])
  print(Assays(obj_list[[i]]))
  print(DefaultAssay(obj_list[[i]]))
  
  obj_list[[i]]@meta.data <- obj_list[[i]]@meta.data %>%
    mutate(site = case_when(str_detect(id, "AU") ~ "WashU",
                            str_detect(id, "GWF") ~ "Barrow",
                            str_detect(id, "GBB") ~ "Georgetown"))
}

for (i in seq_along(obj_list)){
  
  obj <- obj_list[[i]]
  
  tissue <- names(obj_list)[i]
  
  message("Making elbow plot")
  
  p <- ElbowPlot(obj,
                 ndims = 100)
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_elbow_plot.png"),
         units = "in", dpi = 600,
         height = 6, width = 6,
         bg = "white")
  
  # pick the first 30 PCs
  
  message("Making PC loading plots")
  
  Iterate_PC_Loading_Plots(obj,
                           file_path = plots_dir,
                           file_name = paste0(tissue, "_pca_loadings"))
  
  pca_plot2 <- function(dims){
    
    if (length(unique(Idents(obj))) == 3) {
      pal <- JCO_Four()
    } else {
      if (length(unique(Idents(obj))) == 90) {
        pal <- DiscretePalette_scCustomize(90,
                                           palette = "varibow",
                                           shuffle_pal = T)
      } else {
        pal <- DiscretePalette_scCustomize(num_colors = length(unique(Idents(obj))),
                                           palette = "glasbey")
      }
    }
    
    DimPlot_scCustom(obj,
                     reduction = "pca",
                     label = F,
                     dims = dims,
                     colors_use = pal) +
      # guides(color = guide_legend(position = "inside")) +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.justification = "left")
  }
  
  pca_grid <- function(group.by){
    
    # Idents(obj) <- group.by
    
    design <- "
  AA######
  BBCC####
  DDEEFF##
  GGHHIIJJ
  "
    
    pca_plot2(c(1,2)) +
      pca_plot2(c(1,3)) + pca_plot2(c(2,3)) +
      pca_plot2(c(1,4)) + pca_plot2(c(2,4)) + pca_plot2(c(3,4)) +
      pca_plot2(c(1,5)) + pca_plot2(c(2,5)) + pca_plot2(c(3,5)) + pca_plot2(c(4,5)) +
      plot_layout(design = design,
                  guides = "collect")
  }
  
  message("Making gridded PCA DimPlots, colored by variables of interest")
  
  Idents(obj) <- "tissue"
  p <- pca_grid("tissue")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_tissue.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)
  
  # Idents(obj) <- "orig.ident"
  # p <- pca_grid("orig.ident")
  # ggsave(p,
  #        filename = paste0(plots_dir, "pca_dimplot_sample.png"),
  #        units = "in", dpi = 600,
  #        height = 10, width = 24)
  
  Idents(obj) <- "id"
  p <- pca_grid("id")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_id.png"),
         units = "in", dpi = 600,
         height = 10, width = 15)
  
  Idents(obj) <- "Batch"
  p <- pca_grid("Batch")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_batch.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)
  
  Idents(obj) <- "Group"
  p <- pca_grid("Group")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_group.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)
  
  Idents(obj) <- "site"
  p <- pca_grid("site")
  
ggsave(p,
       filename = paste0(plots_dir, tissue, "_pca_dimplot_site.png"),
       units = "in", dpi = 600,
       height = 10, width = 12)
}


# Cluster cells and compute UMAP ------------------------------------------

for (i in seq_along(obj_list)){
  
  obj <- obj_list[[i]]
  
  tissue <- names(obj_list)[i]
  
  obj <- obj %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30)
  
  # obj <- obj %>%
  #   FindClusters(resolution = 0.3,
  #                cluster.name = "seurat_clusters")
  
  # DimPlot_scCustom(obj,
  #                  reduction = "umap",
  #                  group.by = "orig.ident")
  # ggsave(paste0(plots_dir, "umap_dimplot_sample.png"),
  #        units = "in", dpi = 600,
  #        height = 5, width = 10)
  
  DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "tissue",
                   colors_use = JCO_Four())
  ggsave(paste0(plots_dir, tissue, "_umap_dimplot_tissue.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
  
  DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "Batch",
                   colors_use = DiscretePalette_scCustomize(6,
                                                            palette = "ditto_seq"))
  ggsave(paste0(plots_dir, tissue, "_umap_dimplot_batch.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
  
  DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "Group",
                   colors_use = JCO_Four())
  ggsave(paste0(plots_dir, tissue, "_umap_dimplot_group.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)
  
  DimPlot_scCustom(obj,
                   group.by = "id",
                   reduction = "umap")
  ggsave(paste0(plots_dir, tissue, "_umap_dimplot_id.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
  
  DimPlot_scCustom(obj,
                   group.by = "site",
                   reduction = "umap",
                   colors_use = JCO_Four())
  ggsave(paste0(plots_dir, tissue, "_umap_dimplot_site.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
  
}

# DimPlot_scCustom(obj,
#                  reduction = "umap",
#                  group.by = "seurat_clusters")
# ggsave(paste0(plots_dir, "umap_dimplot_cluster.png"),
#        units = "in", dpi = 600,
#        height = 5, width = 6)

# FeaturePlot_scCustom(obj,
#                      reduction = "umap",
#                      features = c("GFAP", "ITGAM", "PLP1",
#                                   "CD3E"))
# ggsave(paste0(plots_dir, "featureplot_gfap_itgam_plp1_cd3e.png"),
#        units = "in", dpi = 600,
#        height = 6, width = 6)
# 
# FeaturePlot_scCustom(obj,
#                      reduction = "umap",
#                      features = c("MS4A1", "DCN", 
#                                   "RBFOX3", "MYOG"))
# ggsave(paste0(plots_dir, "featureplot_ms4a1_dcn_rbfox3_myog.png"),
#        units = "in", dpi = 600,
#        height = 6, width = 6)
# 
# data_out_dir <- paste0(b1042, "data/09_sketch_umap/")
# dir.create(data_out_dir,
#            showWarnings = F,
#            recursive = T)
# 
# saveRDS(obj,
#         paste0(data_out_dir, "obj.rds"))
