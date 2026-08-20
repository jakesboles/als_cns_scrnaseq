library(Seurat)
library(tidyverse)
# library(cca)
library(scCustomize)

b1169 <- "/projects/b1169/boles/als_multitissue_scfrp/"
b1042 <- "/projects/b1042/Gate_Lab/boles/als_multitissue/"
p31535 <- "/projects/p31535/boles/als_multitissue_scfrp/"

message("Loading object")
t0 <- Sys.time()

obj <- readRDS(paste0(b1169, "data/08_sketch_pca/sc_obj.rds"))

Sys.time() - t0

data_out_dir <- paste0(b1169, "data/10_integration/")
dir.create(data_out_dir,
           showWarnings = F,
           recursive = T)

plots_dir <- paste0(b1169, "plots/10_integration/")
dir.create(plots_dir,
           showWarnings = F,
           recursive = T)

# obj0$tissue %>% unique()
# tissue_of_interest <- "Cervical spinal cord"
file <- "sc"
# 
# t0 <- Sys.time()
# 
# obj <- subset(obj0,
#               subset = tissue == tissue_of_interest)
# 
# Sys.time() - t0

message("Integrating samples using CCA")

t0 <- Sys.time()

# theta <- 2

obj <- IntegrateLayers(obj,
                       method = "CCAIntegration",
                       orig.reduction = "pca",
                       new.reduction = "cca_pca",
                       k.anchor = 20,
                       reference = which(Layers(obj, search = "data") %in% 
                                           c("data.GWF21-56_s", "data.GBB-18-13_s")),
                       dims = 1:20)

t1 <- Sys.time()
t1 - t0

message("Computing NN graphs, UMAP, and clustering")

obj <- FindNeighbors(obj,
                     dims = 1:20,
                     reduction = "cca_pca") %>%
  FindClusters(resolution = 0.3,
               cluster.name = "seurat_clusters") %>%
  RunUMAP(dims = 1:20,
          reduction = "cca_pca",
          reduction.name = "cca_umap",
          reduction.key = "cca_umap_")

t2 <- Sys.time()
t2 - t1

message("Making plots")

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "tissue",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_umap_dimplot_tissue.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "Batch",
                      colors_use = DiscretePalette_scCustomize(6,
                                                               palette = "ditto_seq"))
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_umap_dimplot_batch.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "Group",
                      colors_use = JCO_Four())
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_umap_dimplot_group.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- DimPlot_scCustom(obj,
                      group.by = "id",
                      reduction = "cca_umap")
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_umap_dimplot_id.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)

p <- DimPlot_scCustom(obj,
                      reduction = "cca_umap",
                      group.by = "seurat_clusters")
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_umap_dimplot_cluster.png"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p <- FeaturePlot_scCustom(obj,
                          reduction = "cca_umap",
                          features = c("GFAP", "ITGAM", "PLP1",
                                       "CD3E"))
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_featureplot_gfap_itgam_plp1_cd3e.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

p <- FeaturePlot_scCustom(obj,
                          reduction = "cca_umap",
                          features = c("MS4A1", "DCN", 
                                       "RBFOX3", "MYOG"))
ggsave(p,
       filename = paste0(plots_dir, file,
                         "_cca_featureplot_ms4a1_dcn_rbfox3_myog.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)

message("Saving CCA-intergrated data")

t3 <- Sys.time()

saveRDS(obj,
        paste0(data_out_dir, file, "_cca_obj.rds"))

Sys.time() - t3