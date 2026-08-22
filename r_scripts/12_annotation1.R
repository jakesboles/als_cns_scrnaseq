suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(dittoSeq)
  library(BPCells)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_in_dir0 <- "data/09_integration/"

clusters_in_dir0 <- "data/10_clustering/"

tab_out_dir0 <- "tab_data/12_annotation1/"

plots_dir0 <- "plots/12_annotation1/"

combined_plot <- function(s, gene, res){
  p1 <- FeaturePlot_scCustom(s,
                             features = gene,
                             reduction = "umap",
                             raster = T,
                             raster.dpi = c(900, 900),
                             pt.size = 0.05)
  p2 <- VlnPlot_scCustom(s,
                         features = gene,
                         group.by = paste0("res", res, "_clusters"),
                         raster = T,
                         raster.dpi = 900) + 
    NoLegend()
  
  p <- p1 + p2 + 
    plot_layout(nrow = 2,
                heights = c(2, 1))
  
  ggsave(p,
         filename = paste0(plots_dir, gene, ".png"),
         units = "in", dpi = 300,
         height = 9, width = 7)
}

# Annotating brain --------------------------------------------------------

tissue <- "brain"

data_in_dir <- paste0(data_in_dir0, tissue, "/")

clusters_in_dir <- paste0(clusters_in_dir0, tissue, "/")

tab_out_dir <- paste0(tab_out_dir0, tissue, "/")
dir.create(tab_out_dir,
           showWarnings = F,
           recursive = T)

plots_dir <- paste0(plots_dir0, tissue, "/")
dir.create(plots_dir,
           showWarnings = F,
           recursive = T)

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(clusters_in_dir, "metadata.rds"))
umap <- readRDS(paste0(clusters_in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["umap"]] <- umap

resolution <- 1

DimPlot_scCustom(obj,
                 reduction = "umap",
                 group.by = paste0("res", resolution, "_clusters"))
ggsave(filename = paste0(plots_dir, "cluster_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

markers <- read.csv(paste0("tab_data/11_findmarkers/", tissue, "/markers.csv"))

top5 <- markers %>% 
  filter(pct.1 > 0.3) %>% 
  Extract_Top_Markers(num_features = 5,
                      make_unique = T,
                      named_vector = F)

p <- dittoDotPlot(obj,
             vars = top5,
             group.by = paste0("res", resolution, "_clusters"))
ggsave(p,
       filename = paste0(plots_dir, "top5_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 20)

for (i in c("C1QA", 'ITGAM', 'P2RY12', 'GFAP', 'AQP4', 'PDGFRA', 'DCN',
            'CEMIP', 'RBFOX3', 'OLIG1', 'MOBP', 'PLP1', 'CLDN5', 'PODXL', 
            'RGS5', 'ACTA2', 'NKG7', 'GAD1', 'SNAP25', 'CNP', 'OPALIN')){
  
  combined_plot(obj,
                i,
                resolution)
}

# Plotting etc for annotations --------------------------------------------

# pick dataset to focus on
# 1 = brain, 2 = muscle, 3 = sc
i <- 3

{clust <- 30
  markers[[i]] %>%
    Add_Pct_Diff() %>%
    filter(cluster == clust & 
             pct.1 > 0.3) %>%
    dplyr::arrange(desc(avg_log2FC)) %>%
    head(n = 15L)
}

FeaturePlot_scCustom(obj_list[[i]],
                     features = "nCount_RNA",
                     reduction = "cca_umap")

VlnPlot_scCustom(obj_list[[i]],
                 features = "CD36")

DimPlot_scCustom(obj_list[[i]],
                 group.by = "DF.unadj",
                 reduction = "cca_umap")

VlnPlot_scCustom(obj_list[[i]],
                 features = "MOBP",
                 group.by = "DF.adj")

brain_genes <- c("C1QA", 'ITGAM', 'P2RY12', 'GFAP', 'AQP4', 'PDGFRA', 'DCN',
                 'CEMIP', 'RBFOX3', 'OLIG1', 'MOBP', 'PLP1', 'CLDN5', 'PODXL', 
                 'RGS5', 'ACTA2', 'NKG7', 'GAD1', 'SNAP25', 'CNP', 'OPALIN')
muscle_genes <- c('ITGAM', 'ITGAX', 'MS4A1', 'CD3E', 'NKG7', 'S100A8', 'APOE', 
                  'APOD', 'CLDN5', 'PODXL', 'RGS5', 'ACTA2', 'SNAP25', 'RBFOX3',
                  'TNNC2', 'TNNT1', 'MYL2', 'LYZ', 'CD300E', 'TPSB2')
sc_genes <- c("C1QA", 'ITGAM', 'P2RY12', 'GFAP', 'AQP4', 'PDGFRA', 'DCN',
              'CEMIP', 'RBFOX3', 'OLIG1', 'MOBP', 'PLP1', 'CLDN5', 'PODXL', 
              'RGS5', 'ACTA2', 'NKG7', 'GAD1', 'SNAP25', 'CFAP157', 'TRAC', 
              'CNP', 'MPZ', 'OPALIN')

# for (gene in sc_genes){
{
  gene <- str_to_upper("gpr17")
  p1 <- FeaturePlot_scCustom(obj_list[[i]],
                             features = gene,
                             reduction = "cca_umap")
  p2 <- VlnPlot_scCustom(obj_list[[i]],
                         features = gene) +
    NoLegend()
  p <- p1 / p2 +
    plot_layout(heights = c(2, 1))
  
  p
  
  # ggsave(p,
  #        filename = paste0(plots_dir, names(obj_list)[i],
  #                          "_", gene, ".png"),
  #        units = "in", dpi = 600,
  #        bg = "white",
  #        height = 10, width = 10)
}
# }

Cluster_Highlight_Plot(obj_list[[i]], 
                       cluster_name = c("28", "30"),
                       reduction = "cca_umap",
                       highlight_color = c("red", "green2"))

dittoBarPlot(obj_list[[i]],
             group.by = "res4_clusters",
             var = "Group")

VlnPlot(obj_list[[i]],
        group.by = "DF.unadj",
        feature = "nCount_RNA")


# Subclustering dubious clusters ------------------------------------------
j <- 3
check <- c(13, 20, 31, 34)

sub <- subset(obj_list[[j]],
              subset = res2_clusters %in% check) %>%
  RunPCA(npcs = 100) %>%
  FindNeighbors(reduction = "cca_pca",
                dims = 1:30,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "cca_pca",
                dims = 1:30,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          # reduction = "cca_pca",
          nn.name = "sketch.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "cca_umap")

# DimPlot_scCustom(sub,
#                  reduction = "cca_umap")

sub <- FindClusters(sub,
                    resolution = 0.1)

DimPlot_scCustom(sub,
                 reduction = "cca_umap") +
  DimPlot_scCustom(sub,
                   reduction = "cca_umap",
                   group.by = "res2_clusters",
                   colors_use = JCO_Four())

DimPlot_scCustom(sub,
                 reduction = "cca_umap",
                 group.by = "DF.adj")

sub_markers <- FindAllMarkers(sub)

sub_markers %>%
  filter(cluster == 0) %>% 
  filter(pct.1 > 0.4) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)

FeaturePlot_scCustom(sub,
                     reduction = "cca_umap",
                     features = "SLC6A5")

# Running more DE to tease apart similar clusters -------------------------
k <- 1

de <- FindMarkers(obj_list[[k]],
                  ident.1 = 37, 
                  ident.2 = c(3, 10, 25, 26, 30))

head(de)

# AFTER ANNOTATING --------------------------------------------------------


annots <- list.files(tab_in_dir, 
                     full.names = T)[str_detect(list.files(tab_in_dir), "annotation")]

annots <- map(annots, 
              Pull_Cluster_Annotation)

for (i in seq_along(obj_list)){
  obj_list[[i]] <- Rename_Clusters(obj_list[[i]], 
                                   new_idents = annots[[i]]$new_cluster_idents)
}

DimPlot_scCustom(obj_list[[i]])

DimPlot_scCustom(obj_list[[3]],
                 reduction = "cca_umap",
                 label = F)
