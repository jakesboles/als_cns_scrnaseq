library(Seurat)
library(scCustomize)
library(tidyverse)
library(patchwork)
library(dittoSeq)

proj_dir <- "/projects/b1169/boles/als_multitissue_scfrp/"

data_in_dir <- paste0(proj_dir, "data/11_clustering")
files_in <- list.files(data_in_dir,
                       full.names = T)

plots_dir <- paste0(proj_dir, "plots/13_annotation/")
dir.create(plots_dir,
           showWarnings = F,
           recursive = T)

obj_list <- list()
length(obj_list) <- 3

obj_list <- map(files_in, 
                readRDS)

names(obj_list) <- c("brain", "muscle", "sc")

tab_in_dir <- paste0(proj_dir, "tab_data/12_findmarkers/")
markers <- list.files(tab_in_dir, 
                      full.names = T)[str_detect(list.files(tab_in_dir), "markers")]
markers <- markers[str_detect(markers, "OLD", negate = T)]  
markers <- map(markers, read.csv)

names(markers) <- names(obj_list)

Idents(obj_list[[1]]) <- "res1_clusters"
Idents(obj_list[[2]]) <- "res1.2_clusters"
Idents(obj_list[[3]]) <- "res1.4_clusters"


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
