suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(BPCells)
  library(patchwork)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")


# Make object manifest ----------------------------------------------------

subclustering_targets <- list(
  list(tissue = "brain", cell_types = "EN"),
  list(tissue = "brain", cell_types = "IN"),
  list(tissue = "brain", cell_types = "Myeloid"),
  list(tissue = "brain", cell_types = "Lymphocyte"),
  list(tissue = "brain", cell_types = "Vascular"),
  list(tissue = "muscle", cell_types = "Muscle fiber"),
  list(tissue = "muscle", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell")),
  list(tissue = "muscle", cell_types = c("Macrophage", "Monocyte", "Neutrophil")),
  list(tissue = "muscle", cell_types = c("T-cell", "NK cell")),
  list(tissue = "sc", cell_types = "Myeloid"),
  list(tissue = "sc", cell_types = "Lymphocyte"),
  list(tissue = "sc", cell_types = c("IN", "EN", "SN", "MN")), # not changing these labels at all
  list(tissue = "sc", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell"))
)

# Assemble object, load marker genes --------------------------------------

i <- 1

tissue <- subclustering_targets[[i]]$tissue
cell <- paste(subclustering_targets[[i]]$cell_types, collapse = "_")

data_in_dir <- paste0("data/15_subclustering2/", tissue, "/", cell, "/")
plots_dir <- paste0("plots/15_subclustering2/", tissue, "/", cell, "/")
tab_dir <- paste0("tab_data/15_subclustering2/", tissue, "/", cell, "/")

data_mat <- open_matrix_dir(paste0(data_in_dir, "/bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
umap <- readRDS(paste0(data_in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony_umap"]] <- umap

markers <- read.csv(paste0(tab_dir, "markers.csv"))

# Get clustering resolution -----------------------------------------------

graph_modularity <- read.csv(paste0(tab_dir, "graph_modularity.csv"))

res <- graph_modularity %>% 
  arrange(desc(modularity)) %>% 
  pull(resolution) %>% 
  .[1]

Idents(obj) <- paste0("res", res, "_clusters")

# Plotting ----------------------------------------------------------------

DimPlot_scCustom(obj) 

combined_plot <- function(gene){
  p1 <- FeaturePlot_scCustom(obj,
                             features = gene)
  
  p2 <- VlnPlot_scCustom(obj,
                         features = gene) + 
    NoLegend()
  
  p <- p1 + p2 + 
    plot_layout(nrow = 2,
                heights = c(2, 1))
  
  ggsave(p,
         filename = paste0(plots_dir, gene, ".png"),
         units = "in", dpi = 600,
         width = 6, height = 9)
}

gene_lists <- list(
  c("SLC17A7", "SLC17A6", "CUX2", "CBLN2", "LAMP5", "FREM3", "RORB", "C1QL2", "PRSS12", "FOXP2", "TLE4", "SYT6", "NR4A2", "ETV1", "RSPO2", "PCP4",
    "FEZF2", "BCL11B", "POU3F1", "CRYM", "PCP4", "TSHZ2", "NEFH", "VAT1L"),
  c("NPY", "SST", "PVALB", "GAD1", "SLC17A7", "VIP", "LAMP5", "KIT", "CXCL14", "CCK", "RELN", "NEFH"),
  c("CD163", "VCAN", "MRC1", "CCR2", "MOBP", "GFAP", "MKI67", "S100A8", "LYVE1", "CD14", "CSF1R"),
  c("CD3E", "CD8A", "NKG7", "TRAC", "FCGR3A", "ITGAM", "VCAN", "S100A8", "CSF1R", "CD163", "P2RY12"),
  c("PODXL", "CLDN5", "RGS5", "ACTA2", "PECAM1", "PDGFRB", "VWF", "CSPG4", "MCAM"),
  c("MYH7", "TNNT1", "TNNI1", "MYL2", "ATP2A2", "MYH2", "MYH1", "TNNT3", "TNNI2", "MYL1", "ATP2A1", "CHRNA1", "CHRNG", "MYOG"),
  c("PODXL", "CLDN5", "RGS5", "ACTA2", "PECAM1", "PDGFRB", "VWF", "CSPG4", "MCAM", "MYH11"),
  c("CD163", "CD300E", "LYZ", "LYVE1", "ITGAX", "S100A8", "LTF", "TREM2", "VCAN", "CD14", "FCGR3A", "ITGAM", "MS4A1", "TRAC"),
  c("CD3E", "CD8A", "NKG7", "TRAC", "FCGR3A", "ITGAM", "VCAN", "S100A8", "CSF1R", "CD163", "MS4A1", "NCAM1"),
  c("CD163", "VCAN", "MRC1", "CCR2", "MOBP", "GFAP", "MKI67", "S100A8", "LYVE1", "CD14", "CSF1R"),
  c("CD3E", "CD8A", "NKG7", "TRAC", "FCGR3A", "ITGAM", "VCAN", "S100A8", "CSF1R", "CD163", "P2RY12"),
  c("SLC17A7", "SLC17A6", "GAD1", "SST", "LAMP5", "CCK", "NEFH", "CHAT", "SLC18A3", "FAT2"),
  c("PODXL", "CLDN5", "RGS5", "ACTA2", "PECAM1", "PDGFRB", "VWF", "CSPG4", "MCAM", "MYH11")
)

for (j in gene_lists[[i]]) {
  combined_plot(j)
}

# Inspecting marker genes -------------------------------------------------

markers %>% 
  filter(cluster == 6 & 
           pct.1 > 0.3) %>% 
  arrange(desc(avg_log2FC)) %>% 
  head(30)

# Updating annotation sheet -----------------------------------------------

annotations <- read.csv(paste0(tab_dir, "annotations.csv"))

n_clusters <- length(unique(Idents(obj)))

annotations %>% 
  filter(cluster > 0 & cluster <= n_clusters) %>%
  mutate(cell_type = case_when(cluster %in% c(17, 18) ~ "Remove",
                               cluster %in% c(3, 2, 6, 13, 15) ~ "Type II MF",
                               cluster %in% c(4, 7, 8, 9) ~ "Type I MF",
                               cluster %in% c(1, 5, 10, 11, 12, 14) ~ "Denervated MF",
                               cluster %in% c(16) ~ "Proliferating MF")) %>%
  write.csv(file = paste0(tab_dir, "annotations.csv"),
            row.names = F,
            quote = F)


# Extra stuff for trickier subtypes ---------------------------------------

dittoBarPlot(obj,
             var = "DF.unadj",
             group.by = paste0("res", res, "_clusters"))

dittoDotPlot(obj,
             vars = c("MYH7", "TNNT1", "TNNI1", "MYL2", "ATP2A2", "MYH2", "MYH1", "TNNT3", "TNNI2", "MYL1", "ATP2A1"),
             group.by = paste0("res", res, "_clusters"))
ggsave(filename = paste0(plots_dir, "fiber_typing_genes_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

dittoDotPlot(obj,
             vars = c("CHRNA1", "CHRNG", "CHRND", "MYOG", "MYF6", "NCAM1", "FBXO32", "TRIM63", "GDF15"),
             group.by = paste0("res", res, "_clusters"))
ggsave(filename = paste0(plots_dir, "fiber_health_genes_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

obj@meta.data <- obj@meta.data %>% 
  mutate(mf = case_when(res1_clusters %in% c(2, 13, 15, 6) ~ "Type I",
                        res1_clusters %in% c(3, 4, 7, 8) ~ "Type II",
                        res1_clusters %in% c(5, 10, 11, 14, 12, 9, 1) ~ "Denervated",
                        res1_clusters %in% c(16) ~ "Proliferating",
                        res1_clusters %in% c(17, 18) ~ "Remove"))

obj@meta.data <- obj@meta.data %>% 
  mutate(mf = case_when(res1_clusters %in% c(17, 18) ~ "Remove",
                        res1_clusters %in% c(3, 2, 6, 13, 15) ~ "fast",
                        res1_clusters %in% c(4, 7, 8, 9) ~ "slow",
                        res1_clusters %in% c(16) ~ "proliferating",
                        res1_clusters %in% c(1, 5, 10, 11, 12, 14) ~ "injured"))

DimPlot_scCustom(obj,
                 group.by = "mf")

dittoDotPlot(obj,
             vars = c("CUX2", "CBLN2", "LAMP5", "FREM3", "RORB", "C1QL2", "PRSS12", "FOXP2", "TLE4", "SYT6", "NR4A2", "ETV1", "RSPO2", "PCP4",
                      "FEZF2", "BCL11B", "POU3F1", "CRYM", "TSHZ2", "NEFH", "VAT1L"),
             group.by = paste0("res", res, "_clusters"))
ggsave(filename = paste0(plots_dir, "en_subtype_genes_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 12)

table(obj$orig.ident,
      Idents(obj))

Stacked_VlnPlot(obj,
                features = c("nCount_RNA", "nFeature_RNA", "percent_mito"))
