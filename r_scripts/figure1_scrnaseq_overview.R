library(ggplot2)
library(tidyverse)
library(patchwork)
library(Seurat)
library(scCustomize)
library(BPCells)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir0 <- "figures/"

source("r_scripts/plot_theme.R")

# Assemble full 3-way object ----------------------------------------------

in_dir <- "data/18_full_integration/all_tissues/"

mat <- open_matrix_dir(paste0(in_dir, "bpcells_data"))

meta <- readRDS(paste0(in_dir, "metadata.rds"))
harmony <- readRDS(paste0(in_dir, "harmony.rds"))
umap <- readRDS(paste0(in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = mat, meta.data = meta, assay = "RNA")
obj[["harmony"]] <- harmony
obj[["umap"]] <- umap

obj$tissue <- factor(obj$tissue,
                     levels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"))
obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

# Main UMAPs --------------------------------------------------------------

plots_dir <- paste0(plots_dir0, "figure1/")
dir.create(plots_dir,
           recursive = T,
           showWarnings = F)

umap_labs <- labs(y = "UMAP 2",
                  x = "UMAP 1")

umap_theme <- theme(axis.text = element_blank(),
                    axis.ticks = element_blank(),
                    aspect.ratio = 1)

obj@meta.data <- obj@meta.data %>% 
  mutate(cell_type4 = case_when(str_detect(cell_type3, "MF") ~ "Muscle fiber",
                                str_detect(cell_type3, "EN|IN|MN|SN") ~ "Neuron",
                                .default = cell_type3))

umap_celltype <- DimPlot_scCustom(obj,
                 group.by = "cell_type4",
                 pt.size = 2,
                 raster.dpi = c(900, 900),
                 colors_use = "polychrome") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  umap_labs +
  theme_als_project() + 
  umap_theme +
  theme(plot.title = element_blank())
save_figure(umap_celltype,
            filename = paste0(plots_dir, "umap_celltype.pdf"),
            height = 7, width = 10)
# ggsave(filename = paste0(plots_dir, "umap_celltype.png"),
       # units = "in", dpi = 600,
       # height = 7, width = 10)

umap_tissue <- DimPlot_scCustom(obj,
                 group.by = "tissue",
                 pt.size = 2,
                 raster.dpi = c(900, 900),
                 colors_use = JCO_Four()) + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  umap_labs +
  theme_als_project() +
  umap_theme +
  theme(plot.title = element_blank())
save_figure(umap_tissue,
            filename = paste0(plots_dir, "umap_tissue.pdf"),
            height = 3.5, width = 6)
# ggsave(filename = paste0(plots_dir, "umap_tissue.png"),
       # units = "in", dpi = 600,
       # height = 7, width = 10)

umap_group <- DimPlot_scCustom(obj,
                                group.by = "group",
                                pt.size = 2,
                                raster.dpi = c(900, 900),
                                colors_use = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  umap_labs +
  theme_als_project() +
  umap_theme +
  theme(plot.title = element_blank())
save_figure(umap_group,
            filename = paste0(plots_dir, "umap_group.pdf"),
            height = 3.5, width = 5)

# patchworked version to see if this is better for assembling in Illustrator
design <- "
AAAABB
AAAABB
AAAACC
AAAACC
"

umap_celltype + umap_tissue + umap_group + 
  plot_layout(design = design)
save_figure(last_plot(),
            filename = paste0(plots_dir, "umap_3way.pdf"),
            height = 7, width = 15)

# Highlight neurons and MFs since I collapsed their labels ----------------

Idents(obj) <- "cell_type4"

Cluster_Highlight_Plot(obj,
                       cluster_name = "Neuron",
                       highlight_color = "midnightblue",
                       raster.dpi = c(900, 900)) + 
  theme(legend.position = "none")
ggsave(filename = paste0(plots_dir, "neuron_highlight.png"),
       units = "in", dpi = 600,
       height = 4, width = 4.5)

# UMAPs for David 9/16/26 -------------------------------------------------

in_dir <- "data/18_full_integration/brain_sc/"

mat <- open_matrix_dir(paste0(in_dir, "bpcells_data"))

meta <- readRDS(paste0(in_dir, "metadata.rds"))
harmony <- readRDS(paste0(in_dir, "harmony.rds"))
umap <- readRDS(paste0(in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = mat, meta.data = meta, assay = "RNA")
obj[["harmony"]] <- harmony
obj[["umap"]] <- umap

obj$tissue <- factor(obj$tissue,
                     levels = c("Motor cortex", "Cervical spinal cord"))
obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

obj@meta.data <- obj@meta.data %>% 
  mutate(cell_type4 = case_when(tissue == "Cervical spinal cord" &
                                  cell_type3 %in% c("EN", "MN", "SN", "IN") ~ paste0("Spinal ", cell_type3),
                                .default = cell_type3))

p_tissue <- DimPlot_scCustom(obj,
                             group.by = "tissue",
                             reduction = "umap",
                             raster.dpi = c(900, 900),
                             pt.size = 2) + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  ggtitle("Tissue")
ggsave(p_tissue,
       filename = paste0(plots_dir, "cns_umap_tissue.pdf"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p_group <- DimPlot_scCustom(obj,
                            group.by = "group",
                            reduction = "umap",
                            raster.dpi = c(900, 900),
                            pt.siz = 2) + 
  scale_color_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) +
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  ggtitle("Group")
ggsave(p_group,
       filename = paste0(plots_dir, "cns_umap_group.pdf"),
       units = "in", dpi = 600,
       height = 5, width = 6)

p_cell <- DimPlot_scCustom(obj,
                           group.by = "cell_type4",
                           raster.dpi = c(900, 900),
                           pt.size = 2) + 
  guides(color = guide_legend(ncol = 2,
                              override.aes = list(size = 4))) + 
  ggtitle("Cell type")
ggsave(p_cell,
       filename = paste0(plots_dir, "cns_umap_celltype.pdf"),
       units = "in", dpi = 600,
       height = 5, width = 8)

