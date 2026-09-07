library(Seurat)
library(scCustomize)
library(tidyverse)
library(BPCells)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "figures/"


# QC metrics --------------------------------------------------------------

df <- readRDS("data/02_qc1/metadata.rds")

str(df)

df <- df %>%
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")))

theme <- theme_linedraw(base_size = 12) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))

df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = nCount_RNA)) + 
  facet_wrap(. ~ tissue,
             ncol = 1) + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "# UMIs per cell") +
  scale_y_log10() + 
  theme
ggsave(filename = paste0(plots_dir, "ncount.png"),
       units = "in", dpi = 600,
       height = 8, width = 12)

df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = nFeature_RNA)) + 
  facet_wrap(. ~ tissue,
             ncol = 1) + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "# genes per cell") +
  scale_y_log10() + 
  theme
ggsave(filename = paste0(plots_dir, "nfeature.png"),
       units = "in", dpi = 600,
       height = 8, width = 12)

df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = percent_mito)) + 
  facet_wrap(. ~ tissue,
             ncol = 1) + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "% mitochondrial genes per cell") +
  theme
ggsave(filename = paste0(plots_dir, "mito.png"),
       units = "in", dpi = 600,
       height = 8, width = 12)

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

str(obj@meta.data)

obj@meta.data <- obj@meta.data %>% 
  mutate(cell_type4 = case_when(str_detect(cell_type3, "MF") ~ "Muscle fiber",
                                str_detect(cell_type3, "EN|IN|MN|SN") ~ "Neuron",
                                .default = cell_type3))

length(unique((obj$cell_type4)))

DimPlot_scCustom(obj,
                 group.by = "cell_type4",
                 pt.size = 2,
                 raster.dpi = c(900, 900),
                 colors_use = "polychrome") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  theme(plot.title = element_blank())
ggsave(filename = paste0(plots_dir, "umap_celltype.png"),
       units = "in", dpi = 600,
       height = 7, width = 10)

DimPlot_scCustom(obj,
                 group.by = "tissue",
                 pt.size = 2,
                 raster.dpi = c(900, 900),
                 colors_use = JCO_Four()) + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  theme(plot.title = element_blank())
ggsave(filename = paste0(plots_dir, "umap_tissue.png"),
       units = "in", dpi = 600,
       height = 7, width = 10)

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
