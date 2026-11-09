library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(paletteer)
library(BPCells)
library(ggbeeswarm)
library(scCustomize)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

source("r_scripts/plot_theme.R")

in_dir <- paste0("results/wgcna_consensus/Microglia/")
data_dir <- paste0("data/19_subclustering3/microglia/")

plots_dir <- "figures/figure2/"

mois <- c("blue", "turquoise")
moi_names <- c("ALS module 1", "ALS module 2")
# moi_colors <- c("#3B3B98", "#8C510A")
moi_colors <- c("#B5175E", "#007A87")

# kME plots ---------------------------------------------------------------

modules <- read.csv(paste0(in_dir, "modules.csv"))

for (i in seq_along(mois)){
  
  modules %>% 
    filter(color == mois[i]) %>% 
    arrange(desc(!!sym(paste0("kME_", mois[i])))) %>% 
    mutate(gene_name = fct_inorder(gene_name)) %>% 
    slice_head(n = 25) %>%
    ggplot(aes(x = !!sym(paste0("kME_", mois[i])),
               y = gene_name)) + 
    geom_col(color = "black",
             fill = moi_colors[i],
             linewidth = 0.25) + 
    scale_y_discrete(limits = rev) + 
    scale_x_continuous(expand = c(0, 0),
                       n.breaks = 4) + 
    labs(x = "kME") +
    ggtitle(paste0("Module ", i, " members")) +
    theme_linedraw() + 
    theme_als_project(title_size = 7) +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(margin = margin(t = 0.025, unit = "cm")),
          legend.title = element_blank(),
          legend.position = "none",
          axis.ticks.y = element_blank(),
          plot.title = element_text(margin = margin(b = 0.05, unit = "cm")),
          # axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
          strip.text = element_text(face = "bold", color = "black"),
          axis.text.y = element_text(size = 5))
  save_figure(last_plot(),
              filename = paste0(plots_dir, "wgcna_module", i, "_kme_bars.pdf"),
              height = 2, width = 1.25)
}

# Make Seurat object ------------------------------------------------------

data_mat <- open_matrix_dir(paste0(data_dir, "bpcells_data"))
meta <- readRDS(paste0(data_dir, "metadata.rds"))
harmony <- readRDS(paste0(data_dir, "harmony.rds"))
umap <- readRDS(paste0(data_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony
obj[["umap"]] <- umap

obj$group <- factor(obj$group, levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))
obj$tissue <- factor(obj$tissue,
                     levels = c("Motor cortex", "Cervical spinal cord"),
                     labels = c("MCX", "SC"))

scores <- read.csv(paste0(in_dir, "module_scores_ucell.csv"))

# FeaturePlots ------------------------------------------------------------

obj <- AddMetaData(obj,
                   scores)

module_grads <- list(
  c("#E6E6E6", "#E7A3C0", "#B5175E", "#5E0A30"),
  c("#E6E6E6", "#9CCFD4", "#007A87", "#003E45")
)

embeddings <- Embeddings(obj, reduction = "umap")

featureplot_df <- cbind(embeddings, obj@meta.data)

for (i in seq_along(mois)){
  p <- featureplot_df %>% 
    ggplot(aes(x = harmonyumap_1,
               y = harmonyumap_2)) + 
    ggrastr::rasterize(
      geom_point(aes(color = !!sym(paste0(mois[i], "_UCell_kNN"))),
                 size = 0.1),
      dpi = 1600) + 
    scale_color_gradientn(
      colours = module_grads[[i]],
      # values = rescale(breaks, from = c(min_logfc, max_logfc)),
      # limits = c(min_logfc, max_logfc),
      oob = squish) + 
    ggtitle(moi_names[i]) +
    # umap_labs + 
    theme_void() +
    theme_als_project(title_size = 7) + 
    umap_theme +
    theme(legend.text = element_blank(),
          legend.box.spacing = unit(0.1, "cm"),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.ticks = element_blank())
  save_figure(p,
              filename = paste0(plots_dir, "wgcna_module", i, "_umap.pdf"),
              height = 3, width = 2)
}

# Expression plots ---------------------------------------------------

pb <- scores %>% 
  group_by(orig.ident, group, tissue) %>% 
  mutate(across(where(is.numeric), median)) %>% 
  distinct(orig.ident,
           .keep_all = T)

pb %>%
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) + 
  geom_quasirandom(aes(fill = group),
                   shape = 21,
                   size = 4,
                   alpha = 0.7) + 
  stat_summary(fun = mean,
               geom = "crossbar") +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               linewidth = 1.2,
               width = 0.6) +
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  labs(y = paste0(str_to_title(mois[i]), " module score")) +
  theme_linedraw(base_size = 12) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(in_dir, mois[i], "_expression_pb.png"),
       units = "in", dpi = 600,
       height = 3, width = 4)
}

