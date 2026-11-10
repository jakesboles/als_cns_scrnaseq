library(Seurat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(paletteer)
library(BPCells)
library(ggbeeswarm)
library(scCustomize)
library(ggpubr)
library(scales)

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

scores <- scores %>% 
  mutate(group = factor(group, 
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")),
         tissue = factor(tissue,
                         levels = c("Motor cortex", "Cervical spinal cord"),
                         labels = c("MCX", "SC")))



pb <- scores %>% 
  group_by(orig.ident, group, tissue) %>% 
  mutate(across(where(is.numeric), median)) %>% 
  distinct(orig.ident,
           .keep_all = T)

p1 <- pb %>%
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) + 
  geom_quasirandom(aes(color = group),
                   shape = 16,
                   size = 1,
                   alpha = 0.7) + 
  stat_summary(fun = mean,
               geom = "crossbar",
               linewidth = 0.1) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               linewidth = 0.12,
               width = 0.6) +
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  scale_color_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  # labs(y = paste0(str_to_title(mois[i]), " module score")) +
  theme_linedraw() +
  theme_als_project() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p2 <- scores %>% 
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) + 
  geom_violin(aes(fill = group),
              linewidth = 0.25) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  # labs(y = paste0(str_to_title(mois[i]), " module score")) +
  theme_linedraw() +
  theme_als_project() +
  rotate_x_text(45) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.ticks.x = element_blank(),
        strip.text = element_blank())

p1 + p2 + 
  plot_layout(ncol = 1) &
  theme(panel.spacing = unit(0, "pt"),
        plot.margin = margin(t = 0, b =))
save_figure(last_plot(),
            filename = paste0(plots_dir, "wgcna_module", i, "_expression_patch.pdf"),
            height = 2, width = 1.5)

for (i in seq_along(mois)){
  p <- ggplot(data = scores,
              aes(x = group,
                  y = !!sym(paste0(mois[i], "_UCell_kNN")))) + 
    geom_violin(aes(fill = group),
                linewidth = 0.25,
                show.legend = F) + 
    stat_summary(data = pb,
                 fun = mean,
                 geom = "crossbar",
                 linewidth = 0.1) +
    # stat_summary(data = pb,
    #              fun.data = mean_se,
    #              geom = "errorbar",
    #              linewidth = 0.2,
    #              width = 0.6) +
    geom_quasirandom(data = pb,
                     aes(color = group),
                     shape = 16,
                     size = 1,
                     alpha = 0.6,
                     show.legend = F) + 
    facet_wrap(. ~ tissue,
               nrow = 1) +
    geom_pwc(data = pb, 
             method = "tukey_hsd",
             hide.ns = T,
             label = "p.signif",
             label.size = 5 / .pt,
             tip.length = 0,
             vjust = 0.5,
             step.increase = 0.05) +
    scale_fill_manual(values = c("#D6D2CF", "#79D073", "#E373FF")) + 
    scale_color_manual(values = c("#514E4C", "#065500", "#660080")) +
    scale_y_log10() +
    labs(y = "Module score") + 
    ggtitle(paste0("Module ", i, " expression")) +
    theme_linedraw() + 
    theme_als_project(title_size = 7) +
    rotate_x_text(45) +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = "none",
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0, "cm"),
          strip.text = element_text(margin = margin(t = 1, b = 1, r = 2, l = 2)),
          plot.margin = margin(0, 0, 0, 0),
          plot.title = element_text(margin = margin(b = 1, "cm")))
  save_figure(p,
              filename = paste0(plots_dir, "wgcna_module", i, "_expression.pdf"),
              height = 2, width = 1.5)
}

# ORA result plots --------------------------------------------------------

ora <- read.csv(paste0(in_dir, mois[i], "_gsea.csv"))

wrap_n <- c(17, 26)

ora %>%
  filter(Count > 10 & p.adjust < 0.05 &
           str_detect(X, "GOBP|GOMF|GOCC")) %>%
  mutate(db = str_split_i(X, "_", i = 1) %>% 
           str_remove_all("GO"),
         path = str_remove_all(X, "GOBP_|GOMF_|GOCC_")) %>%
  mutate(path = paste0(str_replace_all(path, "_", " "), " (", db, ")")) %>%
  arrange(desc(FoldEnrichment)) %>%
  mutate(path = fct_inorder(path)) %>%
  slice_head(n = 6) %>% 
  ggplot(aes(y = path,
             x = FoldEnrichment)) +
  geom_col(aes(fill = Count),
           linewidth = 0.25,
           color = "black") + 
  scale_fill_gradientn(
    colours = module_grads[[i]][-1]) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(labels = label_wrap(width = wrap_n[i]),
                   limits = rev)  +
  labs(x = "Fold enrichment",
       fill = "# genes\nshared") +
  ggtitle(paste0("Module ", i, " pathway enrichment")) +
  theme_linedraw() + 
  theme_als_project(title_size = 7) +
  theme(plot.margin = margin(0, 0, 0, 0),
        plot.title = element_text(margin = margin(b = 1, "cm")),
        legend.box.spacing = unit(0.05, "cm"),
        axis.title.y = element_blank(),
        axis.title.x = element_text(margin = margin(t = 0.025, unit = "cm")),
        axis.text.y = element_text(size = 5),
        axis.ticks.y = element_blank())
save_figure(last_plot(),
            filename = paste0(plots_dir, "wgcna_module", i, "_ora.pdf"),
            height = 2, width = 2.5)
