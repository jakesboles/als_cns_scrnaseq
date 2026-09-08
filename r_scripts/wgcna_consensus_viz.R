library(tidyverse)
library(ggplot2)
library(scCustomize)
library(paletteer)
library(BPCells)
library(ggbeeswarm)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

celltype <- "Microglia"

in_dir <- paste0("results/wgcna_consensus/", celltype, "/")
data_dir <- paste0("data/19_subclustering3/", str_to_lower(celltype), "/")

scores <- read.csv(paste0(in_dir, "module_scores_ucell.csv"))

scores <- scores %>% 
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")),
         tissue = factor(tissue,
                         levels = c("Motor cortex", "Cervical spinal cord")))

pb <- scores %>% 
  group_by(orig.ident, group, tissue) %>% 
  mutate(across(where(is.numeric), median)) %>% 
  distinct(orig.ident,
           .keep_all = T)

modules <- read.csv(paste0(in_dir, "modules.csv"))

mois <- c("blue", "turquoise", "cyan", "green")

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
                     levels = c("Motor cortex", "Cervical spinal cord"))

obj <- AddMetaData(obj,
                   scores)

# kME plots ---------------------------------------------------------------

for (i in seq_along(mois)){

modules %>% 
  filter(color == mois[i]) %>% 
  arrange(desc(!!sym(paste0("kME_", mois[i])))) %>% 
  mutate(gene_name = fct_inorder(gene_name)) %>% 
  slice_head(n = 30) %>%
  ggplot(aes(x = !!sym(paste0("kME_", mois[i])),
             y = gene_name)) + 
  geom_col(color = "black",
           fill = mois[i]) + 
  scale_y_discrete(limits = rev) + 
  scale_x_continuous(expand = c(0, 0)) + 
  labs(x = "kME") +
  theme_linedraw(base_size = 12) + 
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        # axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(in_dir, mois[i], "_kme_bars.png"),
       units = "in", dpi = 600,
       height = 6, width = 2.5)

# Expression plots (sc) --------------------------------------------------------

scores %>% 
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) + 
  geom_violin(aes(fill = group)) + 
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
ggsave(filename = paste0(in_dir, mois[i], "_expression_sc.png"),
       units = "in", dpi = 600,
       height = 3, width = 4)

# Expression plots (pb) ---------------------------------------------------

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

# FeaturePlots ------------------------------------------------------------

for (i in mois){
  
  FeaturePlot_scCustom(obj,
                       features = paste0(i, "_UCell_kNN"),
                       colors_use = viridis_inferno_light_high) + 
    ggtitle(paste0(str_to_title(i), " module expression")) + 
    theme(plot.title = element_text(hjust = 0.5, face = "plain"))
  ggsave(filename = paste0(in_dir, i, "_expression_umap.png"),
         units = "in", dpi = 600,
         height = 4, width = 4.5)
}
# ORA result plots --------------------------------------------------------

# df <- read.csv(paste0(in_dir, mois[i], "_gsea.csv"))
# 
# df %>% 
#   filter(Count > 10) %>% 
#   arrange(p.adjust) %>% 
#   slice_head(n = 20) %>% 
#   
