library(ggplot2)
library(tidyverse)
library(patchwork)
library(Seurat)
library(scCustomize)
library(BPCells)
library(paletteer)
library(ggrastr)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir0 <- "figures/"

source("r_scripts/plot_theme.R")

# 3-way barchart of DEGs --------------------------------------------------

results_dir <- "results/deseq2"

get_degs <- function(df,
                     pval_cutoff,
                     logfc_cutoff){
  genes <- df %>% 
    filter(padj < pval_cutoff & 
             abs(log2FoldChange) > logfc_cutoff) %>% 
    pull(X)
  
  return(genes)
}

p_thresh <- 0.05
lfc_thresh <- log2(1.5)

degs <- list()

for (i in c("brain", "sc", "muscle")) {
  celltypes <- list.dirs(paste0(results_dir, "/", i),
                         full.names = F,
                         recursive = F)
  
  # check that all output is present, get rid of cells whose job failed
  
  celltypes <- celltypes[file.exists(paste0(results_dir, "/", i, "/", celltypes, "/sALS_vs_Control.csv"))]
  celltypes <- celltypes[file.exists(paste0(results_dir, "/", i, "/", celltypes, "/C9orf72_vs_Control.csv"))]
  
  degs[[i]] <- matrix(nrow = length(celltypes),
                      ncol = 5)
  
  for (j in seq_along(celltypes)){
    sals <- read.csv(paste0(results_dir, "/", i, "/", celltypes[j], "/sALS_vs_Control.csv"))
    
    c9 <- read.csv(paste0(results_dir, "/", i, "/", celltypes[j], "/C9orf72_vs_Control.csv"))
    
    sals_degs <- get_degs(sals,
                          pval_cutoff = p_thresh,
                          logfc_cutoff = lfc_thresh)
    
    c9_degs <- get_degs(c9,
                        pval_cutoff = p_thresh,
                        logfc_cutoff = lfc_thresh)
    
    shared <- intersect(sals_degs,
                        c9_degs)
    
    c9_only <- setdiff(c9_degs,
                       sals_degs)
    
    sals_only <- setdiff(sals_degs,
                         c9_degs)
    
    degs[[i]][j,] <- c(celltypes[j], length(shared), length(c9_only), length(sals_only), length(union(sals_degs, c9_degs)))    
  }
  
  colnames(degs[[i]]) <- c("celltype", "shared", "c9_only", "sals_only", "total")
  
  degs[[i]] <- degs[[i]] %>%
    as.data.frame() %>% 
    mutate(tissue = i)
}

df <- degs %>% 
  list_rbind() 

df <- df %>% 
  mutate_at(c("shared", "c9_only", "sals_only", "total"),
            as.numeric)

df %>% 
  filter(total != 0) %>%
  mutate(celltype = str_replace_all(celltype, "_", " ")) %>%
  mutate(celltype = case_when(celltype == "Endothelial cell" ~ "Endo. cell",
                              celltype == "Proliferating myeloid cell" ~ "Prolif. myeloid cell",
                              celltype == "Smooth muscle cell" ~ "SMC",
                              .default = celltype)) %>%
  mutate(order = paste0(celltype, "_", tissue)) %>%
  arrange(desc(total)) %>% 
  mutate(order = fct_inorder(order)) %>% 
  pivot_longer(!c(celltype, order, tissue, total),
               values_to = "n",
               names_to = "type") %>%
  mutate(tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")),
         type = factor(type, 
                       levels = c("sals_only", "c9_only", "shared"),
                       labels = c("sALS only", "C9orf72-ALS only", "Both"))) %>%
  ggplot(aes(x = order,
             y = n)) + 
  ggrastr::rasterize(
    geom_bar(aes(fill = type),
           stat = "identity",
           color = "black",
           linewidth = 0.25),
    dpi = 600) + 
  facet_wrap(. ~ tissue, ncol = 3,
             scales = "free_x") +
  scale_fill_manual(values = c("#0CAA00", "#CC00FF", "darkslategrey")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    labels = function(x) str_split_i(x, "_", i = 1)    # remove everything through the last “_”
  ) +
  guides(fill = guide_legend(direction = "vertical",
                             override.aes = list(color = NA))) +
  ylab("# DEGs") + 
  labs(fill = "DEGs versus\ncontrol in:") +
  theme_linedraw() +
  theme_als_project() + 
  rotate_x_text(90) + 
  theme(axis.title.x = element_blank(),
        legend.box.spacing = unit(0.05, "cm"),
        legend.title = element_text(margin = margin(b = 0.05, unit = "cm")),
        strip.text = element_text(size = 7))
save_figure(last_plot(),
            filename = paste0(plots_dir0, "figure1/deseq2_bars.pdf"),
            height = 2.2, width = 6.5)

# Rotated axis barchart? --------------------------------------------------

df %>% 
  filter(total != 0) %>%
  mutate(celltype = str_replace_all(celltype, "_", " ")) %>%
  mutate(celltype = case_when(celltype == "Endothelial cell" ~ "Endo. cell",
                              celltype == "Proliferating myeloid cell" ~ "Prolif. myeloid cell",
                              celltype == "Smooth muscle cell" ~ "SMC",
                              .default = celltype)) %>%
  mutate(order = paste0(celltype, "_", tissue)) %>%
  arrange(desc(total)) %>% 
  mutate(order = fct_inorder(order)) %>% 
  pivot_longer(!c(celltype, order, tissue, total),
               values_to = "n",
               names_to = "type") %>%
  mutate(tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")),
         type = factor(type, 
                       levels = c("sals_only", "c9_only", "shared"),
                       labels = c("sALS only", "C9orf72-ALS only", "Both"))) %>%
  ggplot(aes(y = order,
             x = n)) + 
  ggrastr::rasterize(
    geom_bar(aes(fill = type),
             stat = "identity",
             color = "black",
             linewidth = 0.25),
    dpi = 600) + 
  facet_wrap(. ~ tissue, ncol = 3,
             scales = "free_y") +
  scale_fill_manual(values = c("#0CAA00", "#CC00FF", "darkslategrey")) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev,
    labels = function(x) str_split_i(x, "_", i = 1)    # remove everything through the last “_”
  ) +
  guides(fill = guide_legend(direction = "horizontal",
                             override.aes = list(color = NA))) +
  xlab("# DEGs") + 
  labs(fill = "DEGs versus control in:") +
  ggtitle("Differential expression in the ALS motor circuit") +
  theme_linedraw() +
  theme_als_project() + 
  # rotate_x_text(90) + 
  theme(axis.title.y = element_blank(),
        legend.box.spacing = unit(0.01, "cm"),
        # legend.title = element_text(margin = margin(b = 0.05, unit = "cm")),
        # strip.text = element_text(size = 7),
        legend.position = "bottom",
        panel.spacing = unit(0.05, "cm"))
save_figure(last_plot(),
            filename = paste0(plots_dir0, "figure1/deseq2_bars_flipped.pdf"),
            height = 2.75, width = 6.5)

# 3-way dotplot of DEGs? --------------------------------------------------

head(df)

df %>%
  mutate(c9 = shared + c9_only,
         sals = shared + sals_only,
         .keep = "unused") %>%
  filter(total != 0) %>%
  mutate(celltype = str_replace_all(celltype, "_", " ")) %>%
  mutate(order = paste0(celltype, "_", tissue)) %>%
  arrange(desc(total)) %>% 
  mutate(order = fct_inorder(order)) %>% 
  pivot_longer(!c(celltype, order, tissue, total),
               values_to = "n",
               names_to = "type") %>% 
  mutate(tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")),
         type = factor(type, 
                       levels = c("sals", "c9"),
                       labels = c("sALS", "C9orf72-ALS"))) %>%
  ggplot(aes(x = order,
             y = type)) + 
  geom_tile(aes(fill = n),
            color = "black") +
  facet_wrap(. ~ tissue,
             nrow = 1,
             scales = "free_x") + 
  scale_x_discrete(
    labels = function(x) str_split_i(x, "_", i = 1)    # remove everything through the last “_”
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_paletteer_c("ggthemes::Sunset-Sunrise Diverging") +
  labs(fill = "# DEGs\nversus control") +
  theme_bw() +
  theme_als_project() + 
  rotate_x_text() + 
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.box.spacing = unit(0.1, "cm"))
