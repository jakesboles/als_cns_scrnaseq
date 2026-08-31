library(tidyverse)
library(ggrepel)
library(scCustomize)
library(UpSetR)
library(ggplot2)
library(reshape2)
library(paletteer)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

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
  mutate(order = paste0(celltype, "_", tissue)) %>%
  arrange(desc(total)) %>% 
  mutate(order = fct_inorder(order)) %>% 
  pivot_longer(!c(celltype, order, tissue, total),
               values_to = "n",
               names_to = "type") %>%
  mutate(tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Spinal cord", "Skeletal muscle")),
         type = factor(type, 
                       levels = c("sals_only", "c9_only", "shared"),
                       labels = c("sALS only", "C9-ALS only", "Shared"))) %>%
  ggplot(aes(x = order,
             y = n)) + 
  geom_bar(aes(fill = type),
           stat = "identity",
           color = "black",
           linewidth = 0.4) + 
  facet_wrap(. ~ tissue, ncol = 3,
             scales = "free_x") +
  scale_fill_manual(values = c("chartreuse3", "magenta3", "darkslategrey")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    labels = function(x) str_split_i(x, "_", i = 1)    # remove everything through the last “_”
  ) +
  ylab("# DEGs") + 
  labs(fill = "DEGs versus\ncontrol in:") +
  theme_linedraw(base_size = 16) +
  theme(axis.title.x = element_blank(),
        # axis.title.y = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 14),
        # axis.text.x = element_text(angle = , hjust = 1, vjust = 1),
        # axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 20, color = "black"))
