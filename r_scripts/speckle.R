library(speckle)
library(Seurat)
library(ggplot2)
library(scCustomize)
library(ggbeeswarm)
library(tidyverse)
library(paletteer)

setwd("/projects/b1169/boles/als_cns_scrnaseq/")

data_in_dir <- "data/17_obj_reassembly/"

results_dir <- "results/speckle/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

tissues <- list.dirs(data_in_dir,
                    recursive = F,
                    full.names = F)

for (i in seq_along(tissues)) {
  df <- readRDS(paste0(data_in_dir, tissues[i], "/metadata.rds"))
  
  stats <- propeller(clusters = df$cell_type3,
                     sample = df$orig.ident,
                     group = df$group,
                     transform = "logit")
  
  write.csv(stats,
            paste0(results_dir, tissues[i], "_propeller_stats.csv"),
            row.names = T)
  
  props <- df %>%
    group_by(orig.ident) %>%
    mutate(total_cells = n()) %>%
    ungroup() %>%
    group_by(orig.ident, cell_type3) %>%
    mutate(fraction = (n() / total_cells) * 100) %>%
    dplyr::select(c(orig.ident, group, cell_type3, fraction)) %>%
    distinct(orig.ident, cell_type3, .keep_all = T) %>%
    pivot_wider(names_from = "cell_type3",
                values_from = "fraction")
  
  # props[is.na(props)] <- 0
  
  write.csv(props,
            paste0(results_dir, tissues[i], "_cell_frequencies.csv"),
            row.names = F)
  
}
