library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggrastr)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir0 <- "figures/"

source("r_scripts/plot_theme.R")

# QC metrics --------------------------------------------------------------

df <- readRDS("data/02_qc1/metadata.rds")

str(df)

df <- df %>%
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")))

qc_theme <- theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  legend.title = element_blank(),
  axis.title.x = element_blank(),
  )

p_umi <- df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = nCount_RNA)) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  ggrastr::rasterize(
    geom_violin(aes(fill = group),
                linewidth = 0.25),
    dpi = 600) +
  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  labs(y = "# UMIs per cell",
       x = "Donor") +
  ggtitle("Counts per cell") +
  scale_y_log10() + 
  theme_linedraw() +
  theme_als_project(title_size = 7) + 
  qc_theme
# ggsave(filename = paste0(plots_dir, "ncount.png"),
#        units = "in", dpi = 600,
#        height = 8, width = 12)

p_feature <- df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = nFeature_RNA)) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  ggrastr::rasterize(
    geom_violin(aes(fill = group),
                linewidth = 0.25),
    dpi = 600) +  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  labs(y = "# genes per cell",
       x = "Donor") +
  ggtitle("Unique genes per cell") +
  scale_y_log10() + 
  theme_linedraw() +
  theme_als_project(title_size = 7) + 
  qc_theme
# ggsave(filename = paste0(plots_dir, "nfeature.png"),
       # units = "in", dpi = 600,
       # height = 8, width = 12)

p_mito <- df %>% 
  arrange(group) %>%
  mutate(id = fct_inorder(id)) %>%
  ggplot(aes(x = id,
             y = percent_mito)) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  ggrastr::rasterize(
    geom_violin(aes(fill = group),
                linewidth = 0.25),
    dpi = 600) +  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  labs(y = "% mitochondrial\ngenes per cell",
       x = "Donor") +
  ggtitle("Mitochondrial genes per cell") +
  theme_linedraw() +
  theme_als_project(title_size = 7) + 
  qc_theme
# ggsave(filename = paste0(plots_dir, "mito.png"),
       # units = "in", dpi = 600,
       # height = 8, width = 12)

p_qc <- p_umi + p_feature + p_mito + 
  plot_layout(ncol = 1,
              guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box.spacing = unit(0.1, "cm"))

plots_dir <- paste0(plots_dir0, "supp1/")
dir.create(plots_dir, 
           recursive = T,
           showWarnings = F)

save_figure(p_qc,
            filename = paste0(plots_dir, "scrna_qc.pdf"),
            height = 5, width = 6.5)
