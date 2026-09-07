library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

in_dir <- "results/speckle/"

plots_dir <- "figures/"

# Stat figures ------------------------------------------------------------

stats <- list.files(in_dir,
                    pattern = "stats")
tissues <- str_split_i(stats, "_", i = 1)

stats <- map(paste0(in_dir, stats),
             read.csv)

for (i in seq_along(stats)) {
  stats[[i]]$tissue <- tissues[i]
}

stats <- list_rbind(stats)

stats <- stats %>% 
  mutate(sals_fc = log2(PropMean.sALS / PropMean.Control),
         c9_fc = log2(PropMean.C9orf72 / PropMean.Control),
         sig = if_else(FDR < 0.05, "sig", "ns") %>% 
           factor(levels = c("sig", "ns"),
                  labels = c("True", "False")),
         tissue = factor(tissue,
                         levels = c("brain", "sc", "muscle"),
                         labels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")))
  
sig_cells <- stats %>% 
  filter(FDR < 0.05)  

stats  %>% 
  mutate(label = if_else(FDR < 0.05, X, "")) %>%
  ggplot(aes(x = sals_fc,
             y = c9_fc)) + 
  facet_wrap(. ~ tissue,
             nrow = 1) + 
  geom_point(aes(fill = sig,
                 size = BaselineProp),
             shape = 21) + 
  geom_label_repel(aes(label = label, color = sig),
                   show.legend = FALSE,
                   seed = 42,                    # reproducible layout
                   force = 5,
                   force_pull = 0.5,
                   max.overlaps = Inf,            # never drop a label
                   box.padding = 0.6,             # space around label
                   point.padding = 0.4,           # space around the point it's labeling
                   min.segment.length = 0,        # always draw the leader line
                   segment.color = "grey40",
                   # THIS is the part that avoids all points, not just labeled ones:
                   max.time = 5,
                   max.iter = 100000) +
  scale_color_manual(values = c("dodgerblue", "grey")) +
  scale_fill_manual(values = c("dodgerblue", "grey")) + 
  labs(fill = "Adjusted p < 0.05",
       size = "Average proportion of\ntotal tissue's cells",
       y = "log2FC in C9orf72-ALS vs Control",
       x = "log2FC in sALS vs Control") + 
  guides(fill = guide_legend(ncol = 1,
                             override.aes = list(size = 4))) + 
  theme_linedraw(base_size = 12) + 
  theme(strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(plots_dir, "speckle_volcanoes.png"),
       units = "in", dpi = 600,
       height = 3, width = 10)

# Dot plots ------------------------------------------------------------

props <- list.files(in_dir, pattern = "frequencies")
props <- map(paste0(in_dir, props),
             read.csv)
names(props) <- tissues

sig_cells_guide <- stats %>% 
  filter(FDR < 0.05) %>%
  dplyr::select(X, tissue) %>% 
  mutate(tissue2 = case_when(tissue == "Motor cortex" ~ "brain",
                             tissue == "Cervical spinal cord" ~ "sc",
                             tissue == "Skeletal muscle" ~ "muscle"))

for (i in 1:nrow(sig_cells_guide)) {
  tissue <- sig_cells_guide$tissue2[i]
  cell <- sig_cells_guide$X[i]
  cell <- str_replace_all(cell, " ", ".")
  
  props[[tissue]] %>%
    mutate(group = factor(group,
                          levels = c("Control", "sALS", "C9orf72"),
                          labels = c("Control", "sALS", "C9orf72-ALS"))) %>% 
    ggplot(aes(x = group,
               y = !!sym(cell))) + 
    geom_quasirandom(aes(fill = group),
                     show.legend = F,
                     shape = 21,
                     size = 4,
                     alpha = 0.7) + 
    stat_summary(fun = mean,
                 geom = "crossbar") +
    stat_summary(fun.data = mean_se,
                 geom = "errorbar",
                 linewidth = 1.2,
                 width = 0.6) +
    scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
    labs(y = "% of all cells") +
    ggtitle(paste0(str_replace_all(cell, "[.]", " "), " in\n", str_to_lower(sig_cells_guide$tissue[i]))) +
    theme_linedraw(base_size = 12) + 
    theme(strip.text = element_text(face = "bold", color = "black"),
          strip.background = element_rect(fill = "gray", color = "black"),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5))
  ggsave(filename = paste0(plots_dir, "speckle_dots_", cell, "_", tissue, ".png"),
         units = "in", dpi = 600,
         height = 3, width = 3)
}
