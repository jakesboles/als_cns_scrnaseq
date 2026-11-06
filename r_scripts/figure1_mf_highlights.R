# Interactive visualization of one 19_subclustering3.R target:
# cell-type/tissue UMAPs and a top-marker dot plot from its
# FindAllMarkers() output. Not a SLURM job.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(BPCells)
  library(dittoSeq)
  library(patchwork)
})

options(future.globals.maxSize = 250 * 1024^3)

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

source("r_scripts/plot_theme.R")

data_dir <- "data/19_subclustering3/"

results_dir <- "results/19_subclustering3/"

plots_dir <- "figures/supp2/"

target <- "muscle_fiber"

meta <- readRDS(paste0(data_dir, target, "/metadata.rds"))
umap <- readRDS(paste0(data_dir, target, "/harmony_umap.rds"))

mat <- open_matrix_dir(paste0(data_dir, target, "/bpcells_data"))

obj <- CreateSeuratObject(counts = mat,
                          meta.data = meta)
obj[["RNA"]]$data <- mat
obj[["umap"]] <- umap

umap_subtypes <- DimPlot_scCustom(obj,
                 group.by = "cell_type3",
                 raster = T,
                 raster.dpi = c(1200, 1200), 
                 # colors_use = , # neuron color palette
                 colors_use = paletteer_d("ggthemes::hc_default"), # MF color palette
                 pt.size = 4) + 
  ggtitle("Muscle fiber subtype") + 
  guides(color = guide_legend(ncol = 1, # ncol = 1 for MFs
                              override.aes = list(size = 2))) + 
  umap_labs + 
  theme_als_project(title_size = 7) + 
  umap_theme + 
  theme(legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(umap_subtypes,
            filename = paste0(plots_dir, target, "_subtype_umap.pdf"),
            height = 2, width = 3)

# Marker dot plot ---------------------------------------------------------

markers <- read.csv(paste0(results_dir, target, "/fam_output.csv"))

top <- markers %>% 
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 8, 
                      make_unique = T,
                      named_vector = F)

dittoDotPlot(obj,
             vars = top,
             group.by = "cell_type3", 
             legend.color.title = "Relative\nexpression",
             legend.size.title = "Percent\nexpression") + 
  scale_y_discrete(limits = rev) + 
  ggtitle("Muscle fiber subtype markers") +
  scale_size(range = c(1, 4),
             name = "Percent\nexpression") +
  scale_color_gradient2(name = "Relative\nexpression") +
  rotate_x_text(90) +
  theme_als_project(title_size = 7) + 
  theme(axis.title.y = element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.title = element_text(margin = margin(b = 0.05, unit = "cm")))
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_markers.pdf"),
            height = 2, width = 5)

# Stacked bars by group ---------------------------------------------------

obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

dittoBarPlot(obj,
             var = "cell_type3",
             group.by = "group",
             main = "Muscle fiber subtypes\nby ALS group",
             color.panel = paletteer_d("ggthemes::hc_default"),
             x.reorder = c(2, 3, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of muscle fibers") +
  theme_als_project(title_size = 7) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_subtypes_bars.pdf"),
            height = 2, width = 1.5)
