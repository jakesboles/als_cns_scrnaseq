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

target <- "neurons"

meta <- readRDS(paste0(data_dir, target, "/metadata.rds"))
umap <- readRDS(paste0(data_dir, target, "/harmony_umap.rds"))

mat <- open_matrix_dir(paste0(data_dir, target, "/bpcells_data"))

obj <- CreateSeuratObject(counts = mat,
                          meta.data = meta)
obj[["RNA"]]$data <- mat
obj[["umap"]] <- umap

obj$tissue <- factor(obj$tissue, levels = c("Motor cortex", "Cervical spinal cord"))

obj@meta.data <- obj@meta.data %>%
mutate(cell_type4 = if_else(cell_type3 %in% c("MN", "SN", "IN", "EN"), paste0("Spinal ", cell_type3), paste0("Cortical ", cell_type3)))

umap_subtypes <- DimPlot_scCustom(obj,
                                  group.by = "cell_type4",
                                  raster = T,
                                  raster.dpi = c(1200, 1200), 
                                  # colors_use = , # neuron color palette
                                  pt.size = 4) + 
  ggtitle("Neuron subtype") + 
  guides(color = guide_legend(ncol = 1, # ncol = 1 for MFs
                              override.aes = list(size = 2))) + 
  umap_labs + 
  theme_als_project(title_size = 7) + 
  umap_theme + 
  theme(legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(umap_subtypes,
            filename = paste0(plots_dir, target, "_subtype_umap.pdf"),
            height = 3, width = 3)

umap_tissue <- DimPlot_scCustom(obj,
                                group.by = "tissue",
                                raster = T,
                                raster.dpi = c(1200, 1200), 
                                pt.size = 4,
                                colors_use = JCO_Four()) + 
  ggtitle("Tissue") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 2))) + 
  umap_labs + 
  theme_als_project(title_size = 7) + 
  umap_theme + 
  theme(legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(umap_tissue,
            filename = paste0(plots_dir, target, "_tissue.pdf"),
            height = 3, width = 3)

# Marker dot plot ---------------------------------------------------------

markers <- read.csv(paste0(results_dir, target, "/fam_output.csv"))

top <- markers %>% 
  mutate(cluster = if_else(cluster %in% c("MN", "SN", "IN", "EN"), paste0("Spinal ", cluster), paste0("Cortical ", cluster))) %>%
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 3,
                      make_unique = T,
                      named_vector = F)

dittoDotPlot(obj,
             vars = top,
             group.by = "cell_type4", 
             legend.color.title = "Relative\nexpression",
             legend.size.title = "Percent\nexpression") + 
  scale_y_discrete(limits = rev) + 
  ggtitle("Neuron subtype markers") +
  scale_size(range = c(1, 4),
             name = "Percent\nexpression") +
  scale_color_gradient2(name = "Relative\nexpression") +
  theme_als_project(title_size = 7) + 
  rotate_x_text(90) +
  theme(axis.title.y = element_blank(),
        legend.box.spacing = unit(0.1, "cm"),
        legend.spacing.y = unit(0.01, "cm"),
        legend.title = element_text(margin = margin(b = 0.05, unit = "cm")))
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_markers.pdf"),
            height = 2.75, width = 6.5)

# Stacked bars by group ---------------------------------------------------

obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

sc_cells <- obj@meta.data %>%
  filter(tissue == "Cervical spinal cord") %>%
  rownames()

br_cells <- obj@meta.data %>%
  filter(tissue == "Motor cortex") %>%
  rownames()

neuron_pal <- DiscretePalette_scCustomize(num_colors = length(unique(obj$cell_type4)),
                                          "polychrome")

br_pal <- neuron_pal[1:length(unique(obj$cell_type4[str_detect(obj$cell_type4, "Cortical")]))]

sc_pal <- setdiff(neuron_pal, br_pal)

cortical_neuron_bars <- dittoBarPlot(obj,
                                     cells.use = br_cells,
                                     var = "cell_type4",
                                     group.by = "group",
                                     main = "Cortical neuron\nsubtypes by ALS group",
                                     color.panel = br_pal,
                                     x.reorder = c(2, 3, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of cortical neurons") +
  theme_als_project(title_size = 7) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(cortical_neuron_bars,
            filename = paste0(plots_dir, target, "_cortical_subtypes_bars.pdf"),
            height = 2, width = 1.5)

sc_neuron_bars <- dittoBarPlot(obj,
                               cells.use = sc_cells,
                               var = "cell_type4",
                               group.by = "group",
                               main = "Spinal neuron\nsubtypes by ALS group",
                               color.panel = sc_pal,
                               x.reorder = c(2, 3, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of spinal neurons") +
  theme_als_project(title_size = 7) +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
save_figure(sc_neuron_bars,
            filename = paste0(plots_dir, target, "_spinal_subtypes_bars.pdf"),
            height = 2, width = 1.5)

cortical_neuron_bars + sc_neuron_bars +
  plot_layout(nrow = 1)
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_subtypes_bars.pdf"),
            height = 5, width = 10)