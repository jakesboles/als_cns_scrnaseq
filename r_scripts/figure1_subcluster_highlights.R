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

# uncomment and use cell_type4 for plotting if looking at neurons
# obj$tissue <- factor(obj$tissue, levels = c("Motor cortex", "Cervical spinal cord"))

# obj@meta.data <- obj@meta.data %>%
  # mutate(cell_type4 = if_else(cell_type3 %in% c("MN", "SN", "IN", "EN"), paste0("Spinal ", cell_type3), paste0("Cortical ", cell_type3)))

umap_subtypes <- DimPlot_scCustom(obj,
                 group.by = "cell_type3",
                 raster = T,
                 raster.dpi = c(900, 900), 
                 # colors_use = , # neuron color palette
                 colors_use = paletteer_d("ggthemes::hc_default"), # MF color palette
                 pt.size = 3) + 
  ggtitle("Muscle fiber subtype") + 
  guides(color = guide_legend(ncol = 1, # ncol = 1 for MFs
                              override.aes = list(size = 4))) + 
  umap_labs + 
  theme_als_project() + 
  umap_theme
save_figure(umap_subtypes,
            filename = paste0(plots_dir, target, "_subtype_umap.pdf"),
            height = 4, width = 5)
# 4 x 5 for MFs
# ggsave(filename = paste0(results_dir, target, "/umap_celltype.png"),
       # units = "in", dpi = 600,
       # height = 7, width = 10)

# no need to run this for muscle fibers
umap_tissue <- DimPlot_scCustom(obj,
                 group.by = "tissue",
                 raster = T,
                 raster.dpi = c(900, 900), 
                 pt.size = 3,
                 colors_use = JCO_Four()) + 
  ggtitle("Tissue") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4))) + 
  umap_labs + 
  theme_als_project() + 
  umap_theme
save_figure(umap_tissue,
            filename = paste0(plots_dir, target, "_tissue.pdf"),
            height = 4, width = 5)
# ggsave(filename = paste0(results_dir, target, "/umap_tissue.png"),
       # units = "in", dpi = 600,
       # height = 7, width = 10)

umap_subtypes + umap_tissue
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_umap_patchwork.pdf"),
            height = 4, width = 13)

# Marker dot plot ---------------------------------------------------------

markers <- read.csv(paste0(results_dir, target, "/fam_output.csv"))

top <- markers %>% 
  # mutate(cluster = if_else(cluster %in% c("MN", "SN", "IN", "EN"), paste0("Spinal ", cluster), paste0("Cortical ", cluster))) %>%
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 12, # 12 for MFs
                      make_unique = T,
                      named_vector = F)

dittoDotPlot(obj,
             vars = top,
             group.by = "cell_type3", # make sure this is appropriate for the cell type
             legend.color.title = "Relative\nexpression",
             legend.size.title = "Percent\nexpression") + 
  scale_y_discrete(limits = rev) + 
  labs(y = "Subtype") +
  theme_als_project() + 
  theme(axis.title.y = element_blank())
save_figure(last_plot(),
            filename = paste0(plots_dir, target, "_markers.pdf"),
            height = 4, width = 14)
# 4 x 14 for muscle fibers
# ggsave(filename = paste0(results_dir, target, "/marker_dotplot.png"),
       # units = "in", dpi = 600,
       # height = 4, width = 15)

# Stacked bars by group ---------------------------------------------------

# uncomment code below for neurons 

obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

# sc_cells <- obj@meta.data %>% 
#   filter(tissue == "Cervical spinal cord") %>% 
#   rownames()
# 
# br_cells <- obj@meta.data %>% 
#   filter(tissue == "Motor cortex") %>% 
#   rownames()
# 
# cortical_neuron_bars <- dittoBarPlot(obj,
#              var = "cell_type4",
#              group.by = "group",
#              cells.use = br_cells,
#              main = "Cortical neurons",
#              x.reorder = c(2, 3, 1)) + 
#   scale_y_continuous(expand = c(0, 0)) +
#   labs(y = "Proportion of neurons") +
#   theme_als_project() + 
#   theme(axis.title.x = element_blank())
# save_figure(cortical_neuron_bars,
#             filename = paste0(plots_dir, target, "_cortical_subtypes_bars.pdf"),
#             height = 6, width = 6)
# 
# sc_neuron_bars <- dittoBarPlot(obj,
#                                var = "cell_type4",
#                                group.by = "group",
#                                cells.use = sc_cells,
#                                main = "Spinal neurons",
#                                x.reorder = c(2, 3, 1)) + 
#   scale_y_continuous(expand = c(0, 0)) + 
#   labs(y = "Proportion of neurons") + 
#   theme_als_project() + 
#   theme(axis.title.x = element_blank())
# save_figure(sc_neuron_bars,
#             filename = paste0(plots_dir, target, "_spinal_subtypes_bars.pdf"),
#             height = 6, width = 6)
# 
# cortical_neuron_bars + sc_neuron_bars + 
#   plot_layout(nrow = 1)
# save_figure(last_plot(),
#             filename = paste0(plots_dir, target, "_subtypes_bars.pdf"),
#             height = 5, width = 10)

dittoBarPlot(obj,
             var = "cell_type3",
             group.by = "group",
             main = "Muscle fibers",
             color.panel = paletteer_d("ggthemes::hc_default"),
             x.reorder = c(2, 3, 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion of muscle fibers") +
  theme_als_project() +
  theme(axis.title.x = element_blank())
save_figure(sc_neuron_bars,
            filename = paste0(plots_dir, target, "_subtypes_bars.pdf"),
            height = 6, width = 6)
