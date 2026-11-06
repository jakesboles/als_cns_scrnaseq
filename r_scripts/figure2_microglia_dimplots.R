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

plots_dir <- "figures/figure2/"

target <- "microglia"

meta <- readRDS(paste0(data_dir, target, "/metadata.rds"))
umap <- readRDS(paste0(data_dir, target, "/harmony_umap.rds"))

mat <- open_matrix_dir(paste0(data_dir, target, "/bpcells_data"))

obj <- CreateSeuratObject(counts = mat,
                          meta.data = meta)
obj[["RNA"]]$data <- mat
obj[["umap"]] <- umap

obj$tissue <- factor(obj$tissue, levels = c("Motor cortex", "Cervical spinal cord"),
                     labels = c("MCX", "SC"))
obj$group <- factor(obj$group,
                    levels = c("Control", "sALS", "C9orf72"),
                    labels = c("Control", "sALS", "C9orf72-ALS"))

microglia_umap <- function(grouping, title = NA, 
                           color_panel){
  p <- DimPlot_scCustom(obj,
                        group.by = grouping,
                        raster = T,
                        raster.dpi = c(1200, 1200), 
                        colors_use = color_panel,
                        pt.size = 4) + 
    guides(color = guide_legend(ncol = 1, # ncol = 1 for MFs
                                override.aes = list(size = 2))) + 
    umap_labs + 
    theme_als_project(title_size = 7) + 
    umap_theme + 
    theme(legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
          legend.box.spacing = unit(0.1, "cm"))
  
  if (!is.na(title)) p <- p + ggtitle(title)
  
  return(p)
}

microglia_umap("tissue", "Tissue", 
               JCO_Four())
save_figure(last_plot(),
            filename = paste0(plots_dir, "umap_tissue.pdf"),
            height = 1.5, width = 2)

microglia_umap("group", "ALS group",
               c("#b8b0a8", "#0CAA00", "#CC00FF"))
save_figure(last_plot(),
            filename = paste0(plots_dir, "umap_group.pdf"),
            height = 1.5, width = 2)

DimPlot_scCustom(obj,
                 group.by = "tissue",
                 raster = T,
                 raster.dpi = c(1200, 1200), 
                 colors_use = JCO_Four(),
                 pt.size = 4) + 
  ggtitle("Neuron subtype") +
  guides(color = guide_legend(ncol = 1, # ncol = 1 for MFs
                              override.aes = list(size = 2))) + 
  umap_labs + 
  theme_als_project(title_size = 7) + 
  umap_theme + 
  theme(legend.text = element_text(margin = margin(r = 0.01, unit = "cm")),
        legend.box.spacing = unit(0.1, "cm"))
