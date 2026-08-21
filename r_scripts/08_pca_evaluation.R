# Evaluates the sketch PCA from 08_sketch_pca.R -- elbow/loading/DimPlots to
# gauge batch/site/group structure, then a quick FindNeighbors()/RunUMAP()
# pass -- to help pick the number of PCs and judge whether integration is
# needed before committing to those choices for the full analysis. This
# script only reads 08's output and writes plots; it doesn't save any new
# matrices or objects.

# Load libraries
suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "plots/08_pca_evaluation/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

data_in_dir <- "data/07_norm_pca/"

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

# Load data + PCA per tissue ------------------------------------------

message2("Loading objects")

obj_list <- list()

for (i in seq_along(tissues$file)){
  tissue_file <- tissues$file[i]

  message(paste0("Loading ", tissue_file))

  data_mat <- open_matrix_dir(paste0(data_in_dir, tissue_file, "/bpcells_data"))
  meta <- readRDS(paste0(data_in_dir, tissue_file, "/metadata.rds"))
  pca <- readRDS(paste0(data_in_dir, tissue_file, "/pca.rds"))

  # Nothing here actually needs expression values -- ElbowPlot/PC loading
  # plots/PCA DimPlots all read straight from the attached `pca` reduction,
  # and FindNeighbors()/RunUMAP() below run on that reduction's embeddings
  # (via dims=), not on assay data. This is only reconstructed as a full
  # object because CreateSeuratObject() needs some counts matrix to
  # initialize, and it's the same sketch data 08 already wrote, so it's
  # assigned to both "counts" and "data" in case anything later (e.g. a
  # FeaturePlot) expects a normalized "data" layer.
  obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "sketch")
  obj[["RNA"]]$data <- data_mat
  obj[["pca"]] <- pca

  obj@meta.data <- obj@meta.data %>%
    mutate(site = case_when(str_detect(id, "AU") ~ "WashU",
                            str_detect(id, "GWF") ~ "Barrow",
                            str_detect(id, "GBB") ~ "Georgetown"))

  obj_list[[tissue_file]] <- obj
}

# PCA diagnostics -----------------------------------------------------------

for (i in seq_along(tissues$file)){
  tissue <- tissues$file[i]
  obj <- obj_list[[tissue]]

  message2(paste0("Making elbow plot for ", tissue))

  p <- ElbowPlot(obj,
                ndims = 100)
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_elbow_plot.png"),
         units = "in", dpi = 600,
         height = 6, width = 6,
         bg = "white")

  message2(paste0("Making PC loading plots for ", tissue))
  
  obj <- ScaleData(obj)

  Iterate_PC_Loading_Plots(obj,
                           file_path = plots_dir,
                           file_name = paste0(tissue, "_pca_loadings"))

  pca_plot2 <- function(dims){

    if (length(unique(Idents(obj))) == 3) {
      pal <- JCO_Four()
    } else {
      if (length(unique(Idents(obj))) == 90) {
        pal <- DiscretePalette_scCustomize(90,
                                           palette = "varibow",
                                           shuffle_pal = T)
      } else {
        pal <- DiscretePalette_scCustomize(num_colors = length(unique(Idents(obj))),
                                           palette = "glasbey")
      }
    }

    DimPlot_scCustom(obj,
                     reduction = "pca",
                     label = F,
                     dims = dims,
                     colors_use = pal) +
      theme(axis.text = element_text(size = 8),
            axis.title = element_text(size = 10),
            legend.justification = "left")
  }

  pca_grid <- function(group.by){

    design <- "
  AA######
  BBCC####
  DDEEFF##
  GGHHIIJJ
  "

    pca_plot2(c(1,2)) +
      pca_plot2(c(1,3)) + pca_plot2(c(2,3)) +
      pca_plot2(c(1,4)) + pca_plot2(c(2,4)) + pca_plot2(c(3,4)) +
      pca_plot2(c(1,5)) + pca_plot2(c(2,5)) + pca_plot2(c(3,5)) + pca_plot2(c(4,5)) +
      plot_layout(design = design,
                  guides = "collect")
  }

  message2(paste0("Making gridded PCA DimPlots for ", tissue))

  Idents(obj) <- "tissue"
  p <- pca_grid("tissue")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_tissue.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)

  Idents(obj) <- "id"
  p <- pca_grid("id")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_id.png"),
         units = "in", dpi = 600,
         height = 10, width = 15)

  Idents(obj) <- "Batch"
  p <- pca_grid("Batch")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_batch.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)

  Idents(obj) <- "Group"
  p <- pca_grid("Group")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_group.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)

  Idents(obj) <- "site"
  p <- pca_grid("site")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_pca_dimplot_site.png"),
         units = "in", dpi = 600,
         height = 10, width = 12)
  
  # message2(paste0("Making JackStraw plot for ", tissue))
  # 
  # js <- readRDS(paste0(data_in_dir, tissue, "/jackstraw_scores.rds"))
  # 
  # p <- js %>%
  #   ggplot(aes(x = PC, y = Score)) +
  #   geom_point(aes(color = Score < 0.05)) +
  #   geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
  #   scale_color_manual(values = c("gray60", "black")) +
  #   labs(y = "JackStraw p-value",
  #        color = "p < 0.05") +
  #   theme_bw()
  # 
  # ggsave(p,
  #        filename = paste0(plots_dir, tissue, "_jackstraw.png"),
  #        units = "in", dpi = 600,
  #        height = 5, width = 8)
  # 
  # suggested_pcs <- max(js$PC[js$Score < 0.05])
  # message(paste0("Suggested # PCs (largest with JackStraw p < 0.05) for ",
  #                tissue, ": ", suggested_pcs))
  
  message2(paste0("Finding PCA elbow for ", tissue))
  
  find_elbow <- function(stdev){
    n <- length(stdev)
    
    line_vec <- c(n - 1, stdev[n] - stdev[1])
    line_vec <- line_vec / sqrt(sum(line_vec^2))
    
    vecs <- cbind((1:n) - 1, stdev - stdev[1])
    proj_len <- vecs %*% line_vec
    proj <- proj_len %*% line_vec
    perp <- vecs - proj
    dist <- sqrt(rowSums(perp^2))
    
    which.max(dist)
  }
  
  # pca <- readRDS(paste0(data_in_dir, tissue, "/pca.rds"))
  elbow_pc <- find_elbow(Stdev(pca))
  message(paste0("Suggested # PCs (max-distance elbow) for ", tissue, ": ",
                 elbow_pc))
}

# Cluster cells and compute UMAP ---------------------------------------------

for (i in seq_along(tissues$file)){
  tissue <- tissues$file[i]
  obj <- obj_list[[tissue]]

  message2(paste0("Running FindNeighbors/RunUMAP for ", tissue))

  obj <- obj %>%
    FindNeighbors(dims = 1:30) %>%
    RunUMAP(dims = 1:30)

  p <- DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "tissue",
                   colors_use = JCO_Four())
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_umap_dimplot_tissue.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)

  p <- DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "batch",
                   colors_use = DiscretePalette_scCustomize(6,
                                                            palette = "ditto_seq"))
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_umap_dimplot_batch.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)

  p <- DimPlot_scCustom(obj,
                   reduction = "umap",
                   group.by = "group",
                   colors_use = JCO_Four())
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_umap_dimplot_group.png"),
         units = "in", dpi = 600,
         height = 5, width = 6)

  p <- DimPlot_scCustom(obj,
                   group.by = "id",
                   reduction = "umap")
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_umap_dimplot_id.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)

  p <- DimPlot_scCustom(obj,
                   group.by = "site",
                   reduction = "umap",
                   colors_use = JCO_Four())
  ggsave(p,
         filename = paste0(plots_dir, tissue, "_umap_dimplot_site.png"),
         units = "in", dpi = 600,
         height = 5, width = 8)
}
