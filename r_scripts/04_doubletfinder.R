# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library("ggthemes")
  library("ggrepel")
  library("grid")
  library("DoubletFinder")
  library("doMC")
  library("xlsx")
  library("RColorBrewer")
  library(scCustomize)
  library(ggpubr)
  # library(BPCells)
  library(patchwork)
})

# Script again adapted from Anne's projects

proj_dir <- "/projects/b1169/boles/als_multitissue_scfrp/"

plots_dir <- paste0(proj_dir, "plots/04_doubletfinder/")
dir.create(plots_dir,
           showWarnings = F, recursive = T)

data_out_dir <- paste0(proj_dir, "data/04_doubletfinder/")
dir.create(data_out_dir,
           showWarnings = F, recursive = T)

obj <- readRDS(paste0(proj_dir, "data/03_qc2/filtered_obj.rds"))

# Load pre-filter stats to get counts for each sample
# Concatenate stats and order according to the way object layers will be split for DF
# This way, each sample gets a customized doublet formation rate depending on the # of cells detected 

brain_counts <- read.csv("/projects/b1169/boles/als_multitissue_scfrp/tab_data/02_qc1/median_stats_brain.csv")[-31, ] %>%
  mutate(id = paste0(id, "_b"))

sc_counts <- read.csv("/projects/b1169/boles/als_multitissue_scfrp/tab_data/02_qc1/median_stats_sc.csv")[-31, ] %>%
  mutate(id = paste0(id, "_s"))

muscle_counts <- read.csv("/projects/b1169/boles/als_multitissue_scfrp/tab_data/02_qc1/median_stats_muscle.csv")[-31, ] %>%
  mutate(id = paste0(id, "_m"))

counts <- rbind(brain_counts, sc_counts) %>%
  rbind(muscle_counts)

# Define function to run DoubletFinder
run_doubletfinder <- function(s, doublet_rate) {
  
  # Standard normalization and scaling 
  s <- s %>% NormalizeData() %>% FindVariableFeatures() %>% ScaleData()
  
  # Default Seurat clustering and TSNE 
  s <- s %>% RunPCA() %>% FindNeighbors(dims = 1:10) %>% FindClusters()
  s <- RunUMAP(s, dims = 1:10) # had to set check_duplicates to FALSE for some samples
  
  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(s, PCs = 1:10, sct = FALSE) 
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  max_index <- which.max(bcmvn$BCmetric)
  optimal_pK <- as.numeric(as.character(bcmvn[max_index, "pK"]))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- s@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate*nrow(s@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  # Run DoubletFinder without homotypic adjustment 
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.unadj"
  
  # Run DoubletFinder with homotypic adjustment
  pANN <- colnames(s@meta.data)[grep("^pANN", colnames(s@meta.data))]
  print(pANN)
  s <- doubletFinder(s, PCs = 1:10, pN = 0.25, pK = optimal_pK, nExp = nExp_poi.adj, reuse.pANN = pANN, sct = FALSE)
  
  # Rename meta column
  colnames(s@meta.data)[grep("DF.classifications*", colnames(s@meta.data))] <- "DF.adj"
  
  # Plot unadjusted vs. adjusted doublets in TSNE coordinates
  plt <- ggarrange(
    (DimPlot_scCustom(s, reduction = "umap", group.by = "DF.unadj", 
                      pt.size = 1, shuffle = TRUE, alpha = 0.6) + 
       ggtitle("Unadjusted")),
    (DimPlot_scCustom(s, reduction = "umap", group.by = "DF.adj", 
                      pt.size = 1, shuffle = TRUE, alpha = 0.6) +
       ggtitle("Adjusted for Homotypic Proportion")),
    ncol = 2, nrow = 1,
    common.legend = T, legend = "right"
  )
  
  ggsave(plt,
         filename = paste0(plots_dir, unique(s$orig.ident), "_umap.png"),
         units = "in", dpi = 600,
         height = 5, width = 10)
  
  return(s)
}

# Split into sample objects
obj_list <- SplitObject(obj, split.by = "orig.ident")

counts$id %in% names(obj_list)

idx <- match(counts$id, names(obj_list))
counts <- counts[idx, ] %>%
  pull(Cell_count)

doublet_rates <- (counts / 10000) * 0.08

# Run DoubletFinder on all samples
for (s in seq_along(obj_list)) {
  message(paste0(
    paste(rep("~", 40), collapse = ""),
    "Running ", names(obj_list)[s],
    paste(rep("~", 40), collapse = "")
  ))
  
  obj_list[[s]] <- run_doubletfinder(obj_list[[s]], doublet_rates[s])
  saveRDS(obj_list[[s]], 
          paste0(data_out_dir, names(obj_list)[s], ".rds"))
  
  # obj_list[[s]] <- DietSeurat(obj_list[[s]],
  #                             layers = "counts")
  # 
  # obj_list[[s]][["RNA"]]$data <- NULL
  # obj_list[[s]][["RNA"]]$scale.data <- NULL
}

# message("Merging object again") 
# 
# obj <- Merge_Seurat_List(obj_list, 
#                          add.cell.ids = names(obj_list))
# 
# message("Joining layers again")
# 
# obj <- JoinLayers(obj)
# 
# message("Saving full Seurat object as RDS (just in case)")
# 
# message("Saving BPCells matrix")
# 
# write_matrix_dir(mat = obj[["RNA"]]$counts, 
#                  dir = data_out_dir)
# 
# message("Saving metadata")
# 
# saveRDS(obj@meta.data, 
#         file = paste0(data_out_dir, "metadata.rds"))
# 
# message("Done!")