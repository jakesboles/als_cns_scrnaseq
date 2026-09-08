suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(BPCells)
  library(dittoSeq)
})

options(future.globals.maxSize = 250 * 1024^3)

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

data_dir <- "data/19_subclustering3/"

results_dir <- "results/19_subclustering3/"

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

DimPlot_scCustom(obj,
                 group.by = "cell_type3",
                 raster.dpi = c(900, 900), # comment out if running smaller object
                 pt.size = 3, # comment out if running smaller object
                 colors_use = paletteer_d("ggsci::default_locuszoom")) + 
  ggtitle("Muscle fiber subtype") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4)))
ggsave(filename = paste0(results_dir, target, "/umap_celltype.png"),
       units = "in", dpi = 600,
       height = 7, width = 10)

DimPlot_scCustom(obj,
                 group.by = "tissue",
                 colors_use = JCO_Four()) + 
  ggtitle("Tissue") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4)))
ggsave(filename = paste0(results_dir, target, "/umap_tissue.png"),
       units = "in", dpi = 600,
       height = 7, width = 10)


markers <- read.csv(paste0(results_dir, target, "/fam_output.csv"))

top <- markers %>% 
  # mutate(cluster = if_else(cluster %in% c("MN", "SN", "IN", "EN"), paste0("Spinal ", cluster), paste0("Cortical ", cluster))) %>%
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 12,
                           make_unique = T,
                           named_vector = F)


dittoDotPlot(obj,
             vars = top,
             group.by = "cell_type3") + 
  scale_y_discrete(limits = rev) + 
  labs(y = "Subtype",
       size = "pct",
       color = "exp")
ggsave(filename = paste0(results_dir, target, "/marker_dotplot.png"),
       units = "in", dpi = 600,
       height = 4, width = 15)
