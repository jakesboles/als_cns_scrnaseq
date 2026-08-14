# Load libraries
suppressMessages({
  library("plyr")
  library("tidyverse")
  library("Seurat")
  library(scCustomize)
  library(paletteer)
  library(stringr)
  library(janitor)
})

# Function to print clear log progress updates
message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Create directories ------------------------------------------------------
setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "plots/02_qc1/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

csv_dir <- "tab_data/02_qc1/"
dir.create(csv_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/02_qc1/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

message2("Reading in object")
obj <- readRDS("data/01_obj_creation/obj.rds")

# Add metadata ------------------------------------------------------------

message2("Organizing object for plotting")

meta <- read.csv("tab_data/metadata.csv")

meta <- meta %>%
  dplyr::select(-c(Barcode, Sample.Index, Dummy.code.)) %>%
  mutate(Tissue2 = case_when(Tissue == "Muscle" ~ "m",
                            Tissue == "Cervical spinal cord" ~ "s",
                            Tissue == "Motor cortex" ~ "b")) %>%
  unite(c(ID, Tissue2), 
        col = "orig.ident")

obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>% 
  left_join(meta, by = "orig.ident") %>%
  mutate(orig.ident = str_split_i(orig.ident, "_", i = 1) %>%
           factor(orig.ident),
         Tissue = factor(Tissue, 
                         levels = c("b", "s", "m"),
                         labels = c("Motor cortex", 
                                    "Cervical spinal cord",
                                    "Skeletal muscle"))) %>%
  column_to_rownames(var = "cell")

sample_pal <- DiscretePalette_scCustomize(30, palette = "polychrome")

obj <- Store_Palette_Seurat(seurat_object = obj, 
                            palette = sample_pal, palette_name = "sample_pal",
                            overwrite = T)

tissue_pal <- JCO_Four()[1:3]

obj <- Store_Palette_Seurat(obj,
                            palette = tissue_pal,
                            palette_name = "tissue_pal",
                            overwrite = T)

obj <- obj %>%
  Add_Cell_QC_Metrics(species = "human")

message("Removing odd cell with 1 count")

obj <- subset(obj,
              subset = nCount_RNA > 1)

# For lab meeting 10/30/25 ------------------------------------------------

# df <- obj@meta.data
# 
# df <- df %>%
#   mutate(Group = factor(Group,
#                         levels = c("Control", "sALS", "C9orf72"),
#                         labels = c("Control", "sALS", "C9orf72-ALS")))
# 
# df <- df %>%
#   mutate(id = case_when(id == "GWF-19-47" ~ "GWF19-47",
#                         id == "GWF-20-54" ~ "GWF20-54",
#                         id == "GWF-21-56" ~ "GWF21-56",
#                         .default = id))
# 
# df %>%
#   ggplot(aes(x = id,
#              y = nCount_RNA)) + 
#   geom_violin(aes(fill = id),
#               show.legend = F) + 
#   facet_wrap(. ~ tissue,
#              ncol = 1) + 
#   ggtitle("# UMIs per cell") +
#   scale_y_log10() +
#   theme_bw(base_size = 16) + 
#   theme(axis.text = element_text(color = "black"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         plot.title = element_text(hjust = 0.5, size = 20))
# ggsave(filename = paste0(plots_dir, "20251029_ncount.png"),
#        units = "in", dpi = 600,
#        height = 9, width = 13)
# 
# df %>%
#   ggplot(aes(x = id,
#              y = nFeature_RNA)) + 
#   geom_violin(aes(fill = id),
#               show.legend = F) + 
#   facet_wrap(. ~ tissue,
#              ncol = 1) + 
#   ggtitle("# genes per cell") +
#   scale_y_log10() +
#   theme_bw(base_size = 16) + 
#   theme(axis.text = element_text(color = "black"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         plot.title = element_text(hjust = 0.5, size = 20))
# ggsave(filename = paste0(plots_dir, "20251029_nfeature.png"),
#        units = "in", dpi = 600,
#        height = 9, width = 13)
# 
# df %>%
#   ggplot(aes(x = id,
#              y = percent_mito)) + 
#   geom_violin(aes(fill = id),
#               show.legend = F) + 
#   facet_wrap(. ~ tissue,
#              ncol = 1) + 
#   ggtitle("Fraction mito genes per cell") +
#   # scale_y_log10() +
#   theme_bw(base_size = 16) + 
#   theme(axis.text = element_text(color = "black"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         plot.title = element_text(hjust = 0.5, size = 20))
# ggsave(filename = paste0(plots_dir, "20251029_mito.png"),
#        units = "in", dpi = 600,
#        height = 9, width = 13)
# 
# df %>%
#   ggplot(aes(x = id,
#              y = log10GenesPerUMI)) + 
#   geom_violin(aes(fill = id),
#               show.legend = F) + 
#   facet_wrap(. ~ tissue,
#              ncol = 1) + 
#   ggtitle("Cell complexity (log10(# genes) / log10(# UMIs)") +
#   # scale_y_log10() +
#   theme_bw(base_size = 16) + 
#   theme(axis.text = element_text(color = "black"),
#         axis.title = element_blank(),
#         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
#         plot.title = element_text(hjust = 0.5, size = 20))
# ggsave(filename = paste0(plots_dir, "20251029_complex.png"),
#        units = "in", dpi = 600,
#        height = 9, width = 13)

# ----

message2("Subsetting brain")

list <- list()

list[[1]] <- obj %>%
  subset(subset = tissue == "Motor cortex")

message2("Subsetting muscle")

list[[3]] <- obj %>%
  subset(subset = tissue == "Skeletal muscle")

message2("Subsetting spinal cord")

list[[2]] <- obj %>%
  subset(subset = tissue == "Cervical spinal cord")

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

for (i in seq_along(list)){
  message2(paste0("QC plots and stats for ", tissues$title[i]))
  
  p <- QC_Plots_Genes(list[[i]],
                      group.by = "id",
                      colors_use = list[[i]]@misc$sample_pal,
                      plot_boxplot = T,
                      y_axis_log = T) +
    ylab("# of unique genes") +
    ggtitle(tissues$title[i]) +
    theme(plot.title = element_text(color = list[[i]]@misc$tissue_pal[i],
                                    face = "bold"))
  ggsave(p,
         filename = paste0(plots_dir, "nfeature_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)
  
  p <- QC_Plots_UMIs(list[[i]],
                     group.by = "id",
                     colors_use = list[[i]]@misc$sample_pal,
                     plot_boxplot = T,
                     y_axis_log = T) +
    ylab("# of unique UMIs") +
    ggtitle(tissues$title[i]) +
    theme(plot.title = element_text(color = list[[i]]@misc$tissue_pal[i],
                                    face = "bold"))
  
  ggsave(p,
         filename = paste0(plots_dir, "numi_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- QC_Plots_Mito(list[[i]],
                     colors_use = list[[i]]@misc$sample_pal,
                     plot_boxplot = T,
                     group.by = "id") +
    ylab("% mitochondrial\ngene counts") +
    ggtitle(tissues$title[i]) +
    theme(plot.title = element_text(color = list[[i]]@misc$tissue_pal[i],
                                    face = "bold"))
  
  ggsave(p,
         filename = paste0(plots_dir, "mito_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- QC_Plots_Complexity(list[[i]],
                           colors_use = list[[i]]@misc$sample_pal,
                           plot_boxplot = T,
                           group.by = "id") +
    ggtitle(tissues$title[i]) +
    theme(plot.title = element_text(color = list[[i]]@misc$tissue_pal[i],
                                    face = "bold"))
  
  ggsave(p,
         filename = paste0(plots_dir, "complexity_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)
  
  stats <- Median_Stats(list[[i]],
                        group_by = "id")
  
  counts <- list[[i]]@meta.data %>%
    group_by(id) %>%
    dplyr::summarize(Cell_count = n()) %>%
    adorn_totals(name = "Totals (All Cells)")
  
  stats <- stats %>%
    left_join(counts, by = "id")
  
  write.csv(stats, 
            file = paste0(csv_dir, "median_stats_", tissues$file[i], ".csv"),
            row.names = F)
    
}

message2("Saving full QC'd object")

saveRDS(obj,
        file = paste0(data_out_dir,
                      "obj.rds"))


