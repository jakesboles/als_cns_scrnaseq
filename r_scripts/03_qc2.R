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
  library(paletteer)
  library(stringr)
  library(scales)
  library(janitor)
  library(scater)
  library(ggbeeswarm)
})
# Function to print clear log progress updates
message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Create directories ------------------------------------------------------
message2("Creating directories")

proj_dir <- "/projects/b1169/boles/als_multitissue_scfrp/"

data_in_dir <- paste0(proj_dir, "data/02_qc1/")

plots_dir <- paste0(proj_dir, "plots/03_qc2/")
dir.create(plots_dir, showWarnings = F,
           recursive = T)

csv_dir <- paste0(proj_dir, "tab_data/03_qc2/")
dir.create(csv_dir, showWarnings = F,
           recursive = T)

data_out_dir <- paste0(proj_dir, "data/03_qc2/")
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

obj <- readRDS(paste0(data_in_dir, "obj.rds"))

obj@meta.data <- obj@meta.data %>%
  mutate(log_nFeature = log10(nFeature_RNA),
         log_nCount = log10(nCount_RNA))

# Set cutoffs ---------------------------------------
message2("Setting thresholds and drawing new plots")

samples <- levels(obj$orig.ident)

thresh_df <- data.frame(orig.ident = samples,
                        umi_med = c(rep(NA, 90)),
                        umi_mad = c(rep(NA, 90)),
                        feature_med = c(rep(NA, 90)),
                        feature_mad = c(rep(NA, 90)),
                        mito_med = c(rep(NA, 90)),
                        mito_mad = c(rep(NA, 90)))

meta <- obj@meta.data

for (i in seq_along(samples)){
  message(paste0("Getting cutoffs for ", samples[i]))
  
  df <- meta %>%
    filter(orig.ident == samples[i])
  
  thresh_df$umi_med[i] <- median(df$log_nCount)
  
  thresh_df$umi_mad[i] <- stats::mad(df$log_nCount)
  
  thresh_df$feature_med[i] <- median(df$log_nFeature)
  
  thresh_df$feature_mad[i] <- stats::mad(df$log_nFeature)
  
  thresh_df$mito_med[i] <- median(df$percent_mito)
  
  thresh_df$mito_mad[i] <- stats::mad(df$percent_mito)
}

# Pretty "strict" cutoffs of 2 x MAD still discards very few cells based 
# on number of counts and unique genes
thresh_df <- thresh_df %>%
  mutate(umi_lower = umi_med - 2*umi_mad,
         feature_lower = feature_med - 2*feature_mad,
         mito_upper = 5,
         umi_upper = umi_med + 2*umi_mad,
         feature_upper = feature_med + 2*feature_mad)

meta <- meta %>%
  rownames_to_column(var = "cell") %>%
  left_join(thresh_df,
            by = "orig.ident") %>% 
  mutate(mito_discard = if_else(percent_mito > mito_upper, T, F),
         umi_discard = if_else(log_nCount < umi_lower | 
                                 log_nCount > umi_upper, T, F),
         gene_discard = if_else(log_nFeature < feature_lower |
                                  log_nFeature > feature_upper, T, F))

# Make some plots with cells colored by QC status -------------------------
message2("Making QC plots") 

tissues <- levels(meta$tissue)
files <- c("brain", "sc", "muscle")

for (i in seq_along(tissues)){
  
  message2(paste0("Making plots for ", tissues[i]))
  
  p <- meta %>%
    filter(tissue == tissues[i]) %>%
    ggplot(aes(x = id,
               y = log_nCount)) +
    geom_quasirandom(size = 0.1,
                     aes(color = umi_discard)) +
    scale_color_manual(values = c("gray60", "red")) +
    stat_summary(aes(group = id),
                 fun = median,
                 color = "black",
                 geom = "crossbar",
                 width = 0.5) +
    labs(y = "log10(# unique UMIs)",
         color = "Below threshold?") +
    ggtitle(tissues[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          # legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    color = obj@misc$tissue_pal[i],
                                    hjust = 0.5))
  ggsave(filename = paste0(plots_dir, 
                           "numi_", files[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- meta %>%
    filter(tissue == tissues[i]) %>%
    ggplot(aes(x = id,
               y = percent_mito)) +
    geom_quasirandom(size = 0.1,
                     aes(color = mito_discard)) +
    scale_color_manual(values = c("gray60", "red")) +
    stat_summary(aes(group = id),
                 fun = median,
                 color = "black",
                 geom = "crossbar",
                 width = 0.5) +
    labs(y = "% mitochondrial\ngene counts",
         color = "Above threshold?") +
    ggtitle(tissues[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          # legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    color = obj@misc$tissue_pal[i],
                                    hjust = 0.5))
  ggsave(filename = paste0(plots_dir, 
                           "mito_", files[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- meta %>%
    filter(tissue == tissues[i]) %>%
    ggplot(aes(x = id,
               y = log_nFeature)) +
    geom_quasirandom(size = 0.1,
                     aes(color = gene_discard)) +
    scale_color_manual(values = c("gray60", "red")) +
    stat_summary(aes(group = id),
                 fun = median,
                 color = "black",
                 geom = "crossbar",
                 width = 0.5) +
    labs(y = "log10(# unique genes)",
         color = "Below threshold?") +
    ggtitle(tissues[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) + 
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          # legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    color = obj@misc$tissue_pal[i],
                                    hjust = 0.5))
  ggsave(filename = paste0(plots_dir, 
                           "nfeature_", files[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)
}


# Mark cells for removal and save some stats ------------------------------
message2("Marking cells and saving statistics")  

meta <- meta %>%
  mutate(discard = if_else(mito_discard == T |
                             umi_discard == T |
                             gene_discard == T,
                           T, F))
  

stats <- meta %>%
  group_by(orig.ident, discard) %>% 
  dplyr::summarize(n = n()) %>%
  pivot_wider(names_from = "discard",
              values_from = "n") %>%
  mutate(`TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(retained_percent = (`FALSE` / (`FALSE` + `TRUE`)) * 100) %>%
  separate(orig.ident, 
           into = c("id", "tissue"),
           sep = "_") %>%
  mutate(tissue = case_when(tissue == "b" ~ "Motor cortex",
                            tissue == "s" ~ "Cervical spinal cord",
                            tissue == "m" ~ "Skeletal muscle")) %>%
  dplyr::rename("retained_count" = "FALSE",
                "discarded_count" = "TRUE")

write.csv(thresh_df,
          file = paste0(csv_dir, "cutoffs.csv"),
          row.names = F)

write.csv(stats,
          file = paste0(csv_dir, "filtered_counts.csv"),
          row.names = F)

barcodes <- meta %>%
  dplyr::select(c(cell, nCount_RNA, nFeature_RNA, percent_mito,
                  mito_discard, umi_discard, gene_discard, discard))

write.csv(barcodes,
          file = paste0(csv_dir, "barcode_qc.csv"),
          row.names = F)

# Discard low quality cells and save object -------------------------------
message2("Filtering object of low quality cells")

obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>%
  left_join(meta) %>%
  column_to_rownames(var = "cell")

obj_filtered <- subset(obj, discard == F)

message2("Saving new median stats table")

for (i in seq_along(tissues)){
  message(paste0("Getting stats for ", tissues[i]))
  
  s <- subset(obj_filtered, subset = tissue == tissues[i])
  
  med_stats <- Median_Stats(s,
                            group_by_var = "id")
  
  counts <- s@meta.data %>%
    group_by(id) %>%
    dplyr::summarize(Cell_count = n()) %>%
    adorn_totals(name = "Totals (All Cells)")
  
  med_stats <- med_stats %>%
    left_join(counts, by = "id")
  
  write.csv(med_stats, 
            file = paste0(csv_dir, "median_stats_", files[i], ".csv"),
            row.names = F)
}

message2("Saving filtered object")

saveRDS(obj_filtered,
        file = paste0(data_out_dir, "filtered_obj.rds"))
