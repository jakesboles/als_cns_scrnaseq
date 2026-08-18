# Load libraries
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(scCustomize)
  library(stringr)
  library(janitor)
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

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "plots/03_qc2/"
dir.create(plots_dir, showWarnings = F,
           recursive = T)

csv_dir <- "tab_data/03_qc2/"
dir.create(csv_dir, showWarnings = F,
           recursive = T)

data_out_dir <- "data/03_qc2/"
dir.create(data_out_dir, showWarnings = F,
           recursive = T)

# 02_qc1.R didn't generate a new counts matrix -- it only added QC metadata
# (percent_mito, log10GenesPerUMI, etc.) and saved that as its own RDS, so we
# reuse the same on-disk BPCells matrix from 01_obj_creation.R here, paired
# with 02's updated metadata.
message2("Reading in object")
counts <- open_matrix_dir("data/01_obj_creation/bpcells")
meta <- readRDS("data/02_qc1/metadata.rds")
obj <- CreateSeuratObject(counts = counts, meta.data = meta)

# Add sample/tissue identifiers --------------------------------------------

message2("Deriving sample and tissue identifiers")

# orig.ident is set at object creation (01_obj_creation.R) from each
# sample's directory name as "<subject id>_<tissue code>" (e.g. "AU-066_b"),
# with tissue code one of b/s/m for motor cortex, spinal cord, or muscle.
# `id` is redundant with the `id` column 02_qc1.R already saves (recomputed
# here the same way), kept for clarity/robustness. `tissue` is trusted as-is
# from 02_qc1.R's metadata rather than re-derived.
# NOTE: I wasn't able to inspect the live object or metadata.csv directly,
# so please confirm orig.ident is still "<id>_<tissue code>" with exactly
# one underscore once you have the object loaded.
obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>%
  mutate(id = str_split_i(orig.ident, "_", i = 1)) %>%
  column_to_rownames(var = "cell")

obj@meta.data <- obj@meta.data %>%
  mutate(log_nFeature = log10(nFeature_RNA),
         log_nCount = log10(nCount_RNA))

# Set cutoffs ---------------------------------------
message2("Setting thresholds and drawing new plots")

samples <- sort(unique(as.character(obj$orig.ident)))

thresh_df <- data.frame(orig.ident = samples,
                        umi_med = rep(NA, length(samples)),
                        umi_mad = rep(NA, length(samples)),
                        feature_med = rep(NA, length(samples)),
                        feature_mad = rep(NA, length(samples)),
                        mito_med = rep(NA, length(samples)),
                        mito_mad = rep(NA, length(samples)))

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

tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle")
)

for (i in seq_along(tissues$title)){

  message2(paste0("Making plots for ", tissues$title[i]))

  p <- meta %>%
    filter(tissue == tissues$title[i]) %>%
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
    ggtitle(tissues$title[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    hjust = 0.5))
  ggsave(p,
         filename = paste0(plots_dir,
                           "numi_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- meta %>%
    filter(tissue == tissues$title[i]) %>%
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
    ggtitle(tissues$title[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    hjust = 0.5))
  ggsave(p,
         filename = paste0(plots_dir,
                           "mito_", tissues$file[i], ".png"),
         units = "in", dpi = 600,
         height = 6, width = 12)

  p <- meta %>%
    filter(tissue == tissues$title[i]) %>%
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
    ggtitle(tissues$title[i]) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme_bw() +
    theme(axis.text = element_text(color = "black"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(face = "bold",
                                    hjust = 0.5))
  ggsave(p,
         filename = paste0(plots_dir,
                           "nfeature_", tissues$file[i], ".png"),
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
  group_by(id, tissue, discard) %>%
  dplyr::summarize(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = "discard",
              values_from = "n") %>%
  mutate(`TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`)) %>%
  mutate(retained_percent = (`FALSE` / (`FALSE` + `TRUE`)) * 100) %>%
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

for (i in seq_along(tissues$title)){
  message(paste0("Getting stats for ", tissues$title[i]))

  s <- subset(obj_filtered, subset = tissue == tissues$title[i])

  med_stats <- Median_Stats(s,
                            group.by = "id")

  cell_counts <- s@meta.data %>%
    group_by(id) %>%
    dplyr::summarize(Cell_count = n()) %>%
    adorn_totals(name = "Totals (All Cells)")

  med_stats <- med_stats %>%
    left_join(cell_counts, by = "id")

  write.csv(med_stats,
            file = paste0(csv_dir, "median_stats_", tissues$file[i], ".csv"),
            row.names = F)
}

message2("Saving filtered object")

saveRDS(obj_filtered,
        file = paste0(data_out_dir, "filtered_obj.rds"))
