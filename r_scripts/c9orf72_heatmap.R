# Multi-panel heatmap of C9orf72 expression by cell type and disease
# group, one panel per tissue, with cell types included only where raw
# C9orf72 counts reach a minimum threshold in enough samples -- modeled
# on a published PBMC C9orf72 expression heatmap, adapted to this
# project's 3 tissues. Also writes a raw per-(sample, cell type, tissue)
# C9orf72 count diagnostic table, independent of the heatmap's cell type
# filter. Interactive script, not a SLURM job.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(BPCells)
  library(pheatmap)
  library(patchwork)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

figures_dir <- "figures/"
dir.create(figures_dir, showWarnings = F, recursive = T)

results_dir <- "results/c9orf72_heatmap/"
dir.create(results_dir, showWarnings = F, recursive = T)

gene <- "C9orf72"
min_sample_count <- 10 # change as needed -- see header note above
min_sample_fraction <- 1 / 3 # change as needed -- see header note above

tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")
)

# Diverging color ramp for the row-scaled heatmap -- roughly matches the
# reference figure's blue-white-red scale.
heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

make_tissue_heatmap <- function(tissue_file, tissue_title){

  message(paste0("Processing ", tissue_title))

  data_dir <- paste0("data/17_obj_reassembly/", tissue_file, "/")

  mat <- open_matrix_dir(paste0(data_dir, "bpcells_data"))
  meta <- readRDS(paste0(data_dir, "metadata.rds"))

  if (!(gene %in% rownames(mat))){
    stop(paste0(gene, " not found in ", tissue_title, "'s expression ",
                "matrix -- check whether it's included in this project's ",
                "probe panel."))
  }

  # Raw counts -- computed first now, since the cell type filter below
  # and the diagnostic table both derive from the same per-(sample,
  # cell_type3) total_count. Real raw counts, not this tissue's own
  # (normalized) bpcells_data -- see header note above.
  message(paste0("Reading raw counts for ", tissue_title))

  raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
  raw_mat <- raw_mat[, rownames(meta)]

  if (!(gene %in% rownames(raw_mat))){
    stop(paste0(gene, " not found in the whole-cohort raw count matrix -- ",
                "check data/06_obj_reassembly/bpcells."))
  }

  raw_gene_row <- raw_mat[gene, , drop = F]
  raw_counts <- as.numeric(as.matrix(raw_gene_row))
  names(raw_counts) <- colnames(raw_mat)

  raw_df <- data.frame(cell = names(raw_counts), count = raw_counts) %>%
    left_join(meta %>% rownames_to_column("cell") %>%
                dplyr::select(cell, cell_type3, orig.ident, group),
              by = "cell") %>%
    mutate(group = factor(group, levels = c("Control", "sALS", "C9orf72"),
                          labels = c("Control", "sALS", "C9orf72-ALS")))

  raw_summary <- raw_df %>%
    group_by(orig.ident, group, cell_type3) %>%
    summarise(n_cells = n(),
             total_count = sum(count),
             mean_count = mean(count),
             pct_detected = mean(count > 0),
             .groups = "drop") %>%
    mutate(tissue = tissue_title) %>%
    dplyr::rename(sample = orig.ident) %>%
    relocate(tissue, sample, group, cell_type3)

  # Cell types with >= min_sample_count total counts in at least
  # min_sample_fraction of this tissue's samples -- see header note.
  n_samples <- n_distinct(meta$orig.ident)

  keep_types <- raw_summary %>%
    group_by(cell_type3) %>%
    summarise(n_passing = sum(total_count >= min_sample_count), .groups = "drop") %>%
    dplyr::filter(n_passing >= n_samples * min_sample_fraction) %>%
    dplyr::pull(cell_type3)

  if (length(keep_types) == 0){
    stop(paste0("No cell types in ", tissue_title, " have >= ",
                min_sample_count, " counts in at least ",
                round(min_sample_fraction * 100, 1), "% of ", n_samples,
                " samples -- check min_sample_count/min_sample_fraction ",
                "or this tissue's metadata."))
  }

  gene_row <- mat[gene, , drop = F]
  gene_expr <- as.numeric(as.matrix(gene_row))
  names(gene_expr) <- colnames(mat)

  expr_df <- data.frame(cell = names(gene_expr), expr = gene_expr) %>%
    left_join(meta %>% rownames_to_column("cell") %>%
                dplyr::select(cell, cell_type3, group),
              by = "cell") %>%
    dplyr::filter(cell_type3 %in% keep_types) %>%
    mutate(group = factor(group, levels = c("Control", "sALS", "C9orf72"),
                          labels = c("Control", "sALS", "C9orf72-ALS")))

  avg_df <- expr_df %>%
    group_by(cell_type3, group) %>%
    summarise(mean_expr = mean(expr), .groups = "drop")

  mat_wide <- avg_df %>%
    pivot_wider(names_from = group, values_from = mean_expr) %>%
    column_to_rownames("cell_type3") %>%
    as.matrix()

  # Column order fixed regardless of pivot_wider's own ordering.
  mat_wide <- mat_wide[, c("Control", "sALS", "C9orf72-ALS"), drop = F]

  ht <- pheatmap(mat_wide,
                scale = "row",
                cluster_rows = TRUE,
                cluster_cols = FALSE,
                color = heat_colors,
                main = tissue_title,
                silent = TRUE)

  return(list(heatmap = ht, raw_counts = raw_summary))
}

results <- Map(make_tissue_heatmap, tissues$file, tissues$title)

p <- wrap_plots(lapply(results, function(r) wrap_elements(full = r$heatmap$gtable)),
               ncol = length(results))

ggsave(p,
       filename = paste0(figures_dir, "c9orf72_heatmap.png"),
       units = "in", dpi = 600,
       height = 6, width = 12)

message("Saving raw count diagnostic table")

raw_counts_all <- lapply(results, function(r) r$raw_counts) %>%
  list_rbind()

write.csv(raw_counts_all,
          file = paste0(results_dir, "raw_counts_by_sample.csv"),
          row.names = F)
