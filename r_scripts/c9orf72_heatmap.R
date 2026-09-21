# Multi-panel heatmap of C9orf72 expression by cell type and disease
# group, one panel per tissue -- modeled on Fig. 1g of
# https://www.nature.com/articles/s41593-026-02300-5 (their PBMC C9orf72
# expression heatmap), adapted for this project's 3 tissues instead of
# one PBMC panel.
#
# Design notes:
# - Source: data/17_obj_reassembly/<tissue>/{bpcells_data,metadata.rds} --
#   the terminal, fully-annotated per-tissue object (cell_type3, "Remove"
#   cells already dropped), the standard source for an "every cell type,
#   one tissue" figure in this project.
# - Cell type inclusion (changed from the reference figure's flat "> 500
#   cells" note, per the user): a cell_type3 is shown for a given tissue
#   only if its raw C9orf72 count sums to >= min_sample_count (10) in at
#   least min_sample_fraction (1/3) of that tissue's samples. This reuses
#   the same per-(sample, cell_type3) total_count the raw count
#   diagnostic table below already computes -- one raw-count pass now
#   drives both the heatmap's cell type filter and the diagnostic table,
#   rather than the old flat total-cell-count threshold and the raw
#   counts being two unrelated things. "At least 1/3 of samples" is
#   implemented as n_samples_passing >= n_samples * min_sample_fraction
#   (no explicit rounding needed -- an integer count can't satisfy a
#   fractional threshold like 3 >= 3.33, so this already behaves like a
#   ceiling for non-exact-thirds sample counts).
# - Expression value: mean log-normalized ("data" layer) C9orf72
#   expression per (cell_type3, group), computed manually via rowMeans()
#   rather than Seurat's AverageExpression() -- AverageExpression()
#   exponentiates back to linear scale before averaging by default
#   (return.seurat = FALSE), which isn't wanted here; computing directly
#   on the log-normalized layer avoids that surprise and matches how this
#   project handles expression elsewhere (raw log-normalized "data", not
#   Seurat's own un-logging shortcut).
# - Heatmap color scale: row-scaled (z-scored across groups within each
#   cell type), via pheatmap's own scale = "row" -- matches the
#   reference figure's symmetric ~-1 to 1 diverging scale, and shows
#   "which group is relatively high/low for this cell type" rather than
#   raw expression (which would mostly just reflect each cell type's
#   overall expression level, not group differences).
# - Row clustering only (cluster_rows = TRUE, cluster_cols = FALSE) --
#   matches the reference figure's dendrogram on cell types but a fixed,
#   meaningful column order (Control -> sALS -> C9orf72-ALS), not a
#   clustered column order.
# - Uses pheatmap for the per-tissue dendrogram + heatmap -- this
#   project's first use of it (may need installing on the cluster:
#   install.packages("pheatmap")). ggplot2's geom_tile() can't natively
#   draw a dendrogram, and pheatmap is the standard tool for exactly this
#   kind of clustered expression heatmap; flagging the new dependency
#   since every other script in this project sticks to packages already
#   in use. Panels are combined into one multi-panel figure via
#   patchwork::wrap_elements() around each pheatmap's $gtable, per the
#   request for one combined file rather than 3 separate ones.
# - Each tissue panel is scaled/colored independently (its own z-score,
#   its own legend) -- cell type vocabularies differ substantially by
#   tissue, so a single shared color scale across all three wouldn't be
#   meaningful the way it is within one tissue's own set of cell types.
# - No R environment available in this session to run/verify any of
#   this -- first draft, pending your own run.
# - Raw count diagnostic table: per (sample, cell_type3, tissue), the raw
#   (unnormalized) C9orf72 count -- n_cells, total_count, mean_count, and
#   pct_detected (fraction of cells with count > 0), plus that sample's
#   group for convenience. Pulled from data/06_obj_reassembly/bpcells (the
#   whole-cohort real raw counts), not data/17_obj_reassembly's own
#   bpcells_data -- that holds normalized data, not counts (see this
#   project's own CLAUDE.md gotcha on this). "Sample" = orig.ident
#   (donor + tissue), this project's standard per-sample key. Unlike the
#   heatmap, this table is NOT restricted to keep_types (the >= 10
#   counts/1/3-of-samples filter above) -- the whole point of a
#   diagnostic table is to also be able to check the cell types/samples
#   the heatmap filters out, e.g. to confirm a low heatmap value isn't
#   actually an artifact of very few cells or very low detection.

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
