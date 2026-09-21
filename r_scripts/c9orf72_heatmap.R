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
# - Cell type inclusion: a cell_type3 is only shown for a given tissue if
#   it has > 500 cells total in that tissue (summed across all samples/
#   groups), matching the reference figure's "Cell count > 500" note.
#   Read as a whole-tissue total, not a per-group minimum -- the
#   reference figure doesn't specify per-group, and a whole-tissue total
#   is the simpler, more common convention for this kind of filter. Flag
#   if a per-group threshold is actually what's wanted.
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

gene <- "C9orf72"
min_cells <- 500 # change as needed -- see header note above

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

  # Cell types with > min_cells total in this tissue -- see header note.
  keep_types <- meta %>%
    dplyr::count(cell_type3) %>%
    dplyr::filter(n > min_cells) %>%
    dplyr::pull(cell_type3)

  if (length(keep_types) == 0){
    stop(paste0("No cell types in ", tissue_title, " have > ", min_cells,
                " cells -- check min_cells or this tissue's metadata."))
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

  return(ht)
}

heatmaps <- Map(make_tissue_heatmap, tissues$file, tissues$title)

p <- wrap_plots(lapply(heatmaps, function(h) wrap_elements(full = h$gtable)),
               ncol = length(heatmaps))

ggsave(p,
       filename = paste0(figures_dir, "c9orf72_heatmap.png"),
       units = "in", dpi = 600,
       height = 6, width = 12)
