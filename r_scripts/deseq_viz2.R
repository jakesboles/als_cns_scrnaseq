# Upset plot of microglia DEGs (sALS vs. Control and C9orf72 vs. Control,
# in both motor cortex and cervical spinal cord -- 4 sets, all pairwise/
# higher-order overlaps shown) from deseq2.R's output. Standalone script,
# not part of the array job chain (matching demographics_figure.R's
# precedent) -- this only reads a handful of small CSVs and draws one
# plot, no Seurat/BPCells/heavy compute involved.
#
# Reworked from the user's pushed sample script (from a different, older
# project's own DESeq2 visualization). Changes made to fit this project:
# - Input paths point at deseq2.R's actual output structure:
#   results/deseq2/<tissue>/<cell type>/<contrast>.csv, with tissue file
#   names "brain"/"sc" (not "Brain"/"SpinalCord") and contrast file names
#   "sALS_vs_Control.csv"/"C9orf72_vs_Control.csv" (not a
#   "DESeq2_DEGs.csv" per cell type).
# - Reads the raw (non-LFC-shrunk) results CSVs, not
#   *_lfc_shrunk.csv -- matching deseq_viz1.R's own established DEG-
#   counting convention already used elsewhere in this project (same
#   p < 0.05 & abs(log2FoldChange) > log2(1.5) thresholds on the same raw
#   files), for consistency with that existing analysis rather than
#   introducing a second, differently-thresholded definition of "DEG" in
#   parallel.
# - Gene identifiers are pulled from the "X" column, not "gene" --
#   deseq2.R's write.csv() calls don't set row.names = F, so the gene
#   symbols (results()'s rownames) get written as an unnamed first column
#   that read.csv() reads back in as "X" (same column deseq_viz1.R itself
#   reads from).
# - Added a file-existence check before reading, matching deseq_viz1.R's
#   own precedent -- deseq2.R skips a (tissue, cell type) combination
#   entirely if it's too sparse for pseudobulk DESeq2, so a missing file
#   here is a real possibility, not just a typo'd path.
# - Output goes to figures/ (already covered by .gitignore), not
#   results/deseq_viz2/ -- distinct from this project's usual per-script
#   results/ convention, since these are meant as curated, presentation-
#   ready figures rather than per-script diagnostic output.

library(tidyverse)
library(ComplexUpset)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

results_dir <- "results/deseq2/"
plots_dir <- "figures/"
dir.create(plots_dir, showWarnings = F, recursive = T)

# Microglia ---------------------------------------------------------------

celltype <- "Microglia"

file_paths <- c(
  sALS_brain = paste0(results_dir, "brain/", celltype, "/sALS_vs_Control.csv"),
  `C9-ALS_brain` = paste0(results_dir, "brain/", celltype, "/C9orf72_vs_Control.csv"),
  sALS_sc = paste0(results_dir, "sc/", celltype, "/sALS_vs_Control.csv"),
  `C9-ALS_sc` = paste0(results_dir, "sc/", celltype, "/C9orf72_vs_Control.csv")
)

missing <- file_paths[!file.exists(file_paths)]
if (length(missing) > 0){
  stop(paste0("Missing DESeq2 results file(s) for ", celltype, ": ",
              paste(missing, collapse = ", "),
              " -- check whether deseq2.R's abundance filter skipped ",
              celltype, " in this tissue."))
}

list_names <- names(file_paths)

list <- lapply(file_paths, read.csv)
names(list) <- list_names

for (i in seq_along(list)){
  list[[i]] <- list[[i]] %>%
    filter(padj < 0.05 &
             abs(log2FoldChange) > log2(1.5)) %>%
    pull(X)
}

genes <- unique(list_c(list))

intersect_df <- data.frame(gene = genes)
for (i in seq_along(list)){

  col <- list_names[i]

  intersect_df <- intersect_df %>%
    mutate(x = if_else(gene %in% list[[i]], T, F)) %>%
    dplyr::rename(!!sym(col) := "x")
}

pdf(file = paste0(plots_dir, "microglia_degs_upset.pdf"),
    height = 6, width = 7)

upset(intersect_df, list_names,
      set_sizes = F,
      stripes = "white",
      wrap = T,
      encode_sets = F,
      height_ratio = 0.9,
      labeller = ggplot2::as_labeller(c(
        "C9-ALS_brain" = "C9orf72-ALS\nMotor cortex",
        "sALS_brain" = "sALS\nMotor cortex",
        "C9-ALS_sc" = "C9orf72-ALS\nCervical spinal cord",
        "sALS_sc" = "sALS\nCervical spinal cord"
      )),
      matrix = (intersection_matrix(
        geom = geom_point(shape = 18, size = 10),
        segment = geom_segment(linewidth = 1.5),
        outline_color = list(active = "white", inactive = "white")
      )),
      base_annotations = list(
        'Intersection size' = intersection_size(text = list(size = 0)) +
          scale_y_continuous(expand = c(0, 0)) +
          ylab("# DEGs")
      ),
      queries = list(
        upset_query(intersect = c("C9-ALS_brain", "C9-ALS_sc"),
                    color = "magenta3",
                    fill = "magenta3"),
        upset_query(intersect = c("C9-ALS_brain", "sALS_brain"),
                    color = "#0073C2",
                    fill = "#0073C2"),
        upset_query(intersect = c("C9-ALS_sc", "sALS_sc"),
                    color = "#EFC000",
                    fill = "#EFC000")
      ),
      theme = upset_modify_themes(
        list(
          'Intersection size' = theme(axis.text = element_text(color = "black", size = 16),
                                      axis.title = element_text(size = 20),
                                      axis.ticks.y = element_line(),
                                      panel.border = element_rect(color = "black", fill = "transparent")),
          'intersections_matrix' = theme(axis.text = element_text(color = "black", size = 16),
                                         axis.title = element_blank())
        )
      )
) +
  ggtitle("Microglia DEGs") +
  theme(plot.title = element_text(hjust = 0.5, size = 20))

dev.off()
