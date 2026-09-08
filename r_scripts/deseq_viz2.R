# Three DESeq2/GSEA visualization sections for microglia, all reusing the
# same celltype/results_dir/plots_dir set up front. Standalone script,
# not part of the array job chain (matching demographics_figure.R's
# precedent) -- everything here reads a handful of small CSVs and draws
# plots, no Seurat/BPCells/heavy compute involved.
#
# 1. Upset plot of microglia DEGs (sALS vs. Control and C9orf72 vs.
#    Control, in both motor cortex and cervical spinal cord -- 4 sets,
#    all pairwise/higher-order overlaps shown) from deseq2.R's output.
# 2. plot_fc_scatter(): a function (not a fixed block, since the user
#    wants this flexible) plotting any two (tissue, contrast)
#    comparisons' log2FoldChange against each other for genes
#    significant in at least one of the two, colored by which
#    comparison(s) each gene is significant in.
# 3. A lollipop chart of the top dysregulated pathways (by p.adjust)
#    across all 4 comparisons (2 tissues x 2 groups) from
#    deseq2_gsea.R's GSEA output, colored by group and shaped by tissue.
#
# Reworked from the user's pushed sample script (from a different, older
# project's own DESeq2 visualization) for section 1, and from the user's
# own description (sections 2-3 had no sample code). Changes/design
# choices, section 1:
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
#
# Design choices, section 2 (plot_fc_scatter()):
# - Same raw results CSVs and "X" gene column as section 1, for the same
#   consistency reason. Only genes significant (same thresholds as
#   section 1) in at least one of the two chosen comparisons are plotted
#   -- "plots of significant DEGs", not the whole transcriptome.
# - celltype is not a function argument -- reused from the top of the
#   script, since this whole file is scoped to microglia. tissue/contrast
#   are the two axes the user can vary, matching their own example
#   (same group, different tissue).
# - Returns the ggplot object rather than saving it, since the specific
#   pair of comparisons varies by call -- there's no single sensible
#   fixed output filename the way section 1's upset plot has one.
#
# Design choices, section 3 (GSEA lollipop):
# - Reads deseq2_gsea.R's saved
#   results/deseq2/<tissue>/<celltype>/<contrast>_GSEA.csv files (GSEA()
#   results, not DESeq2 results) for microglia, both tissues, both
#   groups. GSEA()'s own default pvalueCutoff already restricts what
#   deseq2_gsea.R saved to nominally significant pathways, so no
#   additional significance filter is applied here.
# - "Top dysregulated" = top_n pathways per comparison by p.adjust
#   (default 10, change as needed); the *union* of those across all 4
#   comparisons is plotted, not just one comparison's top list, so a
#   pathway's behavior can be compared across tissue/group even where it
#   wasn't top-ranked in every comparison.
# - This section builds and saves its own plot (unlike section 2) --
#   unlike the flexible fold-change scatter, this is a complete, specific
#   request (all 4 comparisons, a fixed color/shape encoding), so it gets
#   a real output file the same way section 1 does.

library(tidyverse)
library(ComplexUpset)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

results_dir <- "results/deseq2/"
plots_dir <- "figures/"
dir.create(plots_dir, showWarnings = F, recursive = T)

# celltype is set once here and reused by sections 2 and 3 below -- all
# three sections are scoped to microglia specifically, matching this
# script's established focus (section 1's upset plot).

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

intersect_df <- intersect_df %>%
  mutate(
    highlight_group = case_when(
      `C9-ALS_brain` & `C9-ALS_sc` & !sALS_brain & !sALS_sc ~ "C9orf72-ALS shared",
      `C9-ALS_brain` & sALS_brain & !`C9-ALS_sc` & !sALS_sc ~ "Motor cortex shared",
      `C9-ALS_sc` & sALS_sc & !`C9-ALS_brain` & !sALS_brain ~ "Cervical spinal cord shared",
      TRUE ~ "Other"
    )
  )

png(file = paste0(plots_dir, "microglia_degs_upset.png"),
    height = 6, width = 7,
    units = "in", res = 600)

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
        geom = geom_point(shape = 19, size = 10),
        segment = geom_segment(linewidth = 1.5),
        outline_color = list(active = "white", inactive = "white")
      )),
      base_annotations = list(
        'Intersection size' = intersection_size(
          text = list(size = 0),
          mapping = aes(fill = highlight_group)
        ) +
          scale_fill_manual(
            values = c(
              "C9orf72-ALS shared" = "#CC00FF",
              "Motor cortex shared" = "#0073C2",
              "Cervical spinal cord shared" = "#EFC000",
              "Other" = "grey35"
            ),
            guide = "none"
          ) +
          scale_y_continuous(expand = c(0, 0)) +
          ylab("# DEGs")
      ),
      queries = list(
        upset_query(intersect = c("C9-ALS_brain", "C9-ALS_sc"),
                    color = "#CC00FF", fill = "#CC00FF",
                    only_components = "intersections_matrix"),
        upset_query(intersect = c("C9-ALS_brain", "sALS_brain"),
                    color = "#0073C2", fill = "#0073C2",
                    only_components = "intersections_matrix"),
        upset_query(intersect = c("C9-ALS_sc", "sALS_sc"),
                    color = "#EFC000", fill = "#EFC000",
                    only_components = "intersections_matrix")
      ),
      theme = upset_modify_themes(
        list(
          'Intersection size' = theme(axis.text = element_text(color = "black", size = 16),
                                      axis.title = element_text(size = 20),
                                      axis.ticks.y = element_line(),
                                      panel.border = element_rect(color = "black", fill = "transparent"),
                                      panel.grid = element_line(color = "gray70")),
          'intersections_matrix' = theme(axis.text = element_text(color = "black", size = 16),
                                         axis.title = element_blank())
        )
      )
) +
  ggtitle("Microglia DEGs") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "plain"))

dev.off()

# Fold-change scatter between any two comparisons ---------------------------
# Flexible, per the user -- a function rather than a fixed block, so any
# two (tissue, contrast) comparisons can be plotted against each other
# without editing the body each time. Uses the same raw (non-LFC-shrunk)
# results CSVs and "X" gene-identifier column as the upset plot above,
# for the same reason (consistency with deseq_viz1.R's established DEG
# definition in this project). Unlike the upset plot, this only shows
# genes significant in at least one of the two comparisons ("plots of
# significant DEGs"), not the whole transcriptome background.
#
# Returns the ggplot object rather than saving it -- since the specific
# pair of comparisons varies by call, there's no single sensible fixed
# output filename to save to; ggsave() it yourself under whatever name
# fits the comparison you ran.

plot_fc_scatter <- function(tissue_1, contrast_1, tissue_2, contrast_2,
                            label_1 = paste(tissue_1, contrast_1),
                            label_2 = paste(tissue_2, contrast_2)){

  path_1 <- paste0(results_dir, tissue_1, "/", celltype, "/", contrast_1, ".csv")
  path_2 <- paste0(results_dir, tissue_2, "/", celltype, "/", contrast_2, ".csv")

  missing <- c(path_1, path_2)[!file.exists(c(path_1, path_2))]
  if (length(missing) > 0){
    stop(paste0("Missing DESeq2 results file(s): ",
                paste(missing, collapse = ", ")))
  }

  res_1 <- read.csv(path_1) %>%
    dplyr::select(X, log2FoldChange, padj)
  res_2 <- read.csv(path_2) %>%
    dplyr::select(X, log2FoldChange, padj)

  df <- inner_join(res_1, res_2, by = "X", suffix = c("_1", "_2"))

  df <- df %>%
    mutate(sig_1 = !is.na(padj_1) & padj_1 < 0.05 & abs(log2FoldChange_1) > log2(1.5),
           sig_2 = !is.na(padj_2) & padj_2 < 0.05 & abs(log2FoldChange_2) > log2(1.5),
           sig_group = case_when(
             sig_1 & sig_2 ~ "Both",
             sig_1 & !sig_2 ~ label_1,
             !sig_1 & sig_2 ~ label_2,
             TRUE ~ "Neither"
           )) %>%
    filter(sig_group != "Neither")

  ggplot(df, aes(x = log2FoldChange_1, y = log2FoldChange_2, color = sig_group)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = paste0("log2(fold change) in ", label_1),
        y = paste0("log2(fold change) in ", label_2),
        color = "DEGs versus\ncontrol in:") +
    theme_linedraw(base_size = 12) +
    theme(axis.text = element_text(color = "black"))
}

plot_fc_scatter("sc", "C9orf72_vs_Control",
                "sc", "sALS_vs_Control",
                label_1 = "C9orf72",
                label_2 = "sALS") + 
  scale_color_manual(values = c("darkslategrey", "#0CAA00", "#CC00FF")) + 
  ggtitle("C9orf72-ALS vs sALS DEGs\nin cervical spinal cord") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(plots_dir, "microglia_sc_c9_vs_sals_multivolcano.png"),
       units = "in", dpi = 600,
       height = 4, width = 4.5)

plot_fc_scatter("brain", "C9orf72_vs_Control",
                "sc", "C9orf72_vs_Control",
                label_2 = "Cervical spinal cord",
                label_1 = "Motor cortex") + 
  scale_color_manual(values = c("darkgreen", "#EFC000", "#0073C2")) + 
  ggtitle("Motor cortex vs cervical spinal\ncord DEGs in C9orf72-ALS") + 
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(plots_dir, "microglia_c9_br_vs_sc_multivolcano.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# Example matching the user's own comparison:
# p <- plot_fc_scatter("brain", "C9orf72_vs_Control", "sc", "C9orf72_vs_Control",
#                      label_1 = "C9-ALS motor cortex", label_2 = "C9-ALS spinal cord")
# p

# GSEA lollipop chart across all 4 comparisons ------------------------------
# Brings in deseq2_gsea.R's output (results/deseq2/<tissue>/<celltype>/
# <contrast>_GSEA.csv) for microglia, in both tissues and both groups.
# GSEA()'s own default pvalueCutoff already restricts what deseq2_gsea.R
# saved to nominally significant pathways, so no additional significance
# filtering is applied here beyond picking the top N per comparison.
# "top_n" pathways are taken per comparison, then the union of those
# across all 4 comparisons is plotted so the same pathway's behavior can
# be compared across tissue/group even if it wasn't top-ranked in every
# comparison.

top_n <- 10 # change as needed

gsea_files <- expand.grid(tissue = c("brain", "sc"),
                          contrast = c("sALS_vs_Control", "C9orf72_vs_Control"),
                          stringsAsFactors = F) %>%
  mutate(path = paste0(results_dir, tissue, "/", celltype, "/", contrast,
                       "_GSEA.csv"),
         group = if_else(contrast == "sALS_vs_Control", "sALS", "C9orf72-ALS"),
         tissue_label = if_else(tissue == "brain", "Motor cortex",
                                "Cervical spinal cord"))

missing <- gsea_files$path[!file.exists(gsea_files$path)]
if (length(missing) > 0){
  stop(paste0("Missing GSEA results file(s) for ", celltype, ": ",
              paste(missing, collapse = ", "),
              " -- check whether deseq2_gsea.R has been run for this ",
              "cell type/tissue/contrast."))
}

gsea_all <- gsea_files %>%
  mutate(data = map(path, read.csv)) %>%
  unnest(data)

top_pathways <- gsea_all %>%
  group_by(tissue, contrast) %>%
  slice_min(p.adjust, n = top_n) %>%
  ungroup() %>%
  pull(ID) %>%
  unique()

plot_df <- gsea_all %>%
  filter(ID %in% top_pathways)

pathway_order <- plot_df %>%
  group_by(ID) %>%
  summarize(mean_abs_nes = mean(abs(NES))) %>%
  arrange(mean_abs_nes) %>%
  pull(ID)

plot_df <- plot_df %>%
  mutate(ID = factor(ID, levels = pathway_order))

p <- ggplot(plot_df, aes(x = NES, y = ID)) +
  geom_segment(aes(xend = 0, yend = ID, color = group)) +
  geom_point(aes(color = group, shape = tissue_label), size = 3) +
  geom_vline(xintercept = 0, color = "grey50") +
  labs(x = "Normalized enrichment score", y = NULL,
      color = "Group", shape = "Tissue",
      title = "Microglia: top dysregulated pathways") +
  theme_bw(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

ggsave(p,
       filename = paste0(plots_dir, "microglia_gsea_lollipop.pdf"),
       height = max(6, length(pathway_order) * 0.3), width = 9)
