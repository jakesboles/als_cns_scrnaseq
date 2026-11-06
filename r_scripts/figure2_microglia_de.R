library(tidyverse)
library(ComplexUpset)
library(ggrastr)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

results_dir <- "results/deseq2/"
plots_dir <- "figures/figure2/"
dir.create(plots_dir, showWarnings = F, recursive = T)

source("r_scripts/plot_theme.R")

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

pdf(file = paste0(plots_dir, "microglia_degs_upset.pdf"),
    height = 2, width = 2.25)

upset(intersect_df, list_names,
      set_sizes = F,
      stripes = "white",
      wrap = T,
      encode_sets = F,
      height_ratio = 0.9,
      labeller = ggplot2::as_labeller(c(
        "C9-ALS_brain" = "C9orf72-ALS MCX",
        "sALS_brain" = "sALS MCX",
        "C9-ALS_sc" = "C9orf72-ALS SC",
        "sALS_sc" = "sALS SC"
      )),
      matrix = (intersection_matrix(
        geom = geom_point(shape = 19, size = 2),
        segment = geom_segment(linewidth = 1),
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
          'Intersection size' = theme_linedraw() + theme_als_project() +
            theme(panel.border = element_rect(color = "black", fill = "transparent"),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  plot.margin = margin(t = 0, r = 0, b = 1, l = 0, unit = "pt")),
          'intersections_matrix' = theme_als_project() +
            theme(axis.text = element_text(color = "black"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  panel.border = element_rect(color = "black", fill = "transparent"),
                  plot.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "pt"))
        )
      )
) +
  ggtitle("Microglia DEGs") +
  theme_als_project(title_size = 7)

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
    ggrastr::rasterize(
      geom_point(size = 0.5),
      dpi = 1600) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = paste0("log2(fold change) in ", label_1),
         y = paste0("log2(fold change) in ", label_2),
         color = "DEGs versus\ncontrol in:") +
    theme_linedraw() +
    theme_als_project(title_size = 7) +
    theme(axis.text = element_text(color = "black"),
          legend.box.spacing = unit(0.1, "cm"),
          legend.title = element_text(margin = margin(b = 0.05, unit = "cm")))
}

plot_fc_scatter("sc", "sALS_vs_Control",
                "sc", "C9orf72_vs_Control",
                label_1 = "sALS",
                label_2 = "C9orf72-ALS") + 
  scale_color_manual(values = c("darkslategrey", "#0CAA00", "#CC00FF"),
                     breaks = c("Both", "sALS", "C9orf72-ALS")) + 
  ggtitle("sALS vs C9orf72-ALS DEGs\nin cervical spinal cord")
save_figure(plot = last_plot(),
            filename = paste0(plots_dir, "multivolcano_sc.pdf"),
            height = 2, width = 2.5)

plot_fc_scatter("brain", "C9orf72_vs_Control",
                "sc", "C9orf72_vs_Control",
                label_2 = "SC",
                label_1 = "MCX") + 
  scale_color_manual(values = c("darkslategrey", "#0073C2", "#EFC000")) + 
  ggtitle("Motor cortex vs cervical spinal\ncord DEGs in C9orf72-ALS")
save_figure(plot = last_plot(),
            filename = paste0(plots_dir, "multivolcano_c9.pdf"),
            height = 2, width = 2.25)
