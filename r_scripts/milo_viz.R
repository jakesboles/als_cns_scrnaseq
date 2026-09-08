# Assembles the plotting inputs for milo.R's differential neighborhood
# abundance results as UMAPs colored by log fold-change -- one entry per
# (tissue, disease-group contrast) combination available for a given
# 19_subclustering3.R target, so all of them can be plotted together
# (e.g. with patchwork) and compared side by side. Interactive script,
# not a SLURM job -- run by hand, changing target_name as needed.
#
# The neighborhoods/graph/embedding are the same shared, cross-tissue
# ones milo.R built once per target -- only the logFC coloring changes
# per entry, since that's what actually varies per tissue's testNhoods()
# result. This lets the same neighborhood, in the same spot on the plot,
# be visually compared across tissues/groups to see whether an abundance
# change is shared or tissue-specific, matching milo.R's own design
# rationale.
#
# This only assembles plotting inputs (plot_data below) -- it does not
# build or save any ggplot object. Each plot_data[[...]] entry has
# everything needed to build one panel (layout, min_logfc, max_logfc,
# breaks, tissue, contrast) for however you want to lay it out.
#
# Design notes (carried over from the earlier array-job version of this
# script, still relevant):
# - milo.R's saved Milo object has no UMAP reduction attached (it only
#   attaches "harmony" for neighborhood construction) --
#   19_subclustering3.R's own saved harmony_umap.rds is loaded and
#   attached here instead.
# - Only combinations with an actual results CSV are included --
#   milo.R's min_donors_per_group check can skip a tissue entirely, so a
#   missing file is a real possibility, not a typo.
# - `size` (per-neighborhood cell count) doesn't depend on tissue/
#   contrast, so it's computed once, not per combination.
# - min_logfc/max_logfc/breaks are computed per combination, since the
#   logFC range differs per test.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(miloR)
  library(scater)
  library(igraph)
  library(scales)
  library(ComplexUpset)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

target_name <- "microglia" # change this to switch targets

data_dir <- paste0("data/milo/", target_name, "/")

results_dir <- paste0("results/milo/", target_name, "/")

# Load the Milo object and reattach a UMAP for plotting ---------------------
# milo.R's saved Milo object has no UMAP reduction -- see header note
# above.

message("Reading in Milo object and UMAP embedding")

milo <- readRDS(paste0(data_dir, "milo_obj.rds"))

umap <- readRDS(paste0("data/19_subclustering3/", target_name,
                       "/harmony_umap.rds"))
reducedDim(milo, "harmony_umap") <- Embeddings(umap)[colnames(milo), ]

# Neighborhood size doesn't depend on tissue/contrast -- computed once.

colData(milo)["size"] <- NA
colData(milo)[unlist(nhoodIndex(milo)), "size"] <-
  as.numeric(vertex_attr(nhoodGraph(milo))[["size"]])

tissues_present <- sort(unique(as.character(colData(milo)$tissue)))
contrasts <- c("sALS_vs_Control", "C9orf72_vs_Control")

# Assemble one plotting-input entry per (tissue, contrast) -----------------

plot_data <- list()

for (tissue_title in tissues_present){

  tissue_file <- str_replace_all(tissue_title, " ", "_")

  for (contrast in contrasts){

    results_csv <- paste0("results/milo/", target_name, "/", contrast, "_",
                          tissue_file, "_nhood_results.csv")

    if (!file.exists(results_csv)){
      message(paste0("Skipping ", target_name, " -- ", contrast, " (",
                     tissue_title, ") -- no results CSV found (likely ",
                     "skipped by milo.R's min_donors_per_group check)."))
      next
    }

    message(paste0("Assembling ", target_name, " -- ", contrast, " (",
                   tissue_title, ")"))

    signif_res <- read.csv(results_csv)

    # Handle untested neighborhoods and zero out non-significant logFC.
    signif_res$SpatialFDR[is.na(signif_res$SpatialFDR)] <- 1
    signif_res[signif_res$SpatialFDR >= 0.05, "logFC"] <- 0

    milo2 <- milo
    colData(milo2)["logFC"] <- NA
    colData(milo2)[unlist(nhoodIndex(milo2)[signif_res$Nhood]), "logFC"] <-
      signif_res[, "logFC"]

    nhoodGraph(milo2) <- set_vertex_attr(nhoodGraph(milo2), name = "logFC",
                                         value = signif_res[, "logFC"])

    nh_graph <- nhoodGraph(milo2)
    nh_graph <- permute(nh_graph,
                        order(abs(vertex_attr(nh_graph)[["logFC"]]),
                              decreasing = TRUE))

    layout <- reducedDim(milo2, "harmony_umap")[as.numeric(vertex_attr(nh_graph)$name), ] %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell")

    coldata <- colData(milo2) %>%
      as.data.frame() %>%
      rownames_to_column(var = "cell")

    layout <- layout %>%
      left_join(coldata %>% dplyr::select(c(cell, logFC, size)), by = "cell")

    max_logfc <- max(signif_res$logFC)
    min_logfc <- min(signif_res$logFC)
    breaks <- c(min_logfc, min_logfc / 2, 0, max_logfc / 2, max_logfc)

    combo_name <- paste0(tissue_file, "_", contrast)

    plot_data[[combo_name]] <- list(
      tissue = tissue_title,
      contrast = contrast,
      layout = layout,
      min_logfc = min_logfc,
      max_logfc = max_logfc,
      breaks = breaks
    )

  }
}

# plot_data now has one entry per available (tissue, contrast)
# combination, e.g. plot_data$Motor_cortex_sALS_vs_Control -- each with
# $layout (harmonyumap_1/2, logFC, size), $min_logfc/$max_logfc/$breaks
# for scale_color_gradientn(), and $tissue/$contrast for titling/
# faceting.

min_logfc <- vector(mode = "numeric", length = 4)
max_logfc <- vector(mode = "numeric", length = 4)

for (i in seq_along(plot_data)){
  min_logfc[i] <- plot_data[[i]]$min_logfc
  max_logfc[i] <- plot_data[[i]]$max_logfc
}

min_logfc <- min(min_logfc)
max_logfc <- max(max_logfc)

breaks <- c(min_logfc, min_logfc/2, 0 , max_logfc/2, max_logfc)

cols <- c(
  "#7F0000",  # dark red
  "#FF3030",  # bright red
  "gray80",   # zero
  "#268BFF",  # bright blue
  "#08306B"   # dark blue
)

p_list <- list()

for (i in seq_along(plot_data)){
  p_list[[i]] <- ggplot(plot_data[[i]]$layout,
                   aes(x = harmonyumap_1,
                       y = harmonyumap_2)) + 
    geom_point(aes(color = logFC),
               size = 0.75) + 
    scale_color_gradientn(
      colours = cols,
      values = rescale(breaks, from = c(min_logfc, max_logfc)),
      limits = c(min_logfc, max_logfc),
      oob = squish) + 
    labs(y = "UMAP 2",
         x = "UMAP 1",
         color = "log2FC\nin ALS") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank())
  
}

wrap_plots(p_list[[3]], p_list[[4]], p_list[[1]], p_list[[2]],
           ncol = 2,
           guides = "collect")
ggsave(filename = paste0(results_dir, "milo_lfc_umaps.png"),
       units = "in", dpi = 600,
       height = 8, width = 9)

for (i in seq_along(plot_data)){
  message(names(plot_data)[i])
}

# Upset plot of differentially abundant (DA) neighborhoods ------------------
# Same 4-category design as deseq_viz2.R's microglia DEG upset plot,
# applied to Milo's DA neighborhoods instead of DESeq2's DEGs: sALS vs.
# Control and C9orf72 vs. Control, each in both tissues.
#
# Neighborhoods are identified by their index cell's barcode (plot_data's
# "cell" column), not re-read from the results CSVs -- since the
# neighborhoods/graph were built once on the shared cross-tissue object
# (see milo.R and this script's own header), the same index-cell barcode
# identifies the same neighborhood in every combination, exactly like a
# gene symbol identifies the same gene across DESeq2 comparisons. A
# neighborhood counts as "DA" in a given combination if its logFC wasn't
# already zeroed out above (SpatialFDR >= 0.05).
#
# Hardcoded to the 4 combinations available for the microglia target (2
# tissues x 2 contrasts), matching deseq_viz2.R's own hardcoding of
# "brain"/"sc" -- fails fast with a clear message if target_name is
# switched to something with a different tissue/contrast layout, rather
# than silently building a mismatched plot.

expected_combos <- c("Motor_cortex_sALS_vs_Control", "Motor_cortex_C9orf72_vs_Control",
                     "Cervical_spinal_cord_sALS_vs_Control", "Cervical_spinal_cord_C9orf72_vs_Control")

if (!setequal(names(plot_data), expected_combos)){
  stop("The DA neighborhood upset plot section assumes exactly the 4 ",
       "combinations available for the microglia target (2 tissues x 2 ",
       "contrasts), matching deseq_viz2.R's DEG upset plot -- got: ",
       paste(names(plot_data), collapse = ", "),
       ". Update the set names/labeller/highlight_group logic below if ",
       "you're running a different target_name.")
}

nhood_lists <- list()
for (combo_name in names(plot_data)){
  nhood_lists[[combo_name]] <- plot_data[[combo_name]]$layout %>%
    filter(logFC != 0) %>%
    pull(cell)
}

nhoods <- unique(list_c(nhood_lists))

intersect_df <- data.frame(nhood = nhoods)
for (combo_name in names(nhood_lists)){
  intersect_df <- intersect_df %>%
    mutate(x = if_else(nhood %in% nhood_lists[[combo_name]], T, F)) %>%
    dplyr::rename(!!sym(combo_name) := "x")
}

list_names <- names(plot_data)

intersect_df <- intersect_df %>%
  mutate(
    highlight_group = case_when(
      Motor_cortex_C9orf72_vs_Control & Cervical_spinal_cord_C9orf72_vs_Control &
        !Motor_cortex_sALS_vs_Control & !Cervical_spinal_cord_sALS_vs_Control ~ "C9orf72-ALS shared",
      Motor_cortex_C9orf72_vs_Control & Motor_cortex_sALS_vs_Control &
        !Cervical_spinal_cord_C9orf72_vs_Control & !Cervical_spinal_cord_sALS_vs_Control ~ "Motor cortex shared",
      Cervical_spinal_cord_C9orf72_vs_Control & Cervical_spinal_cord_sALS_vs_Control &
        !Motor_cortex_C9orf72_vs_Control & !Motor_cortex_sALS_vs_Control ~ "Cervical spinal cord shared",
      TRUE ~ "Other"
    )
  )

png(file = paste0(results_dir, "milo_da_nhoods_upset.png"),
    height = 6, width = 6,
    units = "in", res = 600)

upset(intersect_df, list_names,
      set_sizes = F,
      stripes = "white",
      wrap = T,
      encode_sets = F,
      height_ratio = 0.9,
      labeller = ggplot2::as_labeller(c(
        "Motor_cortex_C9orf72_vs_Control" = "C9orf72-ALS\nMotor cortex",
        "Motor_cortex_sALS_vs_Control" = "sALS\nMotor cortex",
        "Cervical_spinal_cord_C9orf72_vs_Control" = "C9orf72-ALS\nCervical spinal cord",
        "Cervical_spinal_cord_sALS_vs_Control" = "sALS\nCervical spinal cord"
      )),
      matrix = (intersection_matrix(
        geom = geom_point(shape = 19, size = 6),
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
          ylab("# DA neighborhoods")
      ),
      queries = list(
        upset_query(intersect = c("Motor_cortex_C9orf72_vs_Control", "Cervical_spinal_cord_C9orf72_vs_Control"),
                    color = "#CC00FF", fill = "#CC00FF",
                    only_components = "intersections_matrix"),
        upset_query(intersect = c("Motor_cortex_C9orf72_vs_Control", "Motor_cortex_sALS_vs_Control"),
                    color = "#0073C2", fill = "#0073C2",
                    only_components = "intersections_matrix"),
        upset_query(intersect = c("Cervical_spinal_cord_C9orf72_vs_Control", "Cervical_spinal_cord_sALS_vs_Control"),
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
  ggtitle("Microglia DA neighborhoods") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "plain"))

dev.off()
