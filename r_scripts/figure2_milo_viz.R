# Visualizes milo.R's differential neighborhood abundance results for a
# given 19_subclustering3.R target: a shared-scale UMAP panel colored by
# logFC for each (tissue, contrast) combination, and an upset plot of
# which combinations' differentially abundant neighborhoods overlap.
# Interactive script, not a SLURM job -- run by hand, changing
# target_name as needed.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(miloR)
  library(scater)
  library(igraph)
  library(scales)
  library(ComplexUpset)
  library(ggrastr)
  library(patchwork)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

source("r_scripts/plot_theme.R")

target_name <- "microglia" # change this to switch targets

data_dir <- paste0("data/milo/", target_name, "/")

results_dir <- "figures/figure2/"

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
  p_list[[i]] <- plot_data[[i]]$layout %>% 
    arrange(abs(logFC)) %>%
    ggplot(aes(x = harmonyumap_1,
                            y = harmonyumap_2)) + 
      ggrastr::rasterize(
        geom_point(aes(color = logFC),
                 size = 0.1),
        dpi = 1600) + 
    scale_color_gradientn(
      colours = cols,
      values = rescale(breaks, from = c(min_logfc, max_logfc)),
      limits = c(min_logfc, max_logfc),
      oob = squish) + 
    labs(y = "UMAP 2",
         x = "UMAP 1",
         color = "Neighborhood log2FC") +
    umap_labs + 
    theme_void() +
    theme_als_project() + 
    umap_theme
  
}

wrap_plots(p_list[[3]], p_list[[4]], p_list[[1]], p_list[[2]],
           nrow = 2,
           guides = "collect") & 
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.spacing = unit(0.1, "cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        theme(plot.margin = margin(0, 0, 0, 0, "pt")),
        legend.text = element_text(margin = margin(t = 0.02, unit = "cm")))
save_figure(last_plot(),
            filename = paste0(results_dir, "milo_umap.pdf"),
            height = 2.2, width = 2)

for (i in seq_along(plot_data)){
  message(names(plot_data)[i])
}

# Upset plot of diff abundant Nhoods --------------------------------------

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

pdf(file = paste0(results_dir, "microglia_nhood_upset.pdf"),
    height = 2, width = 2.25)

upset(intersect_df, list_names,
      set_sizes = F,
      stripes = "white",
      wrap = T,
      encode_sets = F,
      height_ratio = 0.9,
      labeller = ggplot2::as_labeller(c(
        "Motor_cortex_C9orf72_vs_Control" = "C9orf72-ALS MCX",
        "Motor_cortex_sALS_vs_Control" = "sALS MCX",
        "Cervical_spinal_cord_C9orf72_vs_Control" = "C9orf72-ALS SC",
        "Cervical_spinal_cord_sALS_vs_Control" = "sALS SC"
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
  ggtitle("Microglia DA neighborhoods") +
  theme_als_project(title_size = 7)

dev.off()
