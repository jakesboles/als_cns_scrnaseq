# Visualizes milo.R's differential neighborhood abundance results as
# UMAPs colored by log fold-change, one plot per (tissue, disease-group
# contrast) combination for each 19_subclustering3.R target. The
# neighborhoods/graph/embedding plotted are the same shared, cross-tissue
# ones milo.R built once per target -- only the logFC coloring changes
# per plot, since that's what actually varies per tissue's testNhoods()
# result. This is deliberate: the same neighborhood, in the same spot on
# the plot, can be visually compared across tissues/groups to see whether
# an abundance change is shared or tissue-specific, matching milo.R's own
# design rationale. Runs as a SLURM job array (see jobs/milo_viz.sh), one
# task per subdirectory of data/milo/ -- discovered at runtime, matching
# milo.R's own convention.
#
# Reworked from the user's pushed code chunk (meant to run directly after
# milo.R's own code, still in the same R session -- milo/c9_results were
# assumed already in memory). Design notes:
# - Loads milo.R's saved data/milo/<target>/milo_obj.rds directly, since
#   this runs as its own script/session.
# - milo.R's saved Milo object has no UMAP reduction attached (it only
#   attaches "harmony" for neighborhood construction -- see milo.R's
#   header). 19_subclustering3.R's own saved harmony_umap.rds is loaded
#   and attached here instead of the old project's "INTEGRATED_UMAP2",
#   which doesn't exist in this pipeline.
# - Loops over every (tissue, contrast) combination actually present in
#   results/milo/<target>/ (skipped, not fatal, if missing -- milo.R's
#   min_donors_per_group check can skip a tissue entirely), covering both
#   contrasts and every tissue found -- the user's pushed chunk only
#   demonstrated this for one contrast (c9_results) in what was implicitly
#   a single tissue, per their "tissue- and disease-group-specific"
#   request.
# - `size` (per-neighborhood cell count) was computed in the pushed chunk
#   but never actually used in the final ggplot() call -- added
#   aes(size = size) to the plot here, since computing it only makes
#   sense if it's meant to be plotted; flagging this rather than silently
#   leaving it as dead code.
# - min_logfc/max_logfc/breaks are recomputed for each (tissue, contrast)
#   plot, not once globally, since the logFC range differs per test.
# - `size` itself doesn't depend on tissue/contrast (same neighborhoods,
#   same sizes throughout), so it's computed once per target rather than
#   recomputed in the inner loop.
# - Output goes to results/milo_viz/<target>/, matching this project's
#   usual per-script results/ convention, not figures/ -- this produces
#   several exploratory plots per target (up to 3 tissues x 2 contrasts),
#   not one curated figure like deseq_viz2.R's upset plot.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(miloR)
  library(scater)
  library(igraph)
  library(scales)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

project_root <- "/projects/b1169/boles/als_cns_scrnaseq"
setwd(project_root)

# Figure out which target this task handles ------------------------------
# Discovered from data/milo/ (milo.R's own output), matching milo.R's own
# convention.

targets <- sort(list.dirs("data/milo", recursive = F, full.names = F))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/milo_viz.sh), one task per subdirectory ",
       "of data/milo/, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(targets), " targets found in data/milo/ -- check the ",
              "--array range in jobs/milo_viz.sh."))
}

target_name <- targets[task_id]

message2(paste0("Processing ", target_name, ", task ", task_id, "/",
                length(targets)))

data_dir <- paste0("data/milo/", target_name, "/")

results_dir <- paste0(project_root, "/results/milo_viz/", target_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Load the Milo object and reattach a UMAP for plotting ---------------------
# milo.R's saved Milo object has no UMAP reduction -- see header note
# above.

message2("Reading in Milo object and UMAP embedding")

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

cols <- c(
  "#7F0000",  # dark red
  "#FF3030",  # bright red
  "gray80",   # zero
  "#268BFF",  # bright blue
  "#08306B"   # dark blue
)

for (tissue_title in tissues_present){

  tissue_file <- str_replace_all(tissue_title, " ", "_")

  for (contrast in contrasts){

    results_csv <- paste0(project_root, "/results/milo/", target_name, "/",
                          contrast, "_", tissue_file, "_nhood_results.csv")

    if (!file.exists(results_csv)){
      message2(paste0("Skipping ", target_name, " -- ", contrast, " (",
                      tissue_title, ") -- no results CSV found (likely ",
                      "skipped by milo.R's min_donors_per_group check)."))
      next
    }

    message2(paste0("Plotting ", target_name, " -- ", contrast, " (",
                    tissue_title, ")"))

    signif_res <- read.csv(results_csv)

    # Handle untested neighborhoods and zero out non-significant logFC --
    # same as the user's pushed chunk.
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

    p <- ggplot(layout,
               aes(x = harmonyumap_1,
                   y = harmonyumap_2)) +
      geom_point(aes(color = logFC, size = size)) +
      scale_color_gradientn(
        colours = cols,
        values = rescale(breaks, from = c(min_logfc, max_logfc)),
        limits = c(min_logfc, max_logfc),
        oob = squish
      ) +
      ggtitle(paste0(target_name, ": ", contrast, " (", tissue_title, ")")) +
      theme_bw(base_size = 14) +
      theme(axis.text = element_text(color = "black"),
            plot.title = element_text(hjust = 0.5))

    ggsave(p,
           filename = paste0(results_dir, contrast, "_", tissue_file,
                            "_nhood_umap.png"),
           units = "in", dpi = 300,
           height = 6, width = 7)

  }
}
