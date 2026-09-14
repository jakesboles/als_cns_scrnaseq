# Re-clusters the microglia object (data/19_subclustering3/microglia) at
# several resolutions, picks the best one by graph modularity, and
# cross-references the resulting clusters against milo.R's differential
# neighborhood abundance results, so DA-neighborhood content can inform
# manual subcluster labeling before those labels get folded back onto the
# full CNS tissue object -- itself headed for use as the Cell2Location
# reference against Visium data in the same tissues. Interactive script,
# not a SLURM job (no array setup) -- the clustering/Milo sections below
# are meant to run once and be inspected; final label assignment is done
# by hand afterward (not implemented here, per the user).
#
# Design notes:
# - Resolution sweep uses 13_subclustering1.R/15_subclustering2.R's
#   13-value resolution list and graph_modularity() helper (copied
#   verbatim), not 10_clustering.R's 19-value full-tissue list -- this
#   script re-clusters one already-subclustered cell type, the same scope
#   as 13/15, not a whole tissue.
# - RNA_snn already exists from the single FindNeighbors() call above (no
#   return.neighbor = T companion call needed here, unlike 13/15's two-call
#   pattern -- this script reuses 19_subclustering3.R's already-fit UMAP
#   rather than computing a fresh one, so there's no RunUMAP() call that
#   would need a Neighbor object).
# - obj$seurat_clusters is explicitly reset to the best resolution's column
#   after the sweep -- FindClusters() only ever overwrites seurat_clusters
#   with the *last*-tested resolution (5, the top of res_tests), not
#   necessarily the best one, and the DimPlot_scCustom()/dittoBarPlot()
#   calls at the end (left untouched, per the user) reference
#   seurat_clusters directly.
# - Milo cross-reference (design choice, flagged rather than guessed
#   silently -- no prior script in this project cross-references Milo
#   results against a separate clustering, so there's no established
#   pattern to follow): for each cluster and each (tissue, contrast)
#   comparison, computes the fraction of the cluster's cells belonging to
#   at least one significant (SpatialFDR < 0.05) neighborhood, split by
#   direction (enriched, logFC > 0, vs. depleted, logFC < 0). A cluster
#   with a high fraction in either direction is one "more likely to
#   contain Nhoods that are changing in ALS" for that comparison. This is
#   cell-membership-based (via nhoods(milo)'s incidence matrix), not
#   restricted to cells actually in that comparison's tissue -- a
#   neighborhood spans both tissues by construction (see milo.R's header),
#   and its significance in one tissue's test is itself the whole point of
#   sharing neighborhoods across tissues (spotting shared vs.
#   tissue-specific change), so any member cell counts. If a different
#   summary would be more useful once you're looking at it (e.g. mean
#   logFC per cluster regardless of significance, or a direct UMAP
#   overlay instead of a table/bar chart), easy to swap out -- say so.
# - milo.R's Milo object is built from this exact same
#   data/19_subclustering3/microglia source (same script, no extra
#   filtering in either place), so cell barcodes are expected to match
#   directly -- unlike milo_wgcna.R, which had to bridge two different
#   pipelines' microglia populations. Checked at runtime via a coverage
#   message anyway, with a hard stop() only if coverage is implausibly low.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(scater)
  library(igraph)
  library(scales)
  library(BPCells)
  library(dittoSeq)
  library(patchwork)
  library(miloR)
  library(Matrix)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

set.seed(256)

# Assemble object -----------

cell_type <- "microglia"

data_dir <- paste0("data/20_subcluster_phenotyping/", cell_type, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0("results/20_subcluster_phenotyping/", cell_type, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

in_dir <- paste0("data/19_subclustering3/", cell_type)

meta <- readRDS(paste0(in_dir, "/metadata.rds"))
harmony <- readRDS(paste0(in_dir, "/harmony.rds"))
umap <- readRDS(paste0(in_dir, "/harmony_umap.rds"))

mat <- open_matrix_dir(paste0(in_dir, "/bpcells_data"))

obj <- CreateSeuratObject(counts = mat,
                          meta.data = meta)
obj[["RNA"]]$data <- mat
obj[["harmony"]] <- harmony
obj[["umap"]] <- umap

obj$group <- factor(obj$group, levels = c("Control", "sALS", "C9orf72"))

# Find neighbors and cluster ----------------------------------------------

obj <- FindNeighbors(obj,
                     reduction = "harmony",
                     dims = 1:15,
                     k.param = 15,
                     nn.method = "annoy",
                     annoy.metric = "euclidean",
                     compute.SNN = T)

# Cluster at several resolutions, scoring each by graph modularity ----------
# Copied from 13_subclustering1.R -- see header note above.

graph_modularity <- function(obj, clusters, graph_name = "RNA_snn"){
  snn <- as(obj[[graph_name]], "dgCMatrix")
  g <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::modularity(g, membership = as.integer(factor(clusters)))
}

res_tests <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5)

modularity_vec <- numeric(length(res_tests))

for (i in seq_along(res_tests)){
  res <- res_tests[i]
  message2(paste0("Clustering at resolution = ", res))

  obj <- FindClusters(obj,
                      resolution = res,
                      algorithm = 4,
                      graph.name = "RNA_snn",
                      cluster.name = paste0("res", res, "_clusters"),
                      method = "igraph")

  modularity_vec[i] <- graph_modularity(obj,
                                        obj@meta.data[[paste0("res", res, "_clusters")]])
}

message2("Saving graph modularity table and plot")

modularity_df <- data.frame(resolution = res_tests, modularity = modularity_vec)

write.csv(modularity_df,
          file = paste0(results_dir, "graph_modularity.csv"),
          row.names = F)

p <- ggplot(modularity_df, aes(x = resolution, y = modularity)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(p,
       filename = paste0(results_dir, "graph_modularity.png"),
       units = "in", dpi = 300,
       height = 4, width = 6)

# Pick the resolution with the highest graph modularity ----------------------

best_res <- res_tests[which.max(modularity_vec)]
best_res_col <- paste0("res", best_res, "_clusters")

message2(paste0("Best resolution = ", best_res))

Idents(obj) <- best_res_col
# seurat_clusters is reset here since FindClusters() left it pointing at
# the last-tested resolution, not necessarily the best one -- see header
# note above.
obj$seurat_clusters <- Idents(obj)

# Bring in Milo differential neighborhood abundance results -----------------
# See header note above for the full design rationale.

message2("Reading in Milo object and differential abundance results")

milo <- readRDS("data/milo/microglia/milo_obj.rds")

nh_mat <- nhoods(milo)
if (is.null(rownames(nh_mat))) rownames(nh_mat) <- colnames(milo)

coverage <- length(intersect(rownames(nh_mat), colnames(obj))) / ncol(obj)
message2(paste0(round(coverage * 100, 1), "% of this object's cells are ",
                "found in the Milo object"))
if (coverage < 0.5){
  stop("Fewer than 50% of this object's cells were found in ",
       "data/milo/microglia/milo_obj.rds -- the assumption that both come ",
       "from the same, unfiltered data/19_subclustering3/microglia source ",
       "appears to be wrong. Check whether either has been rebuilt since ",
       "the other last ran.")
}

combos <- tribble(
  ~tissue_file,           ~tissue_short, ~contrast,            ~contrast_short,
  "Motor_cortex",         "brain",       "sALS_vs_Control",    "sALS",
  "Motor_cortex",         "brain",       "C9orf72_vs_Control", "C9",
  "Cervical_spinal_cord", "sc",          "sALS_vs_Control",    "sALS",
  "Cervical_spinal_cord", "sc",          "C9orf72_vs_Control", "C9"
) %>%
  mutate(comparison = paste0(tissue_short, "_", contrast_short))

clusters <- obj@meta.data[[best_res_col]]
names(clusters) <- colnames(obj)

milo_summary <- list()

for (i in seq_len(nrow(combos))){

  comparison <- combos$comparison[i]

  results_csv <- paste0("results/milo/microglia/", combos$contrast[i], "_",
                        combos$tissue_file[i], "_nhood_results.csv")

  if (!file.exists(results_csv)){
    message2(paste0("Missing ", results_csv, " -- skipping ", comparison))
    next
  }

  message2(paste0("Comparison: ", comparison))

  res <- read.csv(results_csv)
  res$SpatialFDR[is.na(res$SpatialFDR)] <- 1

  sig_up <- res$Nhood[res$SpatialFDR < 0.05 & res$logFC > 0]
  sig_down <- res$Nhood[res$SpatialFDR < 0.05 & res$logFC < 0]

  cell_up <- if (length(sig_up) > 0){
    Matrix::rowSums(nh_mat[, sig_up, drop = F]) > 0
  } else {
    rep(FALSE, nrow(nh_mat))
  }
  cell_down <- if (length(sig_down) > 0){
    Matrix::rowSums(nh_mat[, sig_down, drop = F]) > 0
  } else {
    rep(FALSE, nrow(nh_mat))
  }
  names(cell_up) <- rownames(nh_mat)
  names(cell_down) <- rownames(nh_mat)

  cell_df <- data.frame(
    cell = colnames(obj),
    cluster = clusters[colnames(obj)],
    in_sig_up = cell_up[colnames(obj)],
    in_sig_down = cell_down[colnames(obj)]
  )

  summary_i <- cell_df %>%
    group_by(cluster) %>%
    summarise(frac_enriched = mean(in_sig_up, na.rm = T),
             frac_depleted = mean(in_sig_down, na.rm = T),
             n_cells = n()) %>%
    mutate(comparison = comparison)

  milo_summary[[comparison]] <- summary_i
}

milo_summary_df <- list_rbind(milo_summary)

message2("Saving cluster/Milo summary table and plot")

write.csv(milo_summary_df,
          file = paste0(results_dir, "cluster_milo_summary.csv"),
          row.names = F)

p <- milo_summary_df %>%
  pivot_longer(c(frac_enriched, frac_depleted),
              names_to = "direction", values_to = "fraction") %>%
  mutate(direction = factor(direction,
                            levels = c("frac_enriched", "frac_depleted"),
                            labels = c("Enriched", "Depleted"))) %>%
  ggplot(aes(x = factor(cluster), y = fraction, fill = direction)) +
  geom_col(position = "dodge") +
  facet_wrap(. ~ comparison) +
  labs(x = "Cluster",
      y = "Fraction of cluster's cells in a\nsignificant (SpatialFDR < 0.05) neighborhood",
      fill = "Direction") +
  scale_fill_manual(values = c("Enriched" = "firebrick", "Depleted" = "steelblue")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(p,
       filename = paste0(results_dir, "cluster_milo_summary.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

p

DimPlot_scCustom(obj,
                 label = F) +

dittoBarPlot(obj,
             var = c("group"),
             group.by = "seurat_clusters")
