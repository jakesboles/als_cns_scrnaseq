# Clusters the microglia object (data/19_subclustering3/microglia) at a
# fixed resolution and cross-references the resulting clusters against
# milo.R's differential neighborhood abundance results, so DA-neighborhood
# content can inform manual subcluster labeling before those labels get
# folded back onto the full CNS tissue object -- itself headed for use as
# the Cell2Location reference against Visium data in the same tissues.
# Interactive script, not a SLURM job (no array setup) -- the clustering/
# Milo sections below are meant to run once and be inspected; final label
# assignment is done by hand afterward (not implemented here, per the
# user).
#
# Design notes:
# - Resolution is fixed at 0.8, per the user -- an earlier version of this
#   script swept 13_subclustering1.R/15_subclustering2.R's 13-value
#   resolution list and picked the best by graph modularity, but the user
#   reviewed that output and settled on 0.8 directly, so the sweep is
#   gone.
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
#   tissue-specific change), so any member cell counts.
# - Second Milo view: a jittered dot plot with one dot per neighborhood
#   (not per cell, unlike the bar chart above), x = the cluster of that
#   neighborhood's own index cell, y = logFC, colored by significance.
#   "The cluster a neighborhood belongs to" is necessarily an
#   approximation -- a neighborhood is ~15 cells around an index cell and
#   can span a cluster boundary -- so this uses the same index-cell
#   convention milo_viz.R/milo_wgcna.R already established for
#   representing one neighborhood by one cell, rather than inventing a
#   new one.
# - Note for later: SpatialFDR NA is coerced to 1 (-> not significant)
#   before splitting into sig_up/sig_down for the bar chart, matching
#   milo_viz.R's convention -- milo_wgcna.R instead preserves NA
#   (untested/degenerate neighborhoods, e.g. ones barely present in a
#   given tissue -- see its own header) as NA rather than collapsing to
#   FALSE. Left as-is here since the bar chart code was to be kept
#   unchanged, but worth knowing the new dot plot will show those
#   degenerate neighborhoods as ordinary non-significant logFC == 0 points
#   rather than flagging them as untested.
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

# Cluster at a fixed resolution -----------------------------------------
# resolution = 0.8, per the user -- see header note above.

obj <- FindClusters(obj,
                    algorithm = 4,
                    method = "igraph",
                    resolution = 0.2,
                    cluster.name = "cluster")

p <- DimPlot_scCustom(obj,
                 reduction = "umap",
                 label = F)
ggsave(p,
       filename = paste0(results_dir, "raw_cluster_umap.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)
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

clusters <- obj$seurat_clusters
names(clusters) <- colnames(obj)

# Neighborhood -> cluster lookup for the dot plot below -- one row per
# neighborhood, cluster of its own index cell (see header note above).
nhood_index_rows <- unlist(nhoodIndex(milo))
nhood_cells <- colnames(milo)[nhood_index_rows]
nhood_cluster_lookup <- data.frame(
  Nhood = seq_along(nhood_index_rows),
  cluster = as.character(clusters[nhood_cells])
)

milo_summary <- list()
nhood_level <- list()

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

  # One row per neighborhood (not per cell) for the dot plot below --
  # reuses this same res/sig computation rather than re-reading the CSV.
  nhood_level[[comparison]] <- res %>%
    mutate(sig = SpatialFDR < 0.05) %>%
    dplyr::select(Nhood, logFC, sig) %>%
    left_join(nhood_cluster_lookup, by = "Nhood") %>%
    mutate(comparison = comparison)
}

milo_summary_df <- list_rbind(milo_summary)
nhood_level_df <- list_rbind(nhood_level)

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

# Jittered dot plot: one dot per neighborhood, x = cluster of its index
# cell, y = logFC -- see header note above.

p2 <- nhood_level_df %>%
  ggplot(aes(x = cluster, y = logFC, color = sig)) +
  geom_jitter(width = 0.2, height = 0, alpha = 0.6) +
  facet_wrap(. ~ comparison) +
  labs(x = "Cluster (of neighborhood's index cell)", y = "logFC",
      color = "Significant\n(SpatialFDR < 0.05)") +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(p2,
       filename = paste0(results_dir, "cluster_milo_nhood_dotplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

p2

# FAM ---------------------------------------------------------------------

markers <- FindAllMarkers(obj)

write.csv(markers,
          file = paste0(results_dir, "cluster_markers.csv"))

# Rewrite labels and add to full meta data --------------------------------

obj@meta.data <- obj@meta.data %>% 
  mutate(cell_type4 = paste0("Microglia", cluster))

DimPlot_scCustom(obj,
                 group.by = "cell_type4")

full_meta <- readRDS("data/18_full_integration/brain_sc/metadata.rds")

microglia_meta <- obj@meta.data

colnames(full_meta)
colnames(microglia_meta)

cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "tissue", "batch", "group", "id",
          "percent_mito", "log10GenesPerUMI", "cell_type3", "cell_type4")

full_meta <- full_meta[, colnames(full_meta) %in% cols]
microglia_meta <- microglia_meta[, colnames(microglia_meta) %in% cols]

full_meta$cell_type4 <- full_meta$cell_type3

idx <- match(rownames(microglia_meta), rownames(full_meta))

full_meta$cell_type4[idx] <- microglia_meta$cell_type4

unique(full_meta$cell_type4)

table(full_meta$cell_type3, full_meta$cell_type4)

saveRDS(full_meta,
        file = paste0(data_dir, "full_metadata.rds"))
saveRDS(microglia_meta,
        file = paste0(data_dir, "microglia_metadata.rds"))

# Extra plotting  ---------------------------------------------------------

# if running after the initial run of the above script:
markers <- read.csv(paste0(results_dir, "cluster_markers.csv"))
new_meta <- readRDS(paste0(data_dir, "microglia_metadata.rds"))
obj <- AddMetaData(obj, new_meta)
# obj$group <- factor(obj$group,
#                     levels = c("Control", "sALS", "C9orf72"),
#                     labels = c("Control", "sALS", "C9orf72-ALS"))

DimPlot_scCustom(obj,
                 group.by = "cell_type4",
                 pt.size = 1, # comment out if running smaller object
                 colors_use = paletteer_d("ggsci::default_locuszoom")[c(1:3, 5:6)]) + 
  ggtitle("Microglia sub-cluster") + 
  guides(color = guide_legend(ncol = 1,
                              override.aes = list(size = 4)))
ggsave(filename = paste0(results_dir, "subcluster_dimplot.png"),
       units = "in", dpi = 300,
       height = 4, width = 5)

top <- markers %>% 
  filter(pct.1 > 0.3 & avg_log2FC > 0) %>% 
  Extract_Top_Markers(num_features = 10,
                      make_unique = F,
                      named_vector = F)

markers %>% 
  filter(cluster == 4 & 
           pct.1 > 0.3) %>% 
  arrange(desc(avg_log2FC)) %>% head(30)

dittoDotPlot(obj,
             vars = top,
             group.by = "cell_type4") + 
  scale_y_discrete(limits = rev) + 
  scale_color_gradient2() +
  labs(y = "Subtype",
       size = "pct",
       color = "exp")
ggsave(filename = paste0(results_dir, target, "/marker_dotplot.png"),
       units = "in", dpi = 600,
       height = 4, width = 15)

br_cells <- obj@meta.data %>% 
  filter(tissue == "Motor cortex") %>% 
  rownames()

sc_cells <- obj@meta.data %>% 
  filter(tissue == "Cervical spinal cord") %>% 
  rownames()

p1 <- dittoBarPlot(obj,
             group.by = "cell_type4",
             var = "group",
             scale = "count",
             cells.use = br_cells,
             var.labels.reorder = c(2, 3, 1),
             color.panel = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Motor cortex") + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

p2 <- dittoBarPlot(obj,
             group.by = "cell_type4",
             var = "group",
             scale = "count",
             cells.use = sc_cells,
             var.labels.reorder = c(2, 3, 1),
             color.panel = c("#b8b0a8", "#0CAA00", "#CC00FF")) + 
  scale_y_continuous(expand = c(0, 0)) +
  ggtitle("Cervical spinal cord") + 
  theme(axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5))

p1 + p2 + 
  plot_layout(ncol = 1,
              guides = "collect")
ggsave(filename = paste0(results_dir, "subcluster_by_group_by_tissue_bars.png"),
       units = "in", dpi = 600,
       height = 8, width = 5)

modules <- read.csv("results/wgcna_consensus/Microglia/module_scores_ucell.csv")
modules <- modules %>% 
  column_to_rownames(var = "X")

obj <- AddMetaData(obj, 
                   modules)

VlnPlot_scCustom(obj,
                 features = c("blue_UCell_kNN", "turquoise_UCell_kNN"),
                 group.by = "cell_type4",
                 num_columns = 1,
                 colors_use = paletteer_d("ggsci::default_locuszoom")[c(1:3, 5:6)]) + 
  theme(axis.title.x = element_blank())
ggsave(filename = paste0(results_dir, "subcluster_wgcna_module_scores_vln.png"),
       units = "in", dpi = 600,
       height = 8, width = 4)