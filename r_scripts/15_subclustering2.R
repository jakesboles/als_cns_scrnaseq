# Two things happen here, both per task:
#
# 1. Rebuilds the full tissue object, applies the round-1 annotation
#    (cell_type1, from results/12_annotation1/), and folds in each cell
#    type's round-2 annotation (cell_type2, from 14_findmarkers2.R's
#    per-cell-type annotations.csv) -- unchanged from the original
#    15_annotation2.R, saved to results/15_annotation2/<tissue>_*_labels.png.
# 2. Subsets that object to one of the cell_type2 groups below and
#    re-processes it from scratch, very similarly to 13_subclustering1.R:
#    fresh PCA/Harmony/neighbors/UMAP/clustering, picks the resolution
#    with the highest graph modularity, runs FindAllMarkers(), and saves
#    markers/a top-5 dot plot/a cluster DimPlot plus the re-processed
#    object (metadata, normalized expression, Harmony embedding, UMAP)
#    to data|results/15_subclustering2/<tissue>/<group>/. Also
#    saves a DimPlot of the input cell_type2 labels on the new embedding,
#    before clustering -- most useful for the multi-cell-type groups, to
#    see whether the cell types that went in land in separate regions or
#    intermix, without having to go dig through earlier scripts' plots.
#
# Runs as a SLURM job array (see jobs/15_subclustering2.sh), one task per
# entry in subclustering_targets below -- some tissues get more than one
# task (e.g. brain has 5), so this is no longer a simple 3-task-per-tissue
# array. Each task independently rebuilds the whole tissue object before
# subsetting, same redundancy 13_subclustering1.R already accepted for
# full task parallelism -- for tissues with multiple targets, the
# round1/round2 label plots get regenerated (identically) by every task
# for that tissue, which is wasted but harmless compute.
#
# Design notes, since the request described this "very similarly to 13"
# without restating every parameter:
# - RunPCA() uses npcs = 50, matching 13_subclustering1.R (only the
#   Harmony dims = 1:20 was explicitly restated in the request).
# - Variable features/scaling are re-run from real raw counts
#   (data/06_obj_reassembly/bpcells, subset to each group's cells), same
#   as 13 and for the same reason -- a specific cell-type group's
#   within-population variation isn't what the broader object's variable
#   features were selected for, and FindVariableFeatures()'s default VST
#   fit needs real counts, not the normalized data this pipeline's
#   "counts" layer has held since 09_integration.
# - "the same resolutions from before" is read as 13_subclustering1.R's
#   list (this stage is explicitly modeled on 13), not 10_clustering.R's.
# - A graph_modularity.csv + resolution-vs-modularity plot are saved too,
#   matching 13's convention, even though not restated in this request --
#   free byproduct of the same clustering loop, useful for the same
#   review purpose.
# - Only ONE final DimPlot is saved (the selected/best resolution's
#   clusters), not one per resolution like 13 -- "save the new clusters
#   in a DimPlot like before" reads as matching 14_findmarkers2.R's
#   single cluster_dimplot.png, not 13's per-resolution DimPlot loop.
# - Multi-cell-type groups (e.g. c("Endothelial cell", "Pericyte",
#   "Smooth muscle cell")) are named for directories/files by joining the
#   cell types with "_" -- no semantic group name (e.g. "Vascular") was
#   given for these, so nothing was invented on the user's behalf.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(igraph)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# The same objective Leiden already optimized when it clustered. Reuses
# the SNN graph from FindNeighbors() below -- it's already a sparse
# Graph/dgCMatrix, so igraph builds from it without densifying. Copied
# from 10_clustering.R/13_subclustering1.R.
graph_modularity <- function(obj, clusters, graph_name = "RNA_snn"){
  snn <- as(obj[[graph_name]], "dgCMatrix")
  g <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::modularity(g, membership = as.integer(factor(clusters)))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Resolution chosen per tissue, matching 11_findmarkers.R/12_annotation1.R/
# 13_subclustering1.R exactly -- this is the resolution the round-1
# "cell_type1" annotation (results/12_annotation1/<tissue>_annotations.csv)
# was built from.
tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  resolution = c(1, 1, 1.2)
)

# Which cell_type2 groups to re-subcluster, one task each -- some entries
# combine multiple cell types into one group to double-check identity/
# separation between them.
subclustering_targets <- list(
  list(tissue = "brain", cell_types = "EN"),
  list(tissue = "brain", cell_types = "IN"),
  list(tissue = "brain", cell_types = "Myeloid"),
  list(tissue = "brain", cell_types = "Lymphocyte"),
  list(tissue = "brain", cell_types = "Vascular"),
  list(tissue = "muscle", cell_types = "Muscle fiber"),
  list(tissue = "muscle", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell")),
  list(tissue = "muscle", cell_types = c("Macrophage", "Monocyte", "Neutrophil")),
  list(tissue = "muscle", cell_types = c("T-cell", "NK cell")),
  list(tissue = "sc", cell_types = "Myeloid"),
  list(tissue = "sc", cell_types = "Lymphocyte"),
  list(tissue = "sc", cell_types = c("IN", "EN", "SN", "MN")),
  list(tissue = "sc", cell_types = c("Endothelial cell", "Pericyte", "Smooth muscle cell"))
)

# Figure out which target this task handles ----------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/15_subclustering2.sh), one task per ",
       "entry in subclustering_targets, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(subclustering_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(subclustering_targets), " subclustering targets -- ",
              "check the --array range in jobs/15_subclustering2.sh."))
}

target <- subclustering_targets[[task_id]]
tissue_file <- target$tissue
cell_types_target <- target$cell_types
group_label <- paste(cell_types_target, collapse = "_")

tissue_row <- tissues[tissues$file == tissue_file, ]
tissue_title <- tissue_row$title
resolution <- tissue_row$resolution
resolution_col <- paste0("res", resolution, "_clusters")

message2(paste0("Processing ", group_label, " (", tissue_title,
                "), task ", task_id, "/", length(subclustering_targets)))

results_dir <- "results/15_annotation2/"
dir.create(results_dir, showWarnings = F, recursive = T)

# Load the full tissue object used to create the subsets --------------------
# Same object 13_subclustering1.R built before subsetting: 09_integration's
# expression data + 10_clustering's cluster labels/UMAP, with the round-1
# annotation applied the same way (Idents() -> Pull_Cluster_Annotation() ->
# Rename_Clusters()) to get the "cell_type1" column each 13/14 subset was
# split on.

message2("Reading in full tissue object")

data_mat <- open_matrix_dir(paste0("data/09_integration/", tissue_file, "/bpcells_data"))
meta <- readRDS(paste0("data/10_clustering/", tissue_file, "/metadata.rds"))
umap <- readRDS(paste0("data/10_clustering/", tissue_file, "/harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony_umap"]] <- umap

message2("Applying round-1 cluster annotations")

Idents(obj) <- resolution_col

annots <- Pull_Cluster_Annotation(
  annotation = paste0("results/12_annotation1/", tissue_file, "_annotations.csv")
)

obj <- Rename_Clusters(obj,
                       new_idents = annots$new_cluster_idents,
                       new_ident_name = "cell_type1",
                       overwrite = T)

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type1",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(results_dir, tissue_file, "_round1_labels.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

# Fold in each subset's round-2 annotation -----------------------------------
# Cells not covered by any round-2 subcluster (e.g. a cell type not yet
# annotated) end up NA in cell_type2 -- fine here, since subset() below
# just naturally excludes them from any group (NA %in% x is FALSE).

message2("Folding in round-2 subcluster annotations")

params <- read.csv("jobs/13_params.txt", header = F,
                   col.names = c("cell_type", "tissue_file"))
cell_types <- params$cell_type[params$tissue_file == tissue_file]

new_labels <- list()

for (cell_type_target in cell_types){
  message2(paste0("  -- ", cell_type_target))

  sub_meta <- readRDS(paste0("data/13_subclustering1/", tissue_file, "/",
                             cell_type_target, "/metadata.rds"))

  modularity_df <- read.csv(paste0("results/13_subclustering1/", tissue_file, "/",
                                   cell_type_target, "/graph_modularity.csv"))
  best_res <- modularity_df$resolution[which.max(modularity_df$modularity)]
  sub_res_col <- paste0("res", best_res, "_clusters")

  # Matched directly by cluster value (not via Rename_Clusters()'s
  # positional alignment to levels(Idents(object))) -- no Seurat object is
  # needed for this join, and matching by value sidesteps any risk of a
  # silent off-by-one if cluster factor levels aren't sorted the way
  # Rename_Clusters() would assume.
  annot_csv <- read.csv(paste0("results/14_findmarkers2/", tissue_file, "/",
                               cell_type_target, "/annotations.csv"))

  new_labels[[cell_type_target]] <- annot_csv$cell_type[match(as.character(sub_meta[[sub_res_col]]),
                                          as.character(annot_csv$cluster))]
  names(new_labels[[cell_type_target]]) <- rownames(sub_meta)

}

new_labels <- list_c(new_labels) %>%
  as.data.frame()

colnames(new_labels) <- "cell_type2"

obj <- AddMetaData(obj,
                   new_labels)

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type2",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(results_dir, tissue_file, "_new_labels.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

# Subset to this task's cell_type2 group -------------------------------

message2(paste0("Subsetting to: ", group_label))

sub_obj <- subset(obj, subset = cell_type2 %in% cell_types_target)

group_data_dir <- paste0("data/15_subclustering2/", tissue_file, "/", group_label, "/")
dir.create(group_data_dir, showWarnings = F, recursive = T)

group_results_dir <- paste0("results/15_subclustering2/", tissue_file, "/", group_label, "/")
dir.create(group_results_dir, showWarnings = F, recursive = T)

# Rebuild from real raw counts for this group's cells ------------------------
# FindVariableFeatures()'s default "vst" method needs real counts -- see
# header note above.

message2("Loading raw counts for this group's cells")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, colnames(sub_obj)]

working_obj <- CreateSeuratObject(counts = raw_mat,
                                  meta.data = sub_obj@meta.data,
                                  assay = "RNA")

# Normalize, find variable features, scale, and run PCA ----------------------
# Small/rare groups may have very few cells in some samples (or none at
# all) -- if split()/IntegrateLayers() below fails or behaves oddly, check
# per-sample cell counts for this group first.

message2("Normalizing, finding variable features, scaling, and running PCA")

working_obj <- NormalizeData(working_obj)
working_obj <- FindVariableFeatures(working_obj)
working_obj <- ScaleData(working_obj)
working_obj <- RunPCA(working_obj, npcs = 50)

message2("Integrating samples using Harmony")

working_obj[["RNA"]] <- split(working_obj[["RNA"]], f = working_obj$orig.ident)

working_obj <- IntegrateLayers(working_obj,
                               method = "HarmonyIntegration",
                               orig.reduction = "pca",
                               new.reduction = "harmony",
                               dims = 1:20)

working_obj[["RNA"]] <- JoinLayers(working_obj[["RNA"]])

# Compute the NN graph and UMAP ------------------------------------------
# Matches 10_clustering.R/13_subclustering1.R's neighbor graph/UMAP block
# exactly.

message2("Computing neighbor graph and UMAP")

working_obj <- working_obj %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:20,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                return.neighbor = T) %>%
  FindNeighbors(reduction = "harmony",
                dims = 1:20,
                k.param = 15,
                nn.method = "annoy",
                annoy.metric = "euclidean",
                compute.SNN = T) %>%
  RunUMAP(umap.method = "uwot",
          nn.name = "RNA.nn",
          metric = "euclidean",
          min.dist = 0.5,
          n_neighbors = 15L,
          reduction.name = "harmony_umap",
          return.model = T)

# DimPlot of the cell_type2 labels that went into this group -----------------
# Useful on its own for single-cell-type targets, but especially for the
# multi-cell-type groups -- shows whether the input labels form distinct
# regions on the new embedding (a sanity check before trusting whatever
# the fresh clustering below finds) or fully intermix. Saved before
# clustering so it's available even if a cell type gets split up by the
# new clusters, without having to dig back through 12/13/14's plots.

message2("Making DimPlot of input cell_type2 labels")

p <- DimPlot_scCustom(working_obj,
                      group.by = "cell_type2",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(group_results_dir, "cell_type2_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

# Cluster at several resolutions, scoring each by graph modularity ----------
# Same resolution list as 13_subclustering1.R.

res_tests <- c(0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5)

modularity_vec <- numeric(length(res_tests))

for (i in seq_along(res_tests)){
  res <- res_tests[i]
  message2(paste0("Clustering at resolution = ", res))

  working_obj <- FindClusters(working_obj,
                              resolution = res,
                              algorithm = 4,
                              graph.name = "RNA_snn",
                              cluster.name = paste0("res", res, "_clusters"),
                              method = "igraph")

  modularity_vec[i] <- graph_modularity(working_obj,
                                        working_obj@meta.data[[paste0("res", res, "_clusters")]])
}

message2("Saving graph modularity table and plot")

modularity_df <- data.frame(resolution = res_tests, modularity = modularity_vec)

write.csv(modularity_df,
          file = paste0(group_results_dir, "graph_modularity.csv"),
          row.names = F)

p <- ggplot(modularity_df, aes(x = resolution, y = modularity)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(p,
       filename = paste0(group_results_dir, "graph_modularity.png"),
       units = "in", dpi = 300,
       height = 4, width = 6)

# Pick the resolution with the highest graph modularity and find markers ----

best_res <- res_tests[which.max(modularity_vec)]
best_res_col <- paste0("res", best_res, "_clusters")

message2(paste0("Best resolution = ", best_res))

Idents(working_obj) <- best_res_col

message2("Transposing matrix to row-major order for FindAllMarkers()")

# Same fix as 11_findmarkers.R/14_findmarkers2.R -- avoids repeated live
# re-transpositions during FindAllMarkers()'s per-gene test loop.
working_obj[["RNA"]]$data <- transpose_storage_order(working_obj[["RNA"]]$data)

message2("Finding cluster markers")

markers <- FindAllMarkers(working_obj)

write.csv(markers,
          file = paste0(group_results_dir, "markers.csv"))

message2("Making top 5 marker dot plot")

top5 <- markers %>%
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 5,
                      make_unique = T,
                      named_vector = F)

p <- dittoDotPlot(working_obj,
                  vars = top5,
                  group.by = best_res_col)
ggsave(p,
       filename = paste0(group_results_dir, "top5_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 20)

message2("Making cluster DimPlot")

p <- DimPlot_scCustom(working_obj,
                      group.by = best_res_col,
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(group_results_dir, "cluster_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

message2("Writing annotation template")

Create_Cluster_Annotation_File(file_path = group_results_dir,
                               file_name = "annotations")

# Save metadata, integrated embedding, normalized expression, and UMAP -----

message2("Saving metadata, count matrix, Harmony embedding, and UMAP")

bpcells_data_dir <- paste0(group_data_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = working_obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(working_obj@meta.data,
        file = paste0(group_data_dir, "metadata.rds"))

saveRDS(working_obj[["harmony"]],
        file = paste0(group_data_dir, "harmony.rds"))

saveRDS(working_obj[["harmony_umap"]],
        file = paste0(group_data_dir, "harmony_umap.rds"))
