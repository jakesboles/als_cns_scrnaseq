# Applies the first-pass cluster annotations from 12_annotation1.R to one
# tissue, then subsets to one annotated cell type and re-clusters it from
# scratch (fresh PCA/Harmony/neighbors/UMAP/Leiden), scoring each tested
# resolution by graph modularity. Runs as a SLURM job array (see
# jobs/13_subclustering1.sh), one task per (tissue, cell type) combination
# listed in jobs/13_params.txt, so all combinations run in parallel
# instead of 3 tissues x their cell types in nested for loops.
#
# Design notes (confirmed with the user before writing this):
# - "as in 11*.R (lines 160-179)" and "graph modularity as in 11*.R" in
#   the original request both matched content that only exists in
#   10_clustering.R (the FindNeighbors()/RunUMAP() block and the
#   graph_modularity() helper) -- 11_findmarkers.R has neither. Both are
#   copied from 10_clustering.R here.
# - Re-clustering a cell type needs its own FindVariableFeatures()/
#   ScaleData(), not the full tissue's (those are dominated by genes that
#   distinguish cell types from each other, not within-cell-type
#   variation). FindVariableFeatures()'s default "vst" selection method
#   needs real counts to fit correctly, but the object loaded from
#   09_integration/10_clustering has normalized data in its "counts"
#   layer (same pattern as 09/11/12) -- so real raw counts are pulled
#   from data/06_obj_reassembly/bpcells for this cell type's cells
#   specifically, and NormalizeData()/FindVariableFeatures()/ScaleData()/
#   RunPCA() are all re-run fresh on those, mirroring 07_norm_pca.R's
#   original recipe rather than mixing data sources.
# - RunPCA() uses npcs = 50 (not the full object's 100) since a single
#   cell type is a smaller, more homogeneous population; dims = 1:20 is
#   still used for Harmony/neighbors/UMAP, matching the full object.
# - Only metadata columns containing "pANN" are dropped, per the literal
#   request -- other DoubletFinder columns (e.g. DF.classifications) are
#   left in place.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
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
# from 10_clustering.R (see header note above).
graph_modularity <- function(obj, clusters, graph_name = "RNA_snn"){
  snn <- as(obj[[graph_name]], "dgCMatrix")
  g <- igraph::graph_from_adjacency_matrix(snn, mode = "undirected",
                                           weighted = TRUE, diag = FALSE)
  igraph::modularity(g, membership = as.integer(factor(clusters)))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Resolution chosen per tissue, matching 11_findmarkers.R/12_annotation1.R
# exactly -- the annotation CSVs read below were built from cluster
# labels at these resolutions.
tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  resolution = c(1, 1, 1.2)
)

# Figure out which (tissue, cell type) combination this task handles ------

params <- read.csv("jobs/13_params.txt", header = F,
                   col.names = c("cell_type", "tissue_file"))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/13_subclustering1.sh), one task per ",
       "tissue/cell-type combination in jobs/13_params.txt, not as a ",
       "standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(params)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(params), " tissue/cell-type combinations in ",
              "jobs/13_params.txt -- check the --array range in ",
              "jobs/13_subclustering1.sh."))
}

cell_type_target <- params$cell_type[task_id]
tissue_file <- params$tissue_file[task_id]

tissue_row <- tissues[tissues$file == tissue_file, ]
if (nrow(tissue_row) != 1){
  stop(paste0("Tissue '", tissue_file, "' from jobs/13_params.txt (row ",
              task_id, ") doesn't match any entry in the tissues table -- ",
              "check for a typo in jobs/13_params.txt."))
}
tissue_title <- tissue_row$title
resolution <- tissue_row$resolution
resolution_col <- paste0("res", resolution, "_clusters")

message2(paste0("Processing ", cell_type_target, " (", tissue_title,
                "), task ", task_id, "/", nrow(params)))

# Directories --------------------------------------------------------------

data_dir <- paste0("data/13_subclustering1/", tissue_file, "/",
                   cell_type_target, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

plots_dir <- paste0("plots/13_subclustering1/", tissue_file, "/",
                    cell_type_target, "/")
dir.create(plots_dir, showWarnings = F, recursive = T)

tab_data_dir <- paste0("tab_data/13_subclustering1/", tissue_file, "/",
                       cell_type_target, "/")
dir.create(tab_data_dir, showWarnings = F, recursive = T)

# Load the object used in 12_annotation1.R -----------------------------------

message2("Reading in expression data and cluster labels")

data_mat <- open_matrix_dir(paste0("data/09_integration/", tissue_file, "/bpcells_data"))
meta <- readRDS(paste0("data/10_clustering/", tissue_file, "/metadata.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat

# Apply first-pass annotations ----------------------------------------------

message2("Applying cluster annotations")

Idents(obj) <- resolution_col

annots <- Pull_Cluster_Annotation(
  annotation = paste0("tab_data/12_annotation1/", tissue_file, "_annotations.csv")
)

obj <- Rename_Clusters(obj,
                       new_idents = annots$new_cluster_idents,
                       new_ident_name = "cell_type",
                       overwrite = T)

# Clean up metadata: drop DoubletFinder pANN columns and every clustering
# column from 10_clustering.R except the one actually used, which gets
# renamed instead of dropped.
meta_clean <- obj@meta.data
cell_names <- rownames(meta_clean)

meta_clean <- meta_clean %>%
  rename(round1_clusters = all_of(resolution_col)) %>%
  select(-matches("pANN"),
         -matches("^res[0-9.]+_clusters$"),
         -matches("^seurat_clusters$"),
         -matches("^RNA_snn_res\\."))

rownames(meta_clean) <- cell_names
obj@meta.data <- meta_clean

# Subset to this task's cell type --------------------------------------------

message2(paste0("Subsetting to cell type: ", cell_type_target))

sub_obj <- subset(obj, subset = cell_type == cell_type_target)

# Rebuild from real raw counts for this cell type's cells --------------------
# FindVariableFeatures()'s default "vst" method needs real counts, not the
# normalized data carried in this object's "counts" layer -- see header
# note above.

message2("Loading raw counts for this cell type's cells")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, colnames(sub_obj)]

working_obj <- CreateSeuratObject(counts = raw_mat,
                                  meta.data = sub_obj@meta.data,
                                  assay = "RNA")

# Normalize, find variable features, scale, and run PCA ----------------------
# Small/rare cell types may have very few cells in some samples (or none
# at all) -- if split()/IntegrateLayers() below fails or behaves oddly,
# check per-sample cell counts for this cell type first.

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
# Matches 10_clustering.R's neighbor graph/UMAP block exactly (see header
# note above).

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

# Cluster at several resolutions, scoring each by graph modularity ----------

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

  p <- DimPlot_scCustom(working_obj,
                        group.by = paste0("res", res, "_clusters"),
                        label = F,
                        reduction = "harmony_umap")
  ggsave(p,
         filename = paste0(plots_dir, "res", res, ".png"),
         units = "in", dpi = 600,
         height = 5, width = 6)

  modularity_vec[i] <- graph_modularity(working_obj,
                                        working_obj@meta.data[[paste0("res", res, "_clusters")]])
}

message2("Saving graph modularity table and plot")

modularity_df <- data.frame(resolution = res_tests, modularity = modularity_vec)

write.csv(modularity_df,
          file = paste0(tab_data_dir, "graph_modularity.csv"),
          row.names = F)

p <- ggplot(modularity_df, aes(x = resolution, y = modularity)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
ggsave(p,
       filename = paste0(plots_dir, "graph_modularity.png"),
       units = "in", dpi = 300,
       height = 4, width = 6)

# Save metadata, count matrix, and UMAP as BPCells-ready data ---------------

message2("Saving metadata, count matrix, and UMAP")

bpcells_data_dir <- paste0(data_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = working_obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(working_obj@meta.data,
        file = paste0(data_dir, "metadata.rds"))

saveRDS(working_obj[["harmony_umap"]],
        file = paste0(data_dir, "harmony_umap.rds"))
