# Folds the round-2 (per-cell-type subcluster) annotations from
# 14_findmarkers2.R back onto the full tissue object, for a first look at
# how the refined labels lay out over the whole tissue's UMAP before
# building anything further on top of them. Runs as a SLURM job array
# (see jobs/15_annotation2.sh), one task per tissue, since all 3 tissues
# are fully independent (same restructuring as 04/07/08/10/11/12/13/14).
# Doesn't save metadata or a BPCells matrix yet -- just the diagnostic
# DimPlot, per the user's request; a later script will build on this one
# to save the finalized, fully-annotated object.
#
# Note: the user's message referenced "tab_data/14_findmarkers/..."; the
# actual script/output directory is 14_findmarkers2 (this pipeline's
# script numbering has never had a plain "14_findmarkers"), so that's
# what's read here.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Resolution chosen per tissue, matching 11_findmarkers.R/12_annotation1.R/
# 13_subclustering1.R exactly -- this is the resolution the round-1
# "cell_type" annotation (tab_data/12_annotation1/<tissue>_annotations.csv)
# was built from.
tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  resolution = c(1, 1, 1.2)
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/15_annotation2.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/15_annotation2.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]
resolution <- tissues$resolution[task_id]
resolution_col <- paste0("res", resolution, "_clusters")

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

plots_dir <- "plots/15_annotation2/"
dir.create(plots_dir, showWarnings = F, recursive = T)

# Load the full tissue object used to create the subsets --------------------
# Same object 13_subclustering1.R built before subsetting: 09_integration's
# expression data + 10_clustering's cluster labels/UMAP, with the round-1
# annotation applied the same way (Idents() -> Pull_Cluster_Annotation() ->
# Rename_Clusters()) to get the "cell_type" column each subset was split on.

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
  annotation = paste0("tab_data/12_annotation1/", tissue_file, "_annotations.csv")
)

obj <- Rename_Clusters(obj,
                       new_idents = annots$new_cluster_idents,
                       new_ident_name = "cell_type1",
                       overwrite = T)

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type1",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_round1_labels.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

# Fold in each subset's round-2 annotation -----------------------------------
# Defaults every cell to its round-1 label, then overwrites cells belonging
# to each cell type with that subset's refined round-2 label -- any cell
# type whose round-2 annotation isn't done yet just keeps showing its
# round-1 label instead of going missing from the plot.

message2("Folding in round-2 subcluster annotations")

# LEFT OFF HERE
cells <- list.dirs(paste0("data/13_subclustering1/", tissue_file),
                   full.names = F,
                   recursive = F)

meta_list <- list()

annotation_list <- list()

for (i in seq_along(cell_types)) {
  meta_list[[i]] <- readRDS(paste0("data/13_subclustering1/", tissue_file, "/", cells[i], "/metadata.rds"))
  
  annotation_list[[i]] <- Pull_Cluster_Annotation(
    annotation = paste0("tab_data/14_findmarkers2/", tissue_file, "/", cells[i], "/annotations.csv")
  )
  
}
# TO HERE

obj$cell_type_round2 <- obj$cell_type

params <- read.csv("jobs/13_params.txt", header = F,
                   col.names = c("cell_type", "tissue_file"))
cell_types <- params$cell_type[params$tissue_file == tissue_file]

for (cell_type_target in cell_types){
  message2(paste0("  -- ", cell_type_target))

  sub_meta <- readRDS(paste0("data/13_subclustering1/", tissue_file, "/",
                             cell_type_target, "/metadata.rds"))

  modularity_df <- read.csv(paste0("tab_data/13_subclustering1/", tissue_file, "/",
                                   cell_type_target, "/graph_modularity.csv"))
  best_res <- modularity_df$resolution[which.max(modularity_df$modularity)]
  sub_res_col <- paste0("res", best_res, "_clusters")

  # Matched directly by cluster value (not via Rename_Clusters()'s
  # positional alignment to levels(Idents(object))) -- no Seurat object is
  # needed for this join, and matching by value sidesteps any risk of a
  # silent off-by-one if cluster factor levels aren't sorted the way
  # Rename_Clusters() would assume.
  annot_csv <- read.csv(paste0("tab_data/14_findmarkers2/", tissue_file, "/",
                               cell_type_target, "/annotations.csv"))

  new_labels <- annot_csv$cell_type[match(as.character(sub_meta[[sub_res_col]]),
                                          as.character(annot_csv$cluster))]
  names(new_labels) <- rownames(sub_meta)

  obj$cell_type_round2[names(new_labels)] <- new_labels
}

# Diagnostic DimPlot ------------------------------------------------------

message2("Making DimPlot of round-2 labels")

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type_round2",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_new_labels.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)
