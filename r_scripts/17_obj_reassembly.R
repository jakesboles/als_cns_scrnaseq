# Produces the final cleaned, integrated, fully annotated object for one
# tissue. Runs as a SLURM job array (see jobs/17_obj_reassembly.sh), one
# task per tissue, since all 3 tissues are fully independent (same
# restructuring as 04/07/08/09/10/11/12/13/14/15).
#
# Per task:
# 1. Rebuilds the full tissue object and folds in cell_type1 (round 1,
#    tab_data/12_annotation1/) and cell_type2 (round 2,
#    14_findmarkers2.R's per-cell-type annotations.csv) -- identical to
#    15_subclustering2.R's own folding logic (the "older 15_annotation2.R"
#    approach the user pointed back to), reused here rather than
#    reimplemented differently.
# 2. Folds in cell_type3 (round 3, 15_subclustering2.R's per-group
#    annotations.csv): for each subclustering_targets group belonging to
#    this tissue, matches each cell's cluster (at that group's
#    modularity-best resolution) to its annotations.csv row. A cell type
#    never selected as a subclustering target keeps its cell_type2 label
#    untouched (cell_type3 defaults to cell_type2 up front). Within a
#    processed group, any cluster left with a blank cell_type entry (the
#    user's way of saying "the existing label already looked right")
#    falls back to that cell's own cell_type2 value instead of being
#    overwritten with an empty string.
# 3. Removes cells where cell_type3 == "Remove".
# 4. Drops metadata clutter that's no longer meaningful now that
#    cell_type1/2/3 exist: all resX_clusters columns from 10_clustering.R's
#    resolution sweep, and DoubletFinder's pANN* columns. DF.unadj/DF.adj
#    are kept -- still meaningful QC info, confirmed with the user.
# 5. Rebuilds from real raw counts (data/06_obj_reassembly/bpcells,
#    subset to the retained cells) and re-runs NormalizeData/
#    FindVariableFeatures/ScaleData/RunPCA(npcs=100, matching
#    07_norm_pca.R's original full-tissue scale, not 13/15's npcs=50
#    cell-type-subset scale) for the same reason 13/15 did: this
#    pipeline's "counts" layer has held normalized data, not real counts,
#    since 09_integration, and FindVariableFeatures()'s default VST fit
#    needs real counts.
# 6. Re-integrates with Harmony (dims=1:20, matching every other
#    integration step) and computes one diagnostic UMAP -- matching
#    09_*_integration.R's scope, not 10_clustering.R's. No reclustering
#    here (confirmed with the user): cell_type3 is already the final
#    annotation, so there's no reason to run a fresh Leiden resolution
#    sweep on top of it.
# 7. Saves metadata, the normalized BPCells matrix, the Harmony embedding,
#    and the UMAP -- the final per-tissue object -- plus one DimPlot of
#    cell_type3 for a sanity check.

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
# 13_subclustering1.R/15_subclustering2.R exactly -- this is the
# resolution the round-1 "cell_type1" annotation was built from.
tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  resolution = c(1, 1, 1.2)
)

# Same list 15_subclustering2.R/16_subclustering2_inspection.R used --
# needed here to know which cell_type2 groups have a round-3 annotation
# to fold in for each tissue.
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

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/17_obj_reassembly.sh), one task per ",
       "tissue, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/17_obj_reassembly.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]
resolution <- tissues$resolution[task_id]
resolution_col <- paste0("res", resolution, "_clusters")

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

data_dir <- paste0("data/17_obj_reassembly/", tissue_file, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

plots_dir <- paste0("plots/17_obj_reassembly/", tissue_file, "/")
dir.create(plots_dir, showWarnings = F, recursive = T)

# Load the full tissue object and fold in round-1/round-2 annotations -------
# Identical to 15_subclustering2.R's own folding logic.

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
                      group.by = "cell_type2",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_cell_type1_dimplot_pre-integration.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

message2("Folding in round-2 subcluster annotations")

params <- read.csv("jobs/13_params.txt", header = F,
                   col.names = c("cell_type", "tissue_file"))
cell_types <- params$cell_type[params$tissue_file == tissue_file]

new_labels <- list()

for (cell_type_target in cell_types){
  sub_meta <- readRDS(paste0("data/13_subclustering1/", tissue_file, "/",
                             cell_type_target, "/metadata.rds"))

  modularity_df <- read.csv(paste0("tab_data/13_subclustering1/", tissue_file, "/",
                                   cell_type_target, "/graph_modularity.csv"))
  best_res <- modularity_df$resolution[which.max(modularity_df$modularity)]
  sub_res_col <- paste0("res", best_res, "_clusters")

  annot_csv <- read.csv(paste0("tab_data/14_findmarkers2/", tissue_file, "/",
                               cell_type_target, "/annotations.csv"))

  new_labels[[cell_type_target]] <- annot_csv$cell_type[match(as.character(sub_meta[[sub_res_col]]),
                                          as.character(annot_csv$cluster))]
  names(new_labels[[cell_type_target]]) <- rownames(sub_meta)
}

new_labels <- list_c(new_labels) %>%
  as.data.frame()
colnames(new_labels) <- "cell_type2"

obj <- AddMetaData(obj, new_labels)

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type2",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_cell_type2_dimplot_pre-integration.png"),
       units = "in", dpi = 300,
       height = 6, width = 8)

# obj <- subset(obj,
#               cell_type2 != "Remove")

# Fold in round-3 subcluster annotations -------------------------------
# Defaults every cell to its round-2 label -- covers cell types that were
# never selected as a subclustering target at all. Within a processed
# group, a blank annotations.csv row (cluster left unfilled because the
# existing label already looked right) falls back to that cell's own
# round-2 label instead of an empty string.

message2("Folding in round-3 subcluster annotations")

obj$cell_type3 <- obj$cell_type2

targets_this_tissue <- Filter(function(t) t$tissue == tissue_file, subclustering_targets)

for (target in targets_this_tissue){
  group_label <- paste(target$cell_types, collapse = "_")
  message2(paste0("  -- ", group_label))

  group_dir <- paste0("data/15_subclustering2/", tissue_file, "/", group_label, "/")
  group_tab_dir <- paste0("tab_data/15_subclustering2/", tissue_file, "/", group_label, "/")

  sub_meta <- readRDS(paste0(group_dir, "metadata.rds"))

  modularity_df <- read.csv(paste0(group_tab_dir, "graph_modularity.csv"))
  best_res <- modularity_df$resolution[which.max(modularity_df$modularity)]
  res_col <- paste0("res", best_res, "_clusters")

  annot_csv <- read.csv(paste0(group_tab_dir, "annotations.csv"))

  new_round3 <- annot_csv$cell_type[match(as.character(sub_meta[[res_col]]),
                                          as.character(annot_csv$cluster))]
  names(new_round3) <- rownames(sub_meta)

  blank <- is.na(new_round3) | trimws(new_round3) == ""
  new_round3[blank] <- obj$cell_type2[names(new_round3)[blank]]

  obj$cell_type3[names(new_round3)] <- new_round3
}

# Accidentally left a space after a remove somewhere, remove it here
obj$cell_type3 <- str_replace_all(obj$cell_type3, "Remove ", "Remove")

p <- DimPlot_scCustom(obj,
                      group.by = "cell_type3",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_cell_type3_dimplot_pre-integration.png"),
       units = "in", dpi = 300,
       height = 6, width = 10)

# Remove cells flagged for exclusion -----------------------------------

message2("Removing cells labeled 'Remove'")

n_before <- ncol(obj)
obj <- subset(obj, subset = cell_type3 != "Remove")
message2(paste0("Removed ", n_before - ncol(obj), " of ", n_before, " cells"))

# Clean up metadata -----------------------------------------------------
# Drops the resolution sweep from 10_clustering.R (superseded by
# cell_type1/2/3) and DoubletFinder's pANN* columns. DF.unadj/DF.adj are
# kept.

meta_clean <- obj@meta.data
cell_names <- rownames(meta_clean)

meta_clean <- meta_clean %>%
  select(-matches("pANN"),
         -matches("^res[0-9.]+_clusters$"),
         -matches("^seurat_clusters$"),
         -matches("^RNA_snn_res\\."))

rownames(meta_clean) <- cell_names
obj@meta.data <- meta_clean

# Rebuild from real raw counts for the cleaned cell set ----------------------
# FindVariableFeatures()'s default "vst" method needs real counts -- see
# header note above.

message2("Loading raw counts for the cleaned cell set")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, colnames(obj)]

working_obj <- CreateSeuratObject(counts = raw_mat,
                                  meta.data = obj@meta.data,
                                  assay = "RNA")

# Normalize, find variable features, scale, and run PCA ----------------------

message2("Normalizing, finding variable features, scaling, and running PCA")

working_obj <- NormalizeData(working_obj)
working_obj <- FindVariableFeatures(working_obj)
working_obj <- ScaleData(working_obj)
working_obj <- RunPCA(working_obj, npcs = 100)

message2("Integrating samples using Harmony")

working_obj[["RNA"]] <- split(working_obj[["RNA"]], f = working_obj$orig.ident)

working_obj <- IntegrateLayers(working_obj,
                               method = "HarmonyIntegration",
                               orig.reduction = "pca",
                               new.reduction = "harmony",
                               dims = 1:20)

working_obj[["RNA"]] <- JoinLayers(working_obj[["RNA"]])

message2("Computing UMAP for diagnostics")

# Matches 09_*_integration.R's scope -- a diagnostic UMAP only, no
# clustering. RunUMAP() searches directly off the Harmony reduction, no
# separate FindNeighbors() call needed first.
working_obj <- RunUMAP(working_obj,
                       dims = 1:20,
                       reduction = "harmony",
                       reduction.name = "harmony_umap",
                       reduction.key = "harmonyumap_")

message2("Making cell_type3 DimPlot")

p <- DimPlot_scCustom(working_obj,
                      group.by = "cell_type3",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, tissue_file, "_cell_type3_dimplot_post-integration.png"),
       units = "in", dpi = 300,
       height = 6, width = 10)

# Save the cleaned, integrated, fully annotated object -----------------

message2("Saving metadata, count matrix, Harmony embedding, and UMAP")

bpcells_data_dir <- paste0(data_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = working_obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(working_obj@meta.data,
        file = paste0(data_dir, "metadata.rds"))

saveRDS(working_obj[["harmony"]],
        file = paste0(data_dir, "harmony.rds"))

saveRDS(working_obj[["harmony_umap"]],
        file = paste0(data_dir, "harmony_umap.rds"))
