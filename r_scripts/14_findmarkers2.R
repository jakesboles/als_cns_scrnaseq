# Second-round marker finding for each of 13_subclustering1.R's cell-type
# subclusterings, to help annotate the finer-grained clusters it produced.
# Runs as a SLURM job array (see jobs/14_findmarkers2.sh), one task per
# (tissue, cell type) combination listed in jobs/13_params.txt, same as
# 13_subclustering1.R. Only markers/plots/an annotation template are
# produced -- no expression or metadata output, since nothing downstream
# needs a modified object from this script (15 will either do a third
# annotation round or assemble the finalized objects, working from 13's
# saved data plus the annotation.csv this script writes).
#
# Design notes (confirmed with the user before writing this):
# - "old_labels" in the original request doesn't match any actual
#   metadata column -- 13_subclustering1.R's saved metadata has
#   "cell_type" (constant within each file, since each file is already
#   subset to one cell type) and "round1_clusters" (the original numeric
#   round-1 cluster labels). Using round1_clusters, which is the one that
#   can actually show something on a DimPlot here.
# - Output directories use tab_data/14_findmarkers2/... and
#   plots/14_findmarkers2/..., matching every other script's convention
#   (11/12/13), not the reversed 14_findmarkers2/tab_data/... path in the
#   original request -- that would have created a new top-level folder
#   outside the existing .gitignore coverage (plots/, tab_data/, data/,
#   logs/ only).

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(dittoSeq)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Canonical marker genes checked for every tissue -- same lists as
# 12_annotation1.R.
brain_genes <- c("C1QA", "ITGAM", "P2RY12", "GFAP", "AQP4", "PDGFRA", "DCN",
                 "CEMIP", "RBFOX3", "OLIG1", "MOBP", "PLP1", "CLDN5", "PODXL",
                 "RGS5", "ACTA2", "NKG7", "GAD1", "SNAP25", "CNP", "OPALIN")
muscle_genes <- c("ITGAM", "ITGAX", "MS4A1", "CD3E", "NKG7", "S100A8", "APOE",
                  "APOD", "CLDN5", "PODXL", "RGS5", "ACTA2", "SNAP25", "RBFOX3",
                  "TNNC2", "TNNT1", "MYL2", "LYZ", "CD300E", "TPSB2")
sc_genes <- c("C1QA", "ITGAM", "P2RY12", "GFAP", "AQP4", "PDGFRA", "DCN",
              "CEMIP", "RBFOX3", "OLIG1", "MOBP", "PLP1", "CLDN5", "PODXL",
              "RGS5", "ACTA2", "NKG7", "GAD1", "SNAP25", "CFAP157", "TRAC",
              "CNP", "MPZ", "OPALIN")
gene_lists <- list(brain = brain_genes, sc = sc_genes, muscle = muscle_genes)

tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")
)

# Figure out which (tissue, cell type) combination this task handles ------
# Same jobs/13_params.txt used by 13_subclustering1.R, so this iterates
# over exactly the same set of objects.

params <- read.csv("jobs/13_params.txt", header = F,
                   col.names = c("cell_type", "tissue_file"))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/14_findmarkers2.sh), one task per ",
       "tissue/cell-type combination in jobs/13_params.txt, not as a ",
       "standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(params)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(params), " tissue/cell-type combinations in ",
              "jobs/13_params.txt -- check the --array range in ",
              "jobs/14_findmarkers2.sh."))
}

cell_type_target <- params$cell_type[task_id]
tissue_file <- params$tissue_file[task_id]
genes <- gene_lists[[tissue_file]]

tissue_row <- tissues[tissues$file == tissue_file, ]
if (nrow(tissue_row) != 1){
  stop(paste0("Tissue '", tissue_file, "' from jobs/13_params.txt (row ",
              task_id, ") doesn't match any entry in the tissues table -- ",
              "check for a typo in jobs/13_params.txt."))
}
tissue_title <- tissue_row$title

message2(paste0("Processing ", cell_type_target, " (", tissue_title,
                "), task ", task_id, "/", nrow(params)))

data_in_dir <- paste0("data/13_subclustering1/", tissue_file, "/",
                      cell_type_target, "/")
modularity_in_dir <- paste0("tab_data/13_subclustering1/", tissue_file, "/",
                            cell_type_target, "/")

tab_data_dir <- paste0("tab_data/14_findmarkers2/", tissue_file, "/",
                       cell_type_target, "/")
dir.create(tab_data_dir, showWarnings = F, recursive = T)

plots_dir <- paste0("plots/14_findmarkers2/", tissue_file, "/",
                    cell_type_target, "/")
dir.create(plots_dir, showWarnings = F, recursive = T)

# Load the object from 13 and prepare it for FindAllMarkers() ---------------

message2("Reading in expression data and metadata")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(data_in_dir, "metadata.rds"))
umap <- readRDS(paste0(data_in_dir, "harmony_umap.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony_umap"]] <- umap

message2("Transposing matrix to row-major order for FindAllMarkers()")

# Same fix as 11_findmarkers.R -- data_mat is stored column-major (cell-
# major), but FindAllMarkers() iterates per-gene (row-major). Transposing
# once up front avoids repeated live re-transpositions during the test
# loop. Writes to a temp directory local to this run, not back to
# data/13_subclustering1/<tissue>/<cell_type>/bpcells_data.
obj[["RNA"]]$data <- transpose_storage_order(obj[["RNA"]]$data)

# Pick the resolution with the highest graph modularity ----------------------

message2("Selecting resolution by graph modularity")

modularity_df <- read.csv(paste0(modularity_in_dir, "graph_modularity.csv"))
best_res <- modularity_df$resolution[which.max(modularity_df$modularity)]
resolution_col <- paste0("res", best_res, "_clusters")

message2(paste0("Best resolution = ", best_res))

# Find markers at that resolution --------------------------------------------

message2("Finding cluster markers")

Idents(obj) <- resolution_col

markers <- FindAllMarkers(obj)

write.csv(markers,
          file = paste0(tab_data_dir, "markers.csv"))

# Feature plots for canonical marker genes -----------------------------------

message2("Making feature plots for canonical marker genes")

for (gene in genes){
  p <- FeaturePlot_scCustom(obj,
                            features = gene,
                            reduction = "harmony_umap",
                            # raster = T,
                            # raster.dpi = c(900, 900),
                            pt.size = 0.05)
  ggsave(p,
         filename = paste0(plots_dir, gene, ".png"),
         units = "in", dpi = 300,
         height = 5, width = 6)
}

# Diagnostic DimPlots ---------------------------------------------------

message2("Making diagnostic DimPlots")

p <- DimPlot_scCustom(obj,
                      group.by = resolution_col,
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, "cluster_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

p <- DimPlot_scCustom(obj,
                      group.by = "round1_clusters",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, "round1_clusters_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

p <- DimPlot_scCustom(obj,
                      group.by = "DF.unadj",
                      reduction = "harmony_umap")
ggsave(p,
       filename = paste0(plots_dir, "df_unadj_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

p <- dittoBarPlot(obj,
                  var = "DF.unadj",
                  group.by = resolution_col)
ggsave(p,
       filename = paste0(plots_dir, "df_unadj_barplot.png"),
       units = "in", dpi = 300,
       height = 5, width = 7)

# Annotation template for this round -----------------------------------------

message2("Writing annotation template")

Create_Cluster_Annotation_File(file_path = tab_data_dir,
                               file_name = "annotations")
