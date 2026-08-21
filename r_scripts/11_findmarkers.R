# Finds cluster markers (FindAllMarkers()) for one tissue at its chosen
# clustering resolution, and writes out a blank cluster-annotation
# template alongside them. Runs as a SLURM job array (see
# jobs/11_findmarkers.sh), one task per tissue, since all 3 tissues are
# fully independent (same restructuring as 04_doubletfinder.R,
# 07_norm_pca.R, 08_sketch_pca.R, and 10_clustering.R).
#
# Only needs the expression matrix (unchanged since 09_*_integration.R --
# 10_clustering.R doesn't touch it) and cluster labels (from
# 10_clustering.R's metadata.rds) -- nothing else is required for
# FindAllMarkers(), which tests on the normalized "data" layer and doesn't
# need counts or scale.data.

# Load libraries
suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Clustering resolution chosen per tissue (via graph modularity -- see
# tab_data/10_clustering/<tissue_file>/graph_modularity.csv), matching the
# choices made in the earlier version of this project.
tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle"),
  resolution = c(1, 1, 1.2)
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/11_findmarkers.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/11_findmarkers.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]
resolution <- tissues$resolution[task_id]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ") at resolution = ", resolution))

data_in_dir <- paste0("data/09_integration/", tissue_file, "/")
clusters_in_dir <- paste0("data/10_clustering/", tissue_file, "/")

tab_data_out_dir <- paste0("tab_data/11_findmarkers/", tissue_file, "/")
dir.create(tab_data_out_dir, showWarnings = F,
           recursive = T)

# Load expression data and cluster labels ------------------------------------

message2("Reading in expression data and cluster labels")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(clusters_in_dir, "metadata.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat

Idents(obj) <- paste0("res", resolution, "_clusters")

message2("Finding cluster markers")

markers <- FindAllMarkers(obj)

message2("Saving markers and annotation template")

write.csv(markers,
          file = paste0(tab_data_out_dir, "markers.csv"))

Create_Cluster_Annotation_File(file_path = tab_data_out_dir,
                               file_name = paste0("res", resolution, "_annotations"))
