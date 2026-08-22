# Makes diagnostic plots for manually annotating one tissue's clusters
# offline: a UMAP colored by cluster, a dot plot of each cluster's top 5
# marker genes (from 11_findmarkers.R), and a feature/violin plot pair for
# each of a tissue-specific set of canonical marker genes. The user
# annotates clusters by hand from these plots plus the marker CSV from 11,
# saving results into the annotation template 11 already writes;
# 13 (to be written) will load that annotation and label the clusters.
# Runs as a SLURM job array (see jobs/12_annotation1.sh), one task per
# tissue, since all 3 tissues are fully independent (same restructuring as
# 04_doubletfinder.R, 07_norm_pca.R, 08_sketch_pca.R, 10_clustering.R, and
# 11_findmarkers.R).

suppressMessages({
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(dittoSeq)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

# Saves a feature plot + violin plot pair for one gene. Reads plots_dir
# from the global environment (set once below, per-task) rather than
# taking it as an argument, matching how this was originally written.
combined_plot <- function(s, gene, res){
  p1 <- FeaturePlot_scCustom(s,
                             features = gene,
                             reduction = "umap",
                             raster = T,
                             raster.dpi = c(900, 900),
                             pt.size = 0.05)
  p2 <- VlnPlot_scCustom(s,
                         features = gene,
                         group.by = paste0("res", res, "_clusters"),
                         raster = T,
                         raster.dpi = 900) +
    NoLegend()

  p <- p1 + p2 +
    plot_layout(nrow = 2,
               heights = c(2, 1))

  ggsave(p,
         filename = paste0(plots_dir, gene, ".png"),
         units = "in", dpi = 300,
         height = 9, width = 7)
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

# Canonical marker genes checked for every tissue, tissue-specific since
# expected cell types differ (e.g. skeletal muscle fibers vs. CNS
# oligodendrocytes/microglia).
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

# Resolution chosen per tissue, matching 11_findmarkers.R exactly -- the
# markers.csv read below was computed at these resolutions, so the
# cluster labels used here have to match or the dot plot/violin plots
# would be grouping cells by the wrong clustering.
tissues <- data.frame(
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
  file = c("brain", "sc", "muscle"),
  resolution = c(1, 1, 1.2)
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/12_annotation1.sh), one task per tissue, ",
       "not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/12_annotation1.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]
resolution <- tissues$resolution[task_id]
genes <- gene_lists[[tissue_file]]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ") at resolution = ", resolution))

data_in_dir <- paste0("data/09_integration/", tissue_file, "/")
clusters_in_dir <- paste0("data/10_clustering/", tissue_file, "/")
markers_in_dir <- paste0("tab_data/11_findmarkers/", tissue_file, "/")

plots_dir <- paste0("plots/12_annotation1/", tissue_file, "/")
dir.create(plots_dir, showWarnings = F,
           recursive = T)

# Load the clustered object and its markers ----------------------------------

message2("Reading in expression data, cluster labels, UMAP, and markers")

data_mat <- open_matrix_dir(paste0(data_in_dir, "bpcells_data"))
meta <- readRDS(paste0(clusters_in_dir, "metadata.rds"))
umap <- readRDS(paste0(clusters_in_dir, "harmony_umap.rds"))
markers <- read.csv(paste0(markers_in_dir, "markers.csv"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["umap"]] <- umap

# Cluster UMAP ----------------------------------------------------------

message2("Making cluster DimPlot")

p <- DimPlot_scCustom(obj,
                      reduction = "umap",
                      group.by = paste0("res", resolution, "_clusters"))
ggsave(p,
       filename = paste0(plots_dir, "cluster_dimplot.png"),
       units = "in", dpi = 300,
       height = 6, width = 7)

# Top 5 marker genes per cluster -----------------------------------------

message2("Making top 5 marker dot plot")

top5 <- markers %>%
  filter(pct.1 > 0.3) %>%
  Extract_Top_Markers(num_features = 5,
                      make_unique = T,
                      named_vector = F)

p <- dittoDotPlot(obj,
                  vars = top5,
                  group.by = paste0("res", resolution, "_clusters"))
ggsave(p,
       filename = paste0(plots_dir, "top5_dotplot.png"),
       units = "in", dpi = 600,
       height = 8, width = 20)

# Feature/violin plots for canonical marker genes -------------------------

message2("Making feature/violin plots for canonical marker genes")

for (gene in genes){
  combined_plot(obj, gene, resolution)
}
