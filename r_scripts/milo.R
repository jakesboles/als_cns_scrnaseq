# Tests for differentially abundant neighborhoods (MiloR) between sALS
# vs. Control and C9orf72 vs. Control, run separately per tissue, for
# each subclustered population saved by 19_subclustering3.R. Runs as a
# SLURM job array (see jobs/milo.sh), one task per subdirectory of
# data/19_subclustering3/ -- discovered at runtime (not a hardcoded
# target list) so this stays correct if targets are ever added to/removed
# from 19_subclustering3.R without needing a matching edit here, per the
# user's "make it amenable to any object in that 19 directory" request.
# Each task internally loops over whichever tissues actually exist in its
# own target's metadata (1 for muscle_fiber, 2 for microglia/astrocyte/
# oligodendrocyte/neurons, 3 for myeloid) -- same "one array task, small
# internal loop over however many sub-units exist" shape as deseq2.R's
# per-tissue loop over cell types.
#
# Reworked from the user's pushed sample script (working code from the
# earlier als_multitissue_scfrp project). Design notes -- changes
# proposed and implemented, flagged for review rather than silently
# assumed:
# - The sample script's obj@graphs$RNA_snn doesn't exist in anything
#   19_subclustering3.R saves -- that script's scope deliberately stops
#   at integration + UMAP, no clustering/FindNeighbors() call at all (see
#   its own header). This script computes the SNN graph itself
#   (FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 15,
#   compute.SNN = T), matching every other script's established neighbor-
#   graph block) on each tissue-subsetted object, immediately before
#   building the Milo object from it -- not once on the whole
#   cross-tissue object, since a graph/neighborhoods built across tissues
#   and then subset afterward would still reflect cross-tissue structure,
#   defeating the point of running brain and spinal cord separately.
# - Doesn't need real raw counts at all -- MiloR's differential-abundance
#   test operates on the neighborhood-by-sample cell-count matrix and the
#   kNN graph/reduced-dims embedding, not on gene expression counts or
#   variable-feature selection, so unlike most other scripts in this
#   project there's no reason to reach back to
#   data/06_obj_reassembly/bpcells here. Loads 19_subclustering3.R's own
#   saved bpcells_data/metadata.rds/harmony.rds directly (same lightweight
#   reload pattern as 14_findmarkers2.R, which similarly only ever reads
#   the "data" layer).
# - reduced_dims = "INTEGRATED_PCA" (the old project's reduction name)
#   doesn't exist here -- changed to "harmony", this project's actual
#   integrated reduction, saved by 19_subclustering3.R.
# - d = 13 in the sample script looks like a leftover, dataset-specific
#   value (possibly from a JackStraw/elbow check on that old microglia
#   object specifically) with no equivalent justification here -- changed
#   to d = 20, matching the dims = 1:20 convention used everywhere else in
#   this project for Harmony/neighbor/UMAP steps.
# - Group/GroupsALS/GroupControl/GroupC9orf72 renamed to this project's
#   actual (lowercase) group column and the resulting
#   groupsALS/groupControl/groupC9orf72 model-matrix term names.
# - Added a min_cells check per (target, tissue) before attempting the
#   pipeline -- same "abundance filter, checked before the expensive/
#   fragile step" lesson already applied in deseq2.R/wgcna_single.R/
#   wgcna_consensus*.R, extended here since nothing in the sample script
#   guarded against a target/tissue combination too sparse for a stable
#   neighborhood-count design matrix (e.g. a donor with very few cells
#   contributing near-empty neighborhoods).
# - Added buildNhoodGraph(milo) before saving -- the sample script never
#   called this, but it's what populates nhoodGraph(milo), which
#   downstream plotting functions like plotNhoodGraphDA() need. Without
#   it the saved Milo object wouldn't actually be ready for the
#   "downstream plotting" the user asked for.
# - The sample script computed sals_results/c9_results but never saved
#   anything. Saves, per (target, tissue): the full Milo object (one
#   deliberate exception to "never a whole-object saveRDS()", matching
#   the same reasoning as hdWGCNA's @misc-object exception -- the Milo
#   object is the standard unit of persistence miloR's own downstream
#   plotting functions expect, and at subclustered-population,
#   single-tissue scale it isn't the whole-cohort expression matrix, so
#   this isn't the "duplicate a huge object per task" problem those
#   scripts were avoiding), the Seurat Graph object used to build it
#   (separately, as a plain dgCMatrix, in case anything needs the raw
#   graph without loading miloR/reconstructing a Milo object), and both
#   contrasts' DA results tables as CSVs.
# - Not implemented, flagged rather than silently added: annotateNhoods().
#   Each target's object is already restricted to one cell type by
#   construction (19_subclustering3.R), so annotating neighborhoods by
#   cell type wouldn't add information; nothing else obviously useful to
#   annotate by suggested itself. Easy to add later if a specific
#   downstream plot needs it.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(miloR)
  library(scater)
  library(igraph)
  library(scales)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

project_root <- "/projects/b1169/boles/als_cns_scrnaseq"
setwd(project_root)

set.seed(256)

# Figure out which 19_subclustering3.R target this task handles -----------
# Discovered from disk, not hardcoded -- see header note above.

subclustering_targets <- sort(list.dirs("data/19_subclustering3",
                                        recursive = F, full.names = F))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/milo.sh), one task per subdirectory of ",
       "data/19_subclustering3/, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(subclustering_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(subclustering_targets), " targets found in ",
              "data/19_subclustering3/ -- check the --array range in ",
              "jobs/milo.sh (it must match the number of subdirectories ",
              "there, which can change if 19_subclustering3.R's targets ",
              "change)."))
}

target_name <- subclustering_targets[task_id]

message2(paste0("Processing ", target_name, ", task ", task_id, "/",
                length(subclustering_targets)))

# Load this target's subclustered object -------------------------------
# No real raw counts needed -- see header note above.

message2("Reading in expression data, metadata, and Harmony embedding")

target_data_dir <- paste0("data/19_subclustering3/", target_name, "/")

data_mat <- open_matrix_dir(paste0(target_data_dir, "bpcells_data"))
meta <- readRDS(paste0(target_data_dir, "metadata.rds"))
harmony <- readRDS(paste0(target_data_dir, "harmony.rds"))

obj <- CreateSeuratObject(counts = data_mat, meta.data = meta, assay = "RNA")
obj[["RNA"]]$data <- data_mat
obj[["harmony"]] <- harmony

obj$group <- factor(obj$group, levels = c("Control", "sALS", "C9orf72"))

# Run the DA analysis separately for each tissue present in this target ---

tissues_present <- sort(unique(as.character(obj$tissue)))

message2(paste0("Tissues found for ", target_name, ": ",
                paste(tissues_present, collapse = ", ")))

min_cells <- 200 # change as needed
min_donors_per_group <- 3 # change as needed

for (tissue_title in tissues_present){

  message2(paste0("Processing ", target_name, " -- ", tissue_title))

  data_dir <- paste0(project_root, "/data/milo/", target_name, "/",
                     str_replace_all(tissue_title, " ", "_"), "/")
  dir.create(data_dir, showWarnings = F, recursive = T)

  results_dir <- paste0(project_root, "/results/milo/", target_name, "/",
                        str_replace_all(tissue_title, " ", "_"), "/")
  dir.create(results_dir, showWarnings = F, recursive = T)

  sub <- subset(obj, tissue == tissue_title)

  # Abundance check before attempting the pipeline -- see header note
  # above.
  if (ncol(sub) < min_cells){
    message2(paste0("Skipping ", target_name, " (", tissue_title, ") -- ",
                    "only ", ncol(sub), " cells, below min_cells = ",
                    min_cells, "."))
    next
  }
  donors_per_group <- sub@meta.data %>%
    distinct(id, group) %>%
    dplyr::count(group) %>%
    tidyr::complete(group, fill = list(n = 0))
  if (any(donors_per_group$n < min_donors_per_group)){
    message2(paste0("Skipping ", target_name, " (", tissue_title, ") -- ",
                    "fewer than ", min_donors_per_group, " donors in at ",
                    "least one group: ",
                    paste(donors_per_group$group, donors_per_group$n,
                          sep = " = ", collapse = ", ")))
    next
  }

  # Compute the neighbor graph on this tissue-subsetted object -------------
  # Not reused from a cross-tissue graph -- see header note above.

  message2("Computing neighbor graph")

  sub <- FindNeighbors(sub,
                       reduction = "harmony",
                       dims = 1:20,
                       k.param = 15,
                       nn.method = "annoy",
                       annoy.metric = "euclidean",
                       compute.SNN = T)

  # Build the Milo object -------------------------------------------------

  message2("Building Milo object")

  sce <- as.SingleCellExperiment(sub)

  milo <- Milo(sce)
  miloR::graph(milo) <- miloR::graph(buildFromAdjacency(sub@graphs$RNA_snn,
                                                        k = 15,
                                                        is.binary = F))

  milo <- makeNhoods(milo,
                     prop = 0.05, # change as needed
                     k = 15, # change as needed
                     d = 20,
                     refined = TRUE,
                     reduced_dims = "harmony")

  p <- plotNhoodSizeHist(milo)
  ggsave(p,
         filename = paste0(results_dir, "nhood_size_hist.png"),
         units = "in", dpi = 300,
         height = 5, width = 6)

  milo <- countCells(milo,
                     meta.data = as.data.frame(colData(milo)),
                     sample = "id")

  design <- data.frame(colData(milo))[, c("id", "group")]
  design <- distinct(design)
  design$group <- factor(design$group, levels = c("Control", "sALS", "C9orf72"))
  rownames(design) <- design$id

  milo <- calcNhoodDistance(milo,
                            d = 20,
                            reduced.dim = "harmony")

  # Populate the neighborhood graph for downstream plotting -----------------
  # The sample script never called this -- see header note above.

  milo <- buildNhoodGraph(milo)

  # Test both contrasts -----------------------------------------------------
  # ~ 0 + group so every level of group is its own column in the model
  # matrix, matching the sample script's approach.

  message2("Testing sALS vs. Control")

  sals_results <- testNhoods(milo,
                             design = ~ 0 + group,
                             design.df = design,
                             model.contrasts = "groupsALS - groupControl",
                             fdr.weighting = "graph-overlap",
                             norm.method = "TMM")

  write.csv(sals_results,
            file = paste0(results_dir, "sALS_vs_Control_nhood_results.csv"),
            row.names = F)

  message2("Testing C9orf72 vs. Control")

  c9_results <- testNhoods(milo,
                           design = ~ 0 + group,
                           design.df = design,
                           model.contrasts = "groupC9orf72 - groupControl",
                           fdr.weighting = "graph-overlap",
                           norm.method = "TMM")

  write.csv(c9_results,
            file = paste0(results_dir, "C9orf72_vs_Control_nhood_results.csv"),
            row.names = F)

  # Save the Milo object and the graph used to build it ---------------------
  # See header note above re: saving the whole Milo object.

  message2("Saving Milo object and neighbor graph")

  saveRDS(milo,
          file = paste0(data_dir, "milo_obj.rds"))

  saveRDS(sub@graphs$RNA_snn,
          file = paste0(data_dir, "snn_graph.rds"))

}
