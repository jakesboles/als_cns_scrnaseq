# Tests for differentially abundant neighborhoods (MiloR) between sALS
# vs. Control and C9orf72 vs. Control, for each subclustered population
# saved by 19_subclustering3.R. The neighbor graph and neighborhoods are
# built once per target on the full cross-tissue object -- neighborhood
# structure is shared across tissues by design, so that a given
# neighborhood's abundance change can be compared as shared vs.
# tissue-specific -- and the two contrasts are then tested separately per
# tissue by restricting testNhoods()'s design.df to that tissue's
# samples, reusing the same shared neighborhoods rather than rebuilding a
# tissue-specific graph. Runs as a SLURM job array (see jobs/milo.sh),
# one task per subdirectory of data/19_subclustering3/ -- discovered at
# runtime (not a hardcoded target list) so this stays correct if targets
# are ever added to/removed from 19_subclustering3.R, per the user's
# "make it amenable to any object in that 19 directory" request.
#
# Reworked from the user's pushed sample script (working code from the
# earlier als_multitissue_scfrp project). Design notes -- changes
# proposed and implemented, flagged for review rather than silently
# assumed:
# - The sample script's obj@graphs$RNA_snn doesn't exist in anything
#   19_subclustering3.R saves -- that script's scope deliberately stops
#   at integration + UMAP, no clustering/FindNeighbors() call at all (see
#   its own header). This script computes the SNN graph itself, once per
#   target on the full cross-tissue object (FindNeighbors(reduction =
#   "harmony", dims = 1:15, k.param = 15, nn.method = "annoy",
#   annoy.metric = "euclidean", compute.SNN = T), per the user -- dims =
#   1:15 here specifically, not this project's usual dims = 1:20).
# - countCells()/testNhoods() key samples by orig.ident, not id -- id is
#   donor-only (deliberately, elsewhere in this pipeline, to give a clean
#   per-subject label), but since neighborhoods now span multiple tissues
#   in one object, a donor contributing samples to both brain and spinal
#   cord would collide under id and have those two samples' cells merged
#   into one counted column. orig.ident (donor + tissue) keeps them
#   separate, which is what the per-tissue design.df subsetting below
#   depends on.
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
#   integrated reduction, saved by 19_subclustering3.R. d = 13 in the
#   sample script is changed to d = 15, matching the dims = 1:15 used for
#   FindNeighbors above (kept consistent with each other, per the user).
# - Group/GroupsALS/GroupControl/GroupC9orf72 renamed to this project's
#   actual (lowercase) group column and the resulting
#   groupsALS/groupControl/groupC9orf72 model-matrix term names.
# - Per the user, both tissues' output goes in the same subfolder (no
#   per-tissue subdirectory) -- results/milo/<target>/ holds one shared
#   nhood_size_hist.png (neighborhoods are shared across tissues now, so
#   there's only one histogram to make) plus tissue-suffixed CSVs for
#   each contrast (e.g. sALS_vs_Control_Motor_cortex_nhood_results.csv),
#   and data/milo/<target>/ holds one shared milo_obj.rds/snn_graph.rds,
#   for exactly the same reason -- there's one graph/Milo object per
#   target now, not one per tissue.
# - testNhoods() accepting a design.df that's a strict subset of
#   nhoodCounts()'s samples (confirmed against miloR's own source/docs,
#   not guessed) is what makes per-tissue testing against shared
#   neighborhoods possible -- each tissue's testNhoods() call gets a
#   design.df restricted to that tissue's orig.ident samples only, and
#   miloR handles subsetting nhoodCounts() to match internally.
# - Added a min_donors_per_group check per (target, tissue) before
#   testNhoods() -- same lesson already applied in deseq2.R/
#   wgcna_single.R/wgcna_consensus*.R, since nothing guarded against a
#   tissue having too few biological replicates in one group for a stable
#   GLM fit. Checked per tissue (not per target overall), since that's
#   the granularity the actual testNhoods() calls run at.
# - Added buildNhoodGraph(milo) before saving -- the sample script never
#   called this, but it's what populates nhoodGraph(milo), which
#   downstream plotting functions like plotNhoodGraphDA() need. Without
#   it the saved Milo object wouldn't actually be ready for the
#   "downstream plotting" the user asked for.
# - The sample script computed sals_results/c9_results but never saved
#   anything. Saves: the full Milo object (one deliberate exception to
#   "never a whole-object saveRDS()", matching the same reasoning as
#   hdWGCNA's @misc-object exception -- the Milo object is the standard
#   unit of persistence miloR's own downstream plotting functions expect,
#   and at subclustered-population scale it isn't the whole-cohort
#   expression matrix, so this isn't the "duplicate a huge object per
#   task" problem those scripts were avoiding), the Seurat Graph object
#   used to build it (separately, as a plain dgCMatrix, in case anything
#   needs the raw graph without loading miloR/reconstructing a Milo
#   object), and both contrasts' DA results tables as CSVs per tissue.
# - Not implemented, flagged rather than silently added: annotateNhoods().
#   Each target's object is already restricted to one cell type by
#   construction (19_subclustering3.R), so annotating neighborhoods by
#   cell type wouldn't add information; tissue composition per
#   neighborhood might be a genuinely useful annotation given neighborhoods
#   are now shared across tissues, but wasn't explicitly requested -- easy
#   to add later if wanted.

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

data_dir <- paste0(project_root, "/data/milo/", target_name, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0(project_root, "/results/milo/", target_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

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

tissues_present <- sort(unique(as.character(obj$tissue)))
message2(paste0("Tissues found for ", target_name, ": ",
                paste(tissues_present, collapse = ", ")))

min_cells <- 200 # change as needed
if (ncol(obj) < min_cells){
  stop(paste0(target_name, " has only ", ncol(obj), " cells, below ",
              "min_cells = ", min_cells, " -- too few for stable ",
              "neighborhood construction."))
}

# Compute the neighbor graph on the full cross-tissue object --------------
# Not per-tissue -- see header note above. Neighborhoods are meant to be
# shared/comparable across tissues by design.

message2("Computing neighbor graph")

obj <- FindNeighbors(obj,
                     reduction = "harmony",
                     dims = 1:15,
                     k.param = 15,
                     nn.method = "annoy",
                     annoy.metric = "euclidean",
                     compute.SNN = T)

# Build the Milo object, once per target -----------------------------------

message2("Building Milo object")

sce <- as.SingleCellExperiment(obj)

milo <- Milo(sce)
miloR::graph(milo) <- miloR::graph(buildFromAdjacency(obj@graphs$RNA_snn,
                                                      k = 15,
                                                      is.binary = F))

milo <- makeNhoods(milo,
                   prop = 0.05, # change as needed
                   k = 15, # change as needed
                   d = 15,
                   refined = TRUE,
                   reduced_dims = "HARMONY")

p <- plotNhoodSizeHist(milo)
ggsave(p,
       filename = paste0(results_dir, "nhood_size_hist.png"),
       units = "in", dpi = 300,
       height = 5, width = 6)

# orig.ident (donor + tissue), not id -- see header note above.
milo <- countCells(milo,
                   meta.data = as.data.frame(colData(milo)),
                   sample = "orig.ident")

milo <- calcNhoodDistance(milo,
                          d = 15,
                          reduced.dim = "HARMONY")

# Populate the neighborhood graph for downstream plotting -----------------
# The sample script never called this -- see header note above.

milo <- buildNhoodGraph(milo)

# Save the Milo object and the graph used to build it, once per target ----
# See header note above re: saving the whole Milo object.

message2("Saving Milo object and neighbor graph")

saveRDS(milo,
        file = paste0(data_dir, "milo_obj.rds"))

saveRDS(obj@graphs$RNA_snn,
        file = paste0(data_dir, "snn_graph.rds"))

# Test both contrasts, separately per tissue, against the shared neighborhoods
# testNhoods() accepts a design.df that's a strict subset of
# nhoodCounts()'s samples -- see header note above.

min_donors_per_group <- 3 # change as needed

for (tissue_title in tissues_present){

  message2(paste0("Testing ", target_name, " -- ", tissue_title))

  tissue_file <- str_replace_all(tissue_title, " ", "_")

  design <- obj@meta.data %>%
    filter(tissue == tissue_title) %>%
    distinct(orig.ident, group) %>%
    mutate(group = factor(group, levels = c("Control", "sALS", "C9orf72")))

  donors_per_group <- design %>%
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

  design <- as.data.frame(design)
  rownames(design) <- design$orig.ident

  message2("Testing sALS vs. Control")

  sals_results <- testNhoods(milo,
                             design = ~ 0 + group,
                             design.df = design,
                             model.contrasts = "groupsALS - groupControl",
                             fdr.weighting = "graph-overlap",
                             norm.method = "TMM",
                             reduced.dim = "HARMONY")

  write.csv(sals_results,
            file = paste0(results_dir, "sALS_vs_Control_", tissue_file,
                         "_nhood_results.csv"),
            row.names = F)

  message2("Testing C9orf72 vs. Control")

  c9_results <- testNhoods(milo,
                           design = ~ 0 + group,
                           design.df = design,
                           model.contrasts = "groupC9orf72 - groupControl",
                           fdr.weighting = "graph-overlap",
                           norm.method = "TMM",
                           reduced.dim = "HARMONY")

  write.csv(c9_results,
            file = paste0(results_dir, "C9orf72_vs_Control_", tissue_file,
                         "_nhood_results.csv"),
            row.names = F)

}
