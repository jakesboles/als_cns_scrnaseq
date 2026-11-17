# Tests module preservation (NetRep, via hdWGCNA::ModulePreservationNetRep())
# across every directed pair among 5 WGCNA results: Microglia in spinal cord,
# brain, and consensus (brain+spinal cord), plus muscle Macrophages and
# Monocytes. wgcna_single.R/wgcna_consensus_cns.R only save each run's
# @misc[[wgcna_name]] experiment list, not a full Seurat object, so each of
# the 5 targets is rebuilt here from its own source script's exact recipe
# (real counts + metadata + Harmony) before reattaching its saved experiment
# and calling SetDatExpr(). Runs as a 20-task SLURM job array (see
# jobs/module_preservation.sh), one task per (reference, query) pair.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(tidyverse)
  library(NetRep)
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
enableWGCNAThreads(nThreads = 16)

# The 5 WGCNA targets to compare, every directed pair ----------------------
# Each entry names: where its saved hdWGCNA experiment (@misc[[wgcna_name]]
# from wgcna_single.R/wgcna_consensus_cns.R) lives, which metadata/Harmony
# source object it was built from, and which cell_type3 value it was
# filtered to -- everything build_wgcna_obj() needs to reconstruct it below.

targets <- list(
  Microglia_sc = list(
    wgcna_path = "data/wgcna_single/sc/Microglia/wgcna_experiment.rds",
    metadata_path = "data/17_obj_reassembly/sc/metadata.rds",
    harmony_path = "data/17_obj_reassembly/sc/harmony.rds",
    cell_type = "Microglia",
    display_name = "Microglia (spinal cord)"
  ),
  Microglia_brain = list(
    wgcna_path = "data/wgcna_single/brain/Microglia/wgcna_experiment.rds",
    metadata_path = "data/17_obj_reassembly/brain/metadata.rds",
    harmony_path = "data/17_obj_reassembly/brain/harmony.rds",
    cell_type = "Microglia",
    display_name = "Microglia (brain)"
  ),
  Microglia_consensus = list(
    wgcna_path = "data/wgcna_consensus/Microglia/wgcna_experiment.rds",
    metadata_path = "data/18_full_integration/brain_sc/metadata.rds",
    harmony_path = "data/18_full_integration/brain_sc/harmony.rds",
    cell_type = "Microglia",
    display_name = "Microglia (consensus)"
  ),
  Macrophage_muscle = list(
    wgcna_path = "data/wgcna_single/muscle/Macrophage/wgcna_experiment.rds",
    metadata_path = "data/17_obj_reassembly/muscle/metadata.rds",
    harmony_path = "data/17_obj_reassembly/muscle/harmony.rds",
    cell_type = "Macrophage",
    display_name = "Macrophage (muscle)"
  ),
  Monocyte_muscle = list(
    wgcna_path = "data/wgcna_single/muscle/Monocyte/wgcna_experiment.rds",
    metadata_path = "data/17_obj_reassembly/muscle/metadata.rds",
    harmony_path = "data/17_obj_reassembly/muscle/harmony.rds",
    cell_type = "Monocyte",
    display_name = "Monocyte (muscle)"
  )
)

pair_targets <- expand.grid(ref = names(targets), test = names(targets),
                            stringsAsFactors = F)
pair_targets <- pair_targets[pair_targets$ref != pair_targets$test, ]
pair_targets <- pair_targets[order(pair_targets$ref, pair_targets$test), ]
rownames(pair_targets) <- NULL

# Figure out which directed pair this task handles --------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/module_preservation.sh), one task per row ",
       "of the in-script pair_targets table (every directed pair among the ",
       "5 WGCNA targets), not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(pair_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(pair_targets), " directed pairs in pair_targets -- check ",
              "the --array range in jobs/module_preservation.sh."))
}

ref_name <- pair_targets$ref[task_id]
test_name <- pair_targets$test[task_id]

message2(paste0("Processing ", ref_name, " -> ", test_name, ", task ",
                task_id, "/", nrow(pair_targets)))

results_dir <- paste0(project_root, "/results/module_preservation/",
                      ref_name, "_into_", test_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Rebuild one of the 5 targets as a full Seurat object -----------------------
# Reproduces wgcna_single.R's/wgcna_consensus_cns.R's own object-build recipe
# exactly (real raw counts from 06_obj_reassembly, that target's metadata/
# Harmony, Normalize/FindVariableFeatures/ScaleData) -- see header note
# above for why this reconstruction is necessary here. role is "ref" or
# "test", matching ModulePreservationNetRep()'s own reference/query roles;
# the saved experiment is reattached under that name instead of its
# original wgcna_name so both objects can be told apart.

build_wgcna_obj <- function(target_name, role){

  t <- targets[[target_name]]

  meta <- readRDS(t$metadata_path)
  meta_sub <- meta[meta$cell_type3 == t$cell_type, ]

  raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
  raw_mat <- raw_mat[, rownames(meta_sub)]

  obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "RNA")
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)

  harmony <- readRDS(t$harmony_path)
  harmony@cell.embeddings <- harmony@cell.embeddings[rownames(meta_sub), ]
  obj[["harmony"]] <- harmony

  obj <- ScaleData(obj)

  wgcna_experiment <- readRDS(t$wgcna_path)
  obj@misc[[role]] <- wgcna_experiment
  obj <- SetActiveWGCNA(obj, role)

  obj <- SetDatExpr(
    obj,
    group_name = t$cell_type,
    group.by = "cell_type3",
    assay = "RNA",
    layer = "data",
    use_metacells = T
  )

  obj
}

message2(paste0("Rebuilding reference object: ", ref_name))
obj_ref <- build_wgcna_obj(ref_name, "ref")

message2(paste0("Rebuilding query object: ", test_name))
obj_test <- build_wgcna_obj(test_name, "test")

# Module preservation testing -----------------------------------------------

message2(paste0("Testing module preservation: ", ref_name, " -> ", test_name))

obj_test <- ModulePreservationNetRep(
  obj_test,
  seurat_ref = obj_ref,
  name = "test",
  n_permutations = 10000,
  TOM_use = "test",
  n_threads = 16,
  wgcna_name = "test",
  wgcna_name_ref = "ref"
)

results <- GetModulePreservation(obj_test, "test", "test")

write.csv(results$p.values,
          file = paste0(results_dir, "pvals.csv"))
write.csv(results$observed,
          file = paste0(results_dir, "observed.csv"))

p <- PlotModulePreservationLollipop(obj_test,
                                    name = "test",
                                    features = "average",
                                    wgcna_name = "test") +
  ggtitle(paste0(targets[[ref_name]]$display_name, " projected\nto ",
                targets[[test_name]]$display_name))

ggsave(p,
       filename = paste0(results_dir, "lollipop.png"),
       units = "in", dpi = 600,
       height = nrow(results$p.values) * 0.7, width = 5)
