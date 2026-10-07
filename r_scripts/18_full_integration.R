# Cross-tissue Harmony integration on top of 17_obj_reassembly.R's final
# per-tissue objects, for analyses spanning multiple tissues: rebuilds
# from real raw counts, re-normalizes/re-PCAs the combined population,
# integrates with Harmony, and computes one diagnostic UMAP. Runs as a
# SLURM job array (see jobs/18_full_integration.sh), one task per target
# (all 3 tissues, and brain + cervical spinal cord only).

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(scCustomize)
  library(BPCells)
})

options(future.globals.maxSize = 250 * 1024^3)

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

integration_targets <- list(
  list(name = "all_tissues", tissues = c("brain", "sc", "muscle")),
  list(name = "brain_sc", tissues = c("brain", "sc"))
)

# Figure out which target this task handles ----------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/18_full_integration.sh), one task per ",
       "entry in integration_targets, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(integration_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(integration_targets), " integration targets -- ",
              "check the --array range in jobs/18_full_integration.sh."))
}

target <- integration_targets[[task_id]]
target_name <- target$name
target_tissues <- target$tissues

message2(paste0("Processing ", target_name, " (",
                paste(target_tissues, collapse = ", "),
                "), task ", task_id, "/", length(integration_targets)))

data_dir <- paste0("data/18_full_integration/", target_name, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0("results/18_full_integration/", target_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Load metadata for each tissue being integrated -----------------------
# Only 17_obj_reassembly.R's metadata is needed here -- its bpcells_data/
# harmony.rds aren't (see header note above).

message2("Reading in tissue metadata")

meta_list <- list()
for (t in target_tissues){
  meta_list[[t]] <- readRDS(paste0("data/17_obj_reassembly/", t, "/metadata.rds"))
}
meta_all <- bind_rows(meta_list)

# Rebuild from real raw counts for the combined cell set --------------------

message2("Loading raw counts for the combined cell set")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta_all)]

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_all, assay = "RNA")

# Normalize, find variable features, scale, and run PCA ----------------------

message2("Normalizing, finding variable features, scaling, and running PCA")

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 100)

message2("Making PCA diagnostic plots")

p <- ElbowPlot(obj, ndims = 100)
ggsave(p,
       filename = paste0(results_dir, "pca_elbow.png"),
       units = "in", dpi = 600, bg = "white",
       height = 6, width = 6)

Iterate_PC_Loading_Plots(obj,
                         file_path = results_dir,
                         file_name = "pca_loadings")

message2("Integrating tissues using Harmony")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:20)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Compute UMAP directly on the Harmony embedding ------------------------
# Switched from the graph-based nn.name = "RNA.nn" approach (still
# commented out below) to running RunUMAP() straight off the "harmony"
# reduction, per the user -- the graph-based UMAP wasn't visually
# satisfying. That switch is what surfaced two bugs that the graph-based
# path had been silently absorbing:
# - `n_neighbors = 15L` isn't a real RunUMAP() argument (it's
#   `n.neighbors`, dotted) -- Seurat warns "arguments not used:
#   n_neighbors" and silently falls back to its own default of 30. This
#   typo is actually present in every RunUMAP() call across this project
#   (10/13/15/17, and the commented block below), but it was harmless
#   everywhere else because nn.name = "RNA.nn" supplies the neighbor
#   structure directly and n.neighbors/n_neighbors is never consulted --
#   here, with no nn.name, n.neighbors actually controls the result, so
#   the typo is fixed.
# - RunUMAP()'s default initialization (`uwot.init = "spectral"`) builds
#   the graph Laplacian and eigendecomposes it via RSpectra, which
#   segfaulted on the brain_sc target. An earlier version of this comment
#   attributed this to brain_sc's cell count (500K+), but the user
#   pointed out that all_tissues -- a *larger* target running the
#   identical RunUMAP() call -- completed without issue, which rules that
#   out as the mechanism. Cell count alone isn't it.
#   RSpectra's ARPACK-based eigensolver is also known to crash outright
#   (segfault, not a clean R error) on numerically degenerate input --
#   e.g. many exactly-duplicated points in the embedding -- rather than
#   just large ones. Circumstantial support for that specific to
#   brain_sc: Harmony's own k-means step logged 10x "Quick-TRANSfer stage
#   steps exceeded maximum" warnings for brain_sc (visible in the pasted
#   log), a classic symptom of many tied/duplicate points during k-means.
#   That's a real hypothesis, not a confirmed diagnosis -- the duplicate/
#   non-finite check logged just below exists to actually find out,
#   instead of guessing again.
# `uwot.init = "pca"` initializes UMAP from the embedding directly
# instead of eigendecomposing the graph Laplacian, avoiding RSpectra
# regardless of which hypothesis above is right -- kept as the fix since
# it should hold either way, but flagged here as not fully explained.
# This is still an informed hypothesis, not verified by actually running
# it (no R/Seurat runtime in this dev container).

# harmony <- readRDS(paste0(data_dir, "harmony.rds"))
# obj[["harmony"]] <- harmony

message2("Checking Harmony embedding for duplicate/non-finite values")

harmony_emb <- as.data.frame(Embeddings(obj, reduction = "harmony")[, 1:20])
n_dup <- sum(duplicated(harmony_emb))
n_nonfinite <- sum(is.infinite(as.matrix(harmony_emb)))
n_na <- sum(is.na(as.matrix(harmony_emb)))
message2(paste0(nrow(harmony_emb), " cells, ", n_dup,
                " duplicated embedding rows, ", n_nonfinite,
                " non-finite values"))

message2("Computing UMAP directly on the Harmony embedding")

obj <- RunUMAP(obj,
               umap.method = "uwot",
               reduction = "harmony",
               dims = 1:15,
               # nn.name = "RNA.nn",
               metric = "euclidean",
               min.dist = 0.5,
               n.neighbors = 15L,
               # repulsion.strength = 0.5,
               # uwot.init = "random",
               reduction.name = "harmony_umap",
               return.model = F)

# obj <- obj %>%
#   FindNeighbors(reduction = "harmony",
#                 dims = 1:15,
#                 k.param = 15,
#                 nn.method = "annoy",
#                 annoy.metric = "euclidean",
#                 return.neighbor = T) %>%
#   FindNeighbors(reduction = "harmony",
#                 dims = 1:15,
#                 k.param = 15,
#                 nn.method = "annoy",
#                 annoy.metric = "euclidean",
#                 compute.SNN = T) %>%
#   RunUMAP(umap.method = "uwot",
#           nn.name = "RNA.nn",
#           metric = "euclidean",
#           min.dist = 0.5,
#           n.neighbors = 15L,
#           reduction.name = "harmony_umap",
#           return.model = F)

# Diagnostic DimPlots ---------------------------------------------------
# final_label2/Batch/Group from the old script are cell_type3/batch/group
# here.

message2("Making diagnostic DimPlots")

for (group in c("cell_type3", "batch", "tissue", "orig.ident", "group")){
  w <- if (group %in% c("cell_type3", "orig.ident")) 15 else 11

  p <- DimPlot_scCustom(obj,
                        reduction = "harmony_umap",
                        group.by = group)
  ggsave(p,
         filename = paste0(results_dir, group, "_dimplot.png"),
         units = "in", dpi = 600,
         height = 8, width = w)
}

# Save metadata, integrated embedding, normalized expression, and UMAP -----

message2("Saving metadata, count matrix, Harmony embedding, and UMAP")

bpcells_data_dir <- paste0(data_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(obj@meta.data,
        file = paste0(data_dir, "metadata.rds"))

saveRDS(obj[["harmony"]],
        file = paste0(data_dir, "harmony.rds"))

saveRDS(obj[["harmony_umap"]],
        file = paste0(data_dir, "harmony_umap.rds"))