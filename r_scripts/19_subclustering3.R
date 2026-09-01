# Subclusters and re-integrates specific cell types that appear across
# multiple tissues (e.g. microglia from both brain and spinal cord),
# using 17_obj_reassembly.R's per-tissue metadata as the source. Runs as
# a SLURM job array (see jobs/19_subclustering3.sh), one task per entry
# in the hardcoded subclustering_targets list below -- targets are
# heterogeneous (different tissue combinations and cell_type3 sets per
# group) and were described inline by the user rather than as a params
# file, matching this project's established convention for that case
# (see 15_subclustering2.R's subclustering_targets).
#
# Outputs feed MiloR, a consensus hdWGCNA, and/or figures -- not another
# annotation round, so unlike 13/15 this script does NOT recluster or run
# FindAllMarkers(); it stops at Harmony integration plus one diagnostic
# UMAP per group, matching 17_obj_reassembly.R/18_full_integration.R's
# scope instead (confirmed with the user).
#
# Design notes:
# - Each target lists the tissue(s) it needs; that tissue's
#   data/17_obj_reassembly/<tissue>/metadata.rds is read and (for multi-
#   tissue targets) bind_rows()'d together -- the same per-tissue
#   metadata-concatenation pattern 18_full_integration.R uses, just
#   scoped to whichever tissues a given target actually needs, rather
#   than routing every target through 18's saved (all_tissues/brain_sc)
#   output. That's a deliberate simplification: every target here rebuilds
#   from real raw counts and refits its own Harmony integration on just
#   its own cells anyway (see next point), so 18's own Harmony fit over a
#   much larger, more heterogeneous population isn't actually used by
#   anything downstream -- only the tissue metadata is needed, and
#   18_full_integration.R's saved artifacts would be an unnecessary
#   dependency for that.
# - Real raw counts are pulled from data/06_obj_reassembly/bpcells (not
#   17_obj_reassembly.R's own saved bpcells_data, which holds normalized
#   data, not counts) -- the same pipeline-wide gotcha as every other
#   script that reruns FindVariableFeatures() from scratch (13/15/17/18).
# - RunPCA() uses npcs = 50, matching 13/15's cell-type/group-subset scale
#   (not 07/17/18's npcs = 100 full-tissue scale) -- every target here is
#   a single cell class or small group of related classes, smaller and
#   more homogeneous than a whole tissue.
# - cell_type3 vocabulary for each target was confirmed against the
#   user's own DimPlot legends (brain_sc and all_tissues) rather than
#   guessed:
#   - "All neuron classes" = every EN/IN subtype label plus MN/SN.
#   - "Myeloid cells" = Microglia/Macrophage/Monocyte/Neutrophil/Mast
#     cell specifically -- NOT "Proliferating myeloid cell", per the
#     user's explicit enumeration.
#   - "Muscle fiber types" = every cell_type3 label containing "MF"
#     (Denervated MF, Proliferating MF, Type I MF, Type II MF).
# - RunUMAP() here is a fresh, unmodified default (no uwot.init override)
#   -- these are much smaller, more homogeneous per-cell-type populations
#   than 18_full_integration.R's whole-tissue objects, and far less likely
#   to reproduce the near-disconnected-graph pathology that forced
#   uwot.init = "pca" there. If a target segfaults the same way, that's
#   the fix to reach for first (see 18_full_integration.R's header for
#   the full diagnosis).

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

subclustering_targets <- list(
  list(name = "microglia", tissues = c("brain", "sc"),
       cell_types = "Microglia"),
  list(name = "astrocyte", tissues = c("brain", "sc"),
       cell_types = "Astrocyte"),
  list(name = "oligodendrocyte", tissues = c("brain", "sc"),
       cell_types = "Oligodendrocyte"),
  list(name = "neurons", tissues = c("brain", "sc"),
       cell_types = c("EN", "L2-3 EN", "L4 EN", "L5 ET EN", "L5 IT EN",
                      "L5-6 NP EN", "L6 CT EN", "L6 IT EN", "L6b EN",
                      "IN", "COL15A1 IN", "CXCL14 IN", "LAMP5 IN",
                      "NPY IN", "PVALB IN", "RELN IN", "SST IN", "VIP IN",
                      "MN", "SN")),
  list(name = "muscle_fiber", tissues = "muscle",
       cell_types = c("Denervated MF", "Proliferating MF", "Type I MF",
                      "Type II MF")),
  list(name = "myeloid", tissues = c("brain", "sc", "muscle"),
       cell_types = c("Microglia", "Macrophage", "Monocyte", "Neutrophil",
                      "Mast cell"))
)

# Figure out which target this task handles ----------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/19_subclustering3.sh), one task per ",
       "entry in subclustering_targets, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(subclustering_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(subclustering_targets), " subclustering targets -- ",
              "check the --array range in jobs/19_subclustering3.sh."))
}

target <- subclustering_targets[[task_id]]
target_name <- target$name
target_tissues <- target$tissues
target_cell_types <- target$cell_types

message2(paste0("Processing ", target_name, " (",
                paste(target_tissues, collapse = ", "),
                "), task ", task_id, "/", length(subclustering_targets)))

data_dir <- paste0("data/19_subclustering3/", target_name, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0("results/19_subclustering3/", target_name, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Load metadata for each tissue this target needs, then subset to the
# target's cell types --------------------------------------------------
# Same per-tissue metadata concatenation pattern as
# 18_full_integration.R, scoped to only the tissue(s) this specific
# target needs -- see header note above.

message2("Reading in tissue metadata")

meta_list <- list()
for (t in target_tissues){
  meta_list[[t]] <- readRDS(paste0("data/17_obj_reassembly/", t, "/metadata.rds"))
}
meta_all <- bind_rows(meta_list)

meta_sub <- meta_all[meta_all$cell_type3 %in% target_cell_types, ]

if (nrow(meta_sub) == 0){
  stop(paste0("No cells matched cell_type3 %in% c(",
              paste(target_cell_types, collapse = ", "),
              ") across tissue(s) ", paste(target_tissues, collapse = ", "),
              " -- check subclustering_targets in this script for a typo ",
              "against the actual cell_type3 labels."))
}

# Rebuild from real raw counts for this cell/tissue selection ---------------
# FindVariableFeatures()'s default "vst" method needs real counts -- see
# header note above.

message2("Loading raw counts for this cell/tissue selection")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta_sub)]

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "RNA")

# Normalize, find variable features, scale, and run PCA ----------------------

message2("Normalizing, finding variable features, scaling, and running PCA")

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 50)

message2("Integrating samples using Harmony")

obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)

obj <- IntegrateLayers(obj,
                       method = "HarmonyIntegration",
                       orig.reduction = "pca",
                       new.reduction = "harmony",
                       dims = 1:20)

obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# Compute a diagnostic UMAP -----------------------------------------------
# Integration + UMAP only, no reclustering -- see header note above.

message2("Computing UMAP")

obj <- RunUMAP(obj,
               umap.method = "uwot",
               reduction = "harmony",
               dims = 1:20,
               metric = "euclidean",
               min.dist = 0.5,
               n.neighbors = 30L,
               reduction.name = "harmony_umap",
               return.model = F)

# Diagnostic DimPlots ---------------------------------------------------

message2("Making diagnostic DimPlots")

for (group in c("cell_type3", "tissue", "batch", "orig.ident", "group")){
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
