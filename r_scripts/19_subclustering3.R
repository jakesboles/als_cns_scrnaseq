# Subclusters and re-integrates specific cell types/cell-type groups that
# appear across multiple tissues (e.g. microglia from both brain and
# spinal cord), using 17_obj_reassembly.R's per-tissue metadata as the
# source, as input for MiloR, consensus hdWGCNA, and related downstream
# figures. Stops at Harmony integration plus one diagnostic UMAP per
# group -- no reclustering or marker finding. Runs as a SLURM job array
# (see jobs/19_subclustering3.sh), one task per entry in the hardcoded
# subclustering_targets list below.

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

fam <- c("tissue", "tissue", "tissue", "cell_type3", "cell_type3", "cell_type3")
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

# Write FindAllMarkers() output for grouping of interest ------------------

Idents(obj) <- fam[task_id]

markers <- FindAllMarkers(obj)

write.csv(markers,
          file = paste0(results_dir, "fam_output.csv"),
          row.names = F)
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
