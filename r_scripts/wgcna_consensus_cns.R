# Runs consensus hdWGCNA (co-expression modules found consistently
# between brain and spinal cord, not fit on one pooled population) for
# three cell types shared across both tissues: microglia,
# oligodendrocyte, astrocyte. Saves module tables, harmonized module
# eigengenes, smoothed UCell module scores, and per-tissue soft-power/
# network diagnostics. Runs as a 3-task SLURM job array (see
# jobs/wgcna_consensus.sh), one task per subclustering_targets entry.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(UCell)
  library(cowplot)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

project_root <- "/projects/b1169/boles/als_cns_scrnaseq"
setwd(project_root)

theme_set(theme_cowplot())
set.seed(256)
enableWGCNAThreads(nThreads = 16)

# Figure out which cell type this task handles ------------------------

subclustering_targets <- c("Microglia", "Oligodendrocyte", "Astrocyte")
tissue_groups <- c("Motor cortex", "Cervical spinal cord")

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/wgcna_consensus.sh), one task per entry ",
       "in subclustering_targets, not as a standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > length(subclustering_targets)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              length(subclustering_targets), " subclustering targets -- ",
              "check the --array range in jobs/wgcna_consensus.sh."))
}

cell_type_target <- subclustering_targets[task_id]
file <- str_replace_all(cell_type_target, " ", "_")

message2(paste0("Processing ", cell_type_target, ", task ", task_id, "/",
                length(subclustering_targets)))

# Absolute paths -- this script (like wgcna_single.R) temporarily
# setwd()s elsewhere mid-run for the ConstructNetwork() TOM-collision
# workaround below.

data_dir <- paste0(project_root, "/data/wgcna_consensus/", file, "/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0(project_root, "/results/wgcna_consensus/", file, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Filter to this cell type before touching raw counts at all ----------------
# See header note above.

message2("Reading in metadata and filtering to this cell type")

meta <- readRDS("data/18_full_integration/brain_sc/metadata.rds")
meta_sub <- meta[meta$cell_type3 == cell_type_target, ]

# Skip cell types too sparsely represented, per tissue, for stable
# per-tissue metacell construction -- see header note above.

min_cells <- 200 # change as needed

cell_counts <- table(droplevels(factor(meta_sub$tissue)))
if (any(cell_counts < min_cells)){
  stop(paste0(cell_type_target, " has too few cells in at least one tissue ",
              "for stable per-tissue metacell construction (min_cells = ",
              min_cells, "): ",
              paste(names(cell_counts), cell_counts, sep = " = ", collapse = ", ")))
}

message2("Reading in raw counts and Harmony embedding")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta_sub)]

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "RNA")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

harmony <- readRDS("data/18_full_integration/brain_sc/harmony.rds")
harmony@cell.embeddings <- harmony@cell.embeddings[rownames(meta_sub), ]
obj[["harmony"]] <- harmony

obj <- ScaleData(obj)

# Identify genes expressed in at least 5% of this cell type's cells ---------

message2("Selecting genes expressed in at least 5% of this cell type")

pe <- rowMeans(GetAssayData(obj, layer = "data", assay = "RNA") > 0)
genes_keep <- names(pe)[pe > 0.05] # change this cutoff as needed

# Set up hdWGCNA -------------------------------------------------------
# obj is already filtered to just this cell type -- see header note above.

message2("Setting up hdWGCNA")

obj <- SetupForWGCNA(obj,
                     gene_select = "custom",
                     features = genes_keep,
                     wgcna_name = file)

message2("Constructing metacells")

# group.by includes tissue (unlike wgcna_single.R) so metacells never mix
# cells across tissue -- see header note above. cell_type3 is dropped
# from group.by (constant post-filtering, same simplification as the
# single script) but kept as ident.group.
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("orig.ident", "tissue", "cell_type3"),
  reduction = "harmony",
  k = 25, # change as needed
  max_shared = 10, # change as needed
  ident.group = "cell_type3"
)

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj, features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj, features = VariableFeatures(obj))
obj <- RunHarmonyMetacells(obj, group.by.vars = "orig.ident")

obj <- SetMultiExpr(
  obj,
  group_name = cell_type_target,
  group.by = "cell_type3",
  multi.group.by = "tissue",
  multi_groups = tissue_groups,
  assay = "RNA",
  layer = "data",
  use_metacells = T
)

# Find soft power per tissue -----------------------------------------------

message2("Testing soft powers")

obj <- TestSoftPowersConsensus(obj)

plot_list <- PlotSoftPowers(obj)

p_list <- lapply(seq_along(tissue_groups), function(i){
  plot_list[[i]][[1]] +
    ggtitle(paste0("Tissue: ", tissue_groups[i])) +
    theme(plot.title = element_text(hjust = 0.5))
})
p <- wrap_plots(p_list, ncol = 2)
ggsave(p,
       filename = paste0(results_dir, "soft_power.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

power_table <- GetPowerTable(obj)
write.csv(power_table,
          file = paste0(results_dir, "soft_powers.csv"),
          row.names = F)

# Build consensus TOM and cluster genes into modules -------------------------
# Letting ConstructNetwork() pick the soft power automatically, matching
# the draft.
#
# Same TOM.rda SLURM-array collision as wgcna_single.R -- see that
# script's header for the full diagnosis (unfixed hdWGCNA bug,
# smorabit/hdWGCNA#182). Same workaround: isolate the working directory
# for just this call.

message2("Constructing consensus network")

tom_dir <- paste0(data_dir, "tom/")
dir.create(tom_dir, showWarnings = F, recursive = T)

setwd(tom_dir)
tryCatch({
  obj <- ConstructNetwork(obj,
                          tom_name = file,
                          consensus = T,
                          overwrite_tom = T)
}, finally = {
  setwd(project_root)
})

png(paste0(results_dir, "dendrogram.png"),
    height = 8, width = 8, units = "in", res = 600)
PlotDendrogram(obj, main = paste0(cell_type_target, " dendrogram"))
dev.off()

# Module stats ------------------------------------------------------------

message2("Computing module eigengenes and connectivity")

obj <- SetActiveWGCNA(obj, file)
obj <- ModuleEigengenes(obj, group.by.vars = "orig.ident")

# ModuleConnectivity() reaches back into the full single-cell "data"
# layer for its corSparse()-based correlation step -- see wgcna_single.R's
# header for the full BPCells/CsparseMatrix diagnosis. Deferred to right
# before this call, not earlier, for the same reason.
obj[["RNA"]]$data <- as(obj[["RNA"]]$data, "dgCMatrix")

obj <- ModuleConnectivity(obj, group_name = cell_type_target, group.by = "cell_type3")

p <- PlotKMEs(obj, ncol = 4, text_size = 4)
ggsave(p,
       filename = paste0(results_dir, "module_connectivity.png"),
       units = "in", dpi = 600,
       height = 12, width = 12)

mods <- obj@misc[[file]][["wgcna_modules"]]
write.csv(mods,
          file = paste0(results_dir, "modules.csv"),
          row.names = F)

# Module expression scores via UCell -----------------------------------------
# Same approach as wgcna_single.R -- see header note above.

message2("Scoring modules with UCell")

module_names <- setdiff(unique(mods$module), "grey")

gene_sets <- lapply(module_names, function(m){
  mods$gene_name[mods$module == m]
})
names(gene_sets) <- module_names

maxrank <- max(lengths(gene_sets))

obj <- AddModuleScore_UCell(obj, features = gene_sets, maxRank = maxrank)
obj <- SmoothKNN(obj,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "harmony")

scores <- obj@meta.data %>%
  dplyr::select(orig.ident, group, tissue, matches("_UCell_kNN$"))

write.csv(scores,
          file = paste0(results_dir, "module_scores_ucell.csv"),
          row.names = T)

# Harmonized module eigengenes -----------------------------------------------

message2("Saving module eigengenes")

hMEs <- GetMEs(obj, harmonized = T)
write.csv(hMEs,
          file = paste0(results_dir, "module_eigengenes.csv"))

png(paste0(results_dir, "module_eigengene_correlogram.png"),
    height = 8, width = 8, units = "in", res = 600)
ModuleCorrelogram(obj, features = "MEs")
dev.off()

# Save the hdWGCNA experiment for further downstream use ---------------------
# Just the hdWGCNA network/module state (@misc[[wgcna_name]]), not the
# whole Seurat object -- see header note above.

message2("Saving hdWGCNA experiment object")

wgcna_experiment <- obj@misc[[file]]
saveRDS(wgcna_experiment,
        file = paste0(data_dir, "wgcna_experiment.rds"))
