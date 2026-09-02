# Runs consensus hdWGCNA between three muscle fiber types (co-expression
# modules found consistently across denervated, Type I, and Type II
# fibers, not fit on one pooled population), sourced from muscle's single
# 17_obj_reassembly.R object. Single-task script, not a SLURM array (see
# jobs/wgcna_consensus_muscle.sh) -- there's exactly one comparison here,
# matching this project's precedent for genuinely single-task work
# (06_obj_reassembly.R has no --array either).
#
# Structurally a mirror image of wgcna_consensus.R (brain vs. spinal cord
# consensus for one fixed cell type each), not an extra task bolted onto
# it: there, cell_type3 was the fixed "which population" identity and
# tissue was the multi-way split; here tissue is fixed (muscle only) and
# cell_type3 is the multi-way split instead. Every cross-cutting
# structural piece (group.by/group_name vs. multi.group.by/multi_groups,
# what's constant vs. what varies for MetacellsByGroups()/the min_cells
# check) has to flip accordingly, which is exactly the kind of thing that
# turns into unreadable branching if crammed into one script -- kept
# separate instead, per the user (a separate script was explicitly fine
# with them if that was cleaner).
#
# Every fix established in wgcna_single.R/wgcna_consensus.R (including
# the two the user found necessary after merging wgcna_consensus.R --
# cell_type3 belongs in MetacellsByGroups()'s group.by even when it looks
# redundant, and FindVariableFeatures() is required right after
# NormalizeData()) is carried forward unchanged:
# - Filters to the three target fiber types before touching raw counts;
#   real raw counts come from data/06_obj_reassembly/bpcells (17's own
#   saved bpcells_data holds normalized data, not counts); muscle's
#   already-fit data/17_obj_reassembly/muscle/harmony.rds is reattached
#   (row-subset to these cells) rather than refit.
# - min_cells is checked per fiber type (not combined), using explicit
#   factor levels so a fiber type with zero matching cells still shows up
#   as a 0 count and fails the check, rather than silently vanishing from
#   the table -- same "abundance filter, checked at the actual comparison
#   granularity" lesson as wgcna_consensus.R's per-tissue check.
# - MetacellsByGroups() groups by c("orig.ident", "cell_type3") -- no
#   "tissue" dimension needed here (constant across this whole script),
#   and cell_type3 is now the dimension that must not get mixed across
#   metacells, matching what tissue was for wgcna_consensus.R.
# - SetMultiExpr()/ModuleConnectivity() use group.by = "tissue",
#   group_name = "Skeletal muscle" (the constant identity, mirroring what
#   cell_type3/the target cell type was in wgcna_consensus.R) with
#   multi.group.by = "cell_type3", multi_groups = fiber_types for the
#   actual three-way comparison. layer = "data" (Seurat v5 naming,
#   already confirmed against hdWGCNA's docs in wgcna_consensus.R).
# - obj[["RNA"]]$data is coerced to a real dgCMatrix right before
#   ModuleConnectivity(), not earlier -- same BPCells/CsparseMatrix
#   reasoning as wgcna_single.R/wgcna_consensus.R.
# - ConstructNetwork() gets the same per-task working-directory isolation
#   to avoid the TOM.rda collision (unfixed hdWGCNA bug,
#   smorabit/hdWGCNA#182); results_dir/data_dir are absolute paths for
#   the same reason.
# - Module scoring uses the same manual AddModuleScore_UCell() +
#   SmoothKNN(reduction = "harmony") approach as the other two scripts.
#   The scores CSV keeps cell_type3 (the fiber type) instead of tissue,
#   since that's the dimension that actually varies here.
# - Saves the same kind of decomposed output (module table, harmonized
#   module eigengenes, smoothed UCell module scores, per-fiber-type
#   soft-power table/plot, dendrogram, KME plot, ME correlogram) under
#   results/wgcna_consensus_muscle/, and only the hdWGCNA
#   @misc[[wgcna_name]] experiment object to
#   data/wgcna_consensus_muscle/wgcna_experiment.rds.

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

fiber_types <- c("Denervated MF", "Type I MF", "Type II MF")
file <- "muscle_fiber_types"

message2(paste0("Processing consensus WGCNA across fiber types: ",
                paste(fiber_types, collapse = ", ")))

# Absolute paths -- this script (like wgcna_single.R/wgcna_consensus.R)
# temporarily setwd()s elsewhere mid-run for the ConstructNetwork()
# TOM-collision workaround below.

data_dir <- paste0(project_root, "/data/wgcna_consensus_muscle/")
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- paste0(project_root, "/results/wgcna_consensus_muscle/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Filter to these fiber types before touching raw counts at all -------------
# See header note above.

message2("Reading in metadata and filtering to these fiber types")

meta <- readRDS("data/17_obj_reassembly/muscle/metadata.rds")
meta_sub <- meta[meta$cell_type3 %in% fiber_types, ]

# Skip if any fiber type is too sparsely represented for stable metacell
# construction -- see header note above.

min_cells <- 200 # change as needed

cell_counts <- table(factor(meta_sub$cell_type3, levels = fiber_types))
if (any(cell_counts < min_cells)){
  stop(paste0("At least one fiber type has too few cells for stable ",
              "metacell construction (min_cells = ", min_cells, "): ",
              paste(names(cell_counts), cell_counts, sep = " = ", collapse = ", ")))
}

message2("Reading in raw counts and Harmony embedding")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta_sub)]

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "RNA")
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

harmony <- readRDS("data/17_obj_reassembly/muscle/harmony.rds")
harmony@cell.embeddings <- harmony@cell.embeddings[rownames(meta_sub), ]
obj[["harmony"]] <- harmony

obj <- ScaleData(obj)

# Identify genes expressed in at least 5% of these cells --------------------

message2("Selecting genes expressed in at least 5% of these cells")

pe <- rowMeans(GetAssayData(obj, layer = "data", assay = "RNA") > 0)
genes_keep <- names(pe)[pe > 0.05] # change this cutoff as needed

# Set up hdWGCNA -------------------------------------------------------
# obj is already filtered to just these three fiber types -- see header
# note above.

message2("Setting up hdWGCNA")

obj <- SetupForWGCNA(obj,
                     gene_select = "custom",
                     features = genes_keep,
                     wgcna_name = file)

message2("Constructing metacells")

# group.by includes cell_type3 (unlike wgcna_consensus.R's tissue-split
# case) so metacells never mix cells across fiber type -- see header note
# above.
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("orig.ident", "cell_type3"),
  reduction = "harmony",
  k = 25, # change as needed
  max_shared = 10, # change as needed
  ident.group = "cell_type3"
)

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj, features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj, features = VariableFeatures(obj))
obj <- RunHarmonyMetacells(obj, group.by.vars = "orig.ident")

# group.by/group_name reference "tissue" (constant = "Skeletal muscle"
# for every cell here) as the fixed identity -- mirroring what
# cell_type3/the target cell type was in wgcna_consensus.R -- while
# multi.group.by/multi_groups do the actual three-way fiber-type
# comparison. See header note above.
obj <- SetMultiExpr(
  obj,
  group_name = "Skeletal muscle",
  group.by = "tissue",
  multi.group.by = "cell_type3",
  multi_groups = fiber_types,
  assay = "RNA",
  layer = "data",
  use_metacells = T
)

# Find soft power per fiber type ---------------------------------------

message2("Testing soft powers")

obj <- TestSoftPowersConsensus(obj)

plot_list <- PlotSoftPowers(obj)

p_list <- lapply(seq_along(fiber_types), function(i){
  plot_list[[i]][[1]] +
    ggtitle(paste0("Fiber type: ", fiber_types[i])) +
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
# wgcna_consensus.R.
#
# Same TOM.rda SLURM-array collision as wgcna_single.R/wgcna_consensus.R
# -- see wgcna_single.R's header for the full diagnosis (unfixed hdWGCNA
# bug, smorabit/hdWGCNA#182). Same workaround: isolate the working
# directory for just this call. (This script isn't itself an array job,
# but the workaround is harmless and kept for consistency in case this
# ever runs alongside other wgcna_* jobs sharing the project's setwd().)

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
PlotDendrogram(obj, main = "Muscle fiber type consensus dendrogram")
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

obj <- ModuleConnectivity(obj, group_name = "Skeletal muscle", group.by = "tissue")

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
# Same approach as wgcna_single.R/wgcna_consensus.R.

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
  dplyr::select(orig.ident, group, cell_type3, matches("_UCell_kNN$"))

write.csv(scores,
          file = paste0(results_dir, "module_scores_ucell.csv"),
          row.names = F)

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
