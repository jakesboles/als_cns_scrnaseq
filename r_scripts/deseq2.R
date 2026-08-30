# Pseudobulk DESeq2 differential expression within each cell type
# (cell_type3) of one tissue, comparing sALS vs Control and C9orf72 vs
# Control. Runs as a SLURM job array (see jobs/deseq2.sh), one task per
# tissue, loading each tissue's final annotated object from
# 17_obj_reassembly.R.
#
# Design notes (fixed while reviewing the user's draft):
# - 17_obj_reassembly.R's saved bpcells_data holds normalized data, not
#   real counts (same pipeline-wide gotcha as 09_integration onward --
#   see 13/17's header notes). AggregateExpression() needs real raw
#   counts to sum correctly, so this script pulls raw counts from
#   data/06_obj_reassembly/bpcells directly, subset to the cell barcodes
#   retained in 17's metadata.rds, same pattern 13/17 use.
# - AggregateExpression()'s counts-layer argument is `layer =` under
#   Seurat v5, not `slot =` (the draft used the pre-v5 name).
# - results()'s independent-filtering toggle is `independentFiltering =`
#   (the draft had `independentFilter`, missing "ing" -- would have
#   errored as an unused argument). Combining independentFiltering = T
#   with filterFun = ihw is the documented IHW usage.
# - The draft's `meta$Genotype <- factor(..., levels = c("WT","HET","HMM"))`
#   block looks like leftover copy-paste from an unrelated mouse study --
#   "Genotype" was never selected into meta, and WT/HET/HMM aren't levels
#   that exist in this human ALS cohort (whose grouping variable is
#   `group`: Control/sALS/C9orf72). Dropped.
# - The draft's second results() contrast used c("group", "C9orf72",
#   "ALS") -- "ALS" isn't a valid factor level. Fixed to c("group",
#   "C9orf72", "Control"), matching the user's stated "C9-ALS vs
#   control" comparison and the already-correctly-named paired
#   lfcShrink(coef = "group_C9orf72_vs_Control") call.
# - The draft selected `sex` and `age` into meta as covariates, but
#   neither column exists anywhere in this project's single-cell
#   metadata pipeline (confirmed against 00_cellbender_plotting.R's
#   explicit column list). Per the user, both are needed for this
#   analysis, so they're now joined in from the separate
#   tab_data/target_als_demographics_compiled.csv (donor-level, keyed by
#   `case_number`, used previously only by demographics_figure.R) onto
#   `id` (the donor-only label added in 02_qc1.R, orig.ident minus the
#   tissue suffix). `case_number` and `id` disagree on hyphen placement
#   for at least the Barrow ("GWF") site -- 02_qc1.R has a commented-out
#   block recoding "GWF-19-47" -> "GWF19-47" for exactly this reason --
#   so the join key strips all hyphens from both sides rather than
#   hardcoding the known mismatches. NOTE: the demographics CSV isn't
#   readable from this dev container (it's an ungitted reference input,
#   only present on the HPC filesystem), so the exact column name for
#   sex after janitor::clean_names() is unverified -- `age_at_death` is
#   confirmed correct (demographics_figure.R already uses it), but
#   double-check the `sex` column name against the actual file and fix
#   the dplyr::select() below if it differs. Both covariates are added to
#   the DESeq2 design (design = ~ sex + age_at_death + group) since the
#   user called them out as important for this analysis; the group
#   contrasts/coefficient names below are unaffected by this since group
#   stays a separate term with the same factor levels.
# - Output directory restructured to results/deseq2/<tissue>/, matching
#   every other script's convention; the draft had no tissue segment and
#   several typo'd path variables (resuls_dir, bare "results").
# - The draft also had no object-loading or SLURM array setup at all (the
#   "load BPCells object" step was left as a comment, `tissue` was
#   referenced in a message() but never defined, and `celltypes[i]` was
#   used once instead of the actually-defined `cell_types[i]`) -- all
#   added/fixed here, following the tissues-table + SLURM_ARRAY_TASK_ID
#   fail-fast pattern used by every other per-tissue array script
#   (09/17/18).
# - The draft matched pseudobulk sample metadata to AggregateExpression()'s
#   output via `match(colnames(exp), rownames(meta))`, but rownames(meta)
#   are still per-cell barcodes at that point (meta was built from
#   sub@meta.data, not yet collapsed to one row per sample), while
#   colnames(exp) are sample names (orig.ident) -- that match() would
#   never hit and every sample's covariates would come out NA. Fixed to
#   match on the orig.ident column directly and set matching rownames,
#   since DESeqDataSetFromMatrix() requires colData's rownames to equal
#   countData's colnames exactly.

suppressMessages({
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(IHW)
  library(Seurat)
  library(BPCells)
  library(janitor)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

tissues <- data.frame(
  file = c("brain", "sc", "muscle"),
  title = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle")
)

# Figure out which tissue this task handles ---------------------------------

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if (task_id == ""){
  stop("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run as a ",
       "SLURM job array (see jobs/deseq2.sh), one task per tissue, not as a ",
       "standalone Rscript call.")
}
task_id <- as.integer(task_id)

if (task_id < 1 | task_id > nrow(tissues)){
  stop(paste0("SLURM_ARRAY_TASK_ID (", task_id, ") is out of range for ",
              nrow(tissues), " tissues -- check the --array range in ",
              "jobs/deseq2.sh."))
}

tissue_file <- tissues$file[task_id]
tissue_title <- tissues$title[task_id]

message2(paste0("Processing ", tissue_title, " (task ", task_id, "/",
                nrow(tissues), ")"))

results_dir <- paste0("results/deseq2/", tissue_file, "/")
dir.create(results_dir, showWarnings = F, recursive = T)

# Load the final annotated object from 17_obj_reassembly.R ------------------
# Real raw counts come from data/06_obj_reassembly/bpcells, subset to the
# cells 17_obj_reassembly.R retained -- see header note above.

message2("Reading in metadata and raw counts")

meta <- readRDS(paste0("data/17_obj_reassembly/", tissue_file, "/metadata.rds"))

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta)]

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta, assay = "RNA")

# Join in donor demographics (sex, age at death) -----------------------------
# See header note above re: the case_number/id hyphen mismatch and the
# unverified "sex" column name.

message2("Joining donor demographics (sex, age at death)")

demo <- read.csv("tab_data/target_als_demographics_compiled.csv") %>%
  clean_names() %>%
  mutate(join_key = str_remove_all(case_number, "-")) %>%
  dplyr::select(join_key, sex, age_at_death) %>% # CONFIRM: "sex" is this file's actual column name post-clean_names()
  distinct()

obj_ids <- obj@meta.data %>%
  distinct(id) %>%
  mutate(join_key = str_remove_all(id, "-"))

unmatched <- obj_ids$id[!obj_ids$join_key %in% demo$join_key]
if (length(unmatched) > 0){
  stop(paste0("No demographics match found for donor(s): ",
              paste(unmatched, collapse = ", "),
              " -- check tab_data/target_als_demographics_compiled.csv's ",
              "case_number values against these ids before rerunning."))
}

obj@meta.data <- obj@meta.data %>%
  rownames_to_column(var = "cell") %>%
  mutate(join_key = str_remove_all(id, "-")) %>%
  left_join(demo, by = "join_key") %>%
  dplyr::select(-join_key) %>%
  column_to_rownames(var = "cell")

obj@meta.data$sex <- factor(obj@meta.data$sex)

# Factorize grouping variable, Control as the reference level ---------------

obj@meta.data$group <- factor(obj@meta.data$group,
                              levels = c("Control", "sALS", "C9orf72"))

# Extract cell type labels to establish for loop -----------------------------

cell_types <- unique(obj$cell_type3)

for (i in seq_along(cell_types)){

  message(paste0(cell_types[i], " in ", tissue_title))

  sub <- subset(obj,
                cell_type3 == cell_types[i])

  file <- str_replace_all(cell_types[i], " ", "_")

  bulk <- AggregateExpression(sub,
                              assays = "RNA",
                              return.seurat = F,
                              # layer = "counts",
                              group.by = c("orig.ident"))

  exp <- bulk$RNA

  meta_ct <- sub@meta.data %>%
    dplyr::select(c(orig.ident, group, sex, age_at_death)) %>%
    distinct() %>%
    mutate(orig.ident = str_replace_all(orig.ident, "_", "-"),
           age_scale = scale(age_at_death, center = T, scale = T)[,1])

  idx <- match(colnames(exp), meta_ct$orig.ident)
  meta_ct <- meta_ct[idx, ]
  rownames(meta_ct) <- meta_ct$orig.ident

  dds <- DESeqDataSetFromMatrix(countData = exp,
                                colData = meta_ct,
                                design = ~ sex + age_scale + group) # change this as needed

  keep <- rowSums(counts(dds) >= 10) >= 10 # change these cutoffs as needed

  dds <- dds[keep, ]

  dds <- DESeq(dds)

  saveRDS(dds,
          file = paste0(results_dir, file, "_dds.rds"))

  # resultsNames(dds)

  ct_results_dir <- paste0(results_dir, file, "/")
  dir.create(ct_results_dir, showWarnings = F, recursive = T)

  res <- results(dds,
                 contrast = c("group", "sALS", "Control"),
                 filterFun = ihw,
                 independentFiltering = T)

  res <- as.data.frame(res)

  write.csv(res,
            file = paste0(ct_results_dir, "sALS_vs_Control.csv"))

  res_shrunk <- lfcShrink(dds,
                          coef = "group_sALS_vs_Control",
                          type = "apeglm")

  write.csv(res_shrunk,
            file = paste0(ct_results_dir, "sALS_vs_Control_lfc_shrunk.csv"))

  res <- results(dds,
                 contrast = c("group", "C9orf72", "Control"),
                 filterFun = ihw,
                 independentFiltering = T)

  res <- as.data.frame(res)

  write.csv(res,
            file = paste0(ct_results_dir, "C9orf72_vs_Control.csv"))

  res_shrunk <- lfcShrink(dds,
                          coef = "group_C9orf72_vs_Control",
                          type = "apeglm")

  write.csv(res_shrunk,
            file = paste0(ct_results_dir, "C9orf72_vs_Control_lfc_shrunk.csv"))

}
