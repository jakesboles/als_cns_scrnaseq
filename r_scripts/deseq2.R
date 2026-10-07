# Pseudobulk DESeq2 differential expression within each cell type
# (cell_type3) of one tissue, comparing sALS vs. Control and C9orf72 vs.
# Control, with per-cell-type-per-sample abundance filtering to avoid
# DESeq2 failures on sparsely represented cell types, and age/sex as
# covariates in the model. Runs as a SLURM job array (see jobs/deseq2.sh),
# one task per tissue, loading each tissue's final annotated object from
# 17_obj_reassembly.R.

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

# A pseudobulk sample built from too few cells is mostly zero, which can
# make every gene contain a zero in some sample -- DESeq2's default
# median-of-ratios size factor estimation then fails outright
# ("every gene contains at least one zero, cannot compute log geometric
# means"). Below this per-sample cell count, that sample is dropped from
# the cell type's pseudobulk; if too few samples remain in any group
# after dropping, the whole cell type is skipped rather than run on an
# unreliable/unbalanced design.
min_cells_per_sample <- 10 # change as needed
min_samples_per_group <- 3 # change as needed

for (i in seq_along(cell_types)){

  message(paste0(cell_types[i], " in ", tissue_title))

  sub <- subset(obj,
                cell_type3 == cell_types[i])

  file <- str_replace_all(cell_types[i], " ", "_")

  ct_results_dir <- paste0(results_dir, file, "/")
  dir.create(ct_results_dir, showWarnings = F, recursive = T)

  bulk <- AggregateExpression(sub,
                              assays = "RNA",
                              return.seurat = F,
                              # layer = "counts",
                              group.by = c("orig.ident"))

  exp <- bulk$RNA

  cell_counts <- sub@meta.data %>%
    dplyr::count(orig.ident, name = "n_cells")

  # Every sample present for this cell type, whether or not it survives
  # the min_cells_per_sample filter below -- saved so a skipped/thinned
  # cell type's sample composition can be checked later without rerunning
  # anything.
  sample_table <- sub@meta.data %>%
    dplyr::select(orig.ident, group) %>%
    distinct() %>%
    left_join(cell_counts, by = "orig.ident") %>%
    mutate(retained = n_cells >= min_cells_per_sample) %>%
    arrange(group, orig.ident)

  write.csv(sample_table,
            file = paste0(ct_results_dir, "sample_filtering.csv"),
            row.names = F)

  meta_ct <- sub@meta.data %>%
    dplyr::select(c(orig.ident, group, sex, age_at_death)) %>%
    distinct() %>%
    left_join(cell_counts, by = "orig.ident") %>%
    filter(n_cells >= min_cells_per_sample)

  if (any(table(meta_ct$group) < min_samples_per_group)){
    message(paste0("Skipping ", cell_types[i], " (", tissue_title, ") -- ",
                   "fewer than ", min_samples_per_group, " samples per ",
                   "group have >= ", min_cells_per_sample, " ",
                   cell_types[i], " cells."))
    next
  }

  meta_ct <- meta_ct %>%
    dplyr::select(-n_cells) %>%
    mutate(orig.ident = str_replace_all(orig.ident, "_", "-"),
           age_scale = scale(age_at_death, center = T, scale = T)[,1])

  exp <- exp[, meta_ct$orig.ident, drop = F]

  idx <- match(colnames(exp), meta_ct$orig.ident)
  meta_ct <- meta_ct[idx, ]
  rownames(meta_ct) <- meta_ct$orig.ident

  # The abundance filter above catches the most common cause of DESeq2's
  # "every gene contains at least one zero" size factor error, but not
  # every case (e.g. a gene that's zero in every retained sample even
  # though each sample individually cleared min_cells_per_sample). Wrap
  # the whole DESeq2 pipeline so a failure on one cell type is logged and
  # skipped instead of killing the rest of this tissue's task.
  tryCatch({

    dds <- DESeqDataSetFromMatrix(countData = exp,
                                  colData = meta_ct,
                                  design = ~ sex + age_scale + group) # change this as needed

    keep <- rowSums(counts(dds) >= 10) >= 10 # change these cutoffs as needed

    dds <- dds[keep, ]

    dds <- DESeq(dds)

    saveRDS(dds,
            file = paste0(results_dir, file, "_dds.rds"))

    # resultsNames(dds)

    res <- results(dds,
                   contrast = c("group", "sALS", "Control"),
                   filterFun = ihw,
                   independentFiltering = T)

    res <- as.data.frame(res)

    write.csv(res,
              file = paste0(ct_results_dir, "sALS_vs_Control.csv"))

    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "group_sALS_vs_Control",
                              type = "apeglm")
    })

    write.csv(res_shrunk,
              file = paste0(ct_results_dir, "sALS_vs_Control_lfc_shrunk.csv"))

    res <- results(dds,
                   contrast = c("group", "C9orf72", "Control"),
                   filterFun = ihw,
                   independentFiltering = T)

    res <- as.data.frame(res)

    write.csv(res,
              file = paste0(ct_results_dir, "C9orf72_vs_Control.csv"))

    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "group_C9orf72_vs_Control",
                              type = "apeglm")
    })

    write.csv(res_shrunk,
              file = paste0(ct_results_dir, "C9orf72_vs_Control_lfc_shrunk.csv"))

  }, error = function(e){
    message(paste0("Skipping ", cell_types[i], " (", tissue_title, ") -- ",
                   "DESeq2 pipeline failed: ", conditionMessage(e)))
  })

}
