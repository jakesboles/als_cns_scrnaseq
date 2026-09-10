# Integrates consensus hdWGCNA module scores (wgcna_consensus_cns.R) with
# Milo differential neighborhood abundance results (milo.R) for microglia
# into one wide data frame: one row per neighborhood, with its UMAP
# embedding (for plotting), logFC/significance in both tissues x both
# disease-group contrasts, and each WGCNA module's mean score across the
# neighborhood's member cells. Assembly only, no plotting -- the user will
# build plots interactively from this data frame.
#
# Design notes:
# - Module score aggregation: for each neighborhood, module scores are
#   averaged across every cell miloR assigned to it (via nhoods(milo)'s
#   cell x neighborhood incidence matrix), not just its index cell's own
#   score. Confirmed with the user -- this matches Milo's own established
#   convention for overlaying a continuous per-cell covariate onto
#   neighborhoods (e.g. plotNhoodExpressionDA()), and keeps it consistent
#   with how testNhoods()'s logFC/significance are themselves computed
#   over the whole neighborhood, not one cell. milo_viz.R's index-cell-only
#   convention doesn't transfer here -- logFC/size/embedding there are
#   already neighborhood-level values computed by miloR itself, not
#   per-cell values needing aggregation.
# - Module scores come from wgcna_consensus_cns.R's
#   results/wgcna_consensus/Microglia/module_scores_ucell.csv (kNN-smoothed
#   UCell scores, one row per cell), not module_eigengenes.csv (the
#   harmonized module eigengenes) -- matching wgcna_consensus_viz.R's own
#   established use of the UCell scores as "the" module scores in this
#   project.
# - wgcna_consensus_cns.R's Microglia population comes from
#   data/18_full_integration/brain_sc (filtered to cell_type3 ==
#   "Microglia"), while milo.R's microglia population comes from
#   data/19_subclustering3/microglia -- two different pipelines, but both
#   ultimately trace back to the same 17_obj_reassembly.R per-tissue
#   cell_type3 annotations with no further cell exclusion in between, so
#   the cell barcode sets are expected to match closely. Checked at
#   runtime (not assumed) via a coverage message, with a hard stop() only
#   if coverage is implausibly low (<50%), since that would mean the
#   assumption is actually wrong rather than a handful of incidental
#   stragglers.
# - Neighborhood ID is the same integer Nhood ID testNhoods() itself uses
#   (1:ncol(nhoods(milo))) -- the natural join key across nhoods(milo)'s
#   incidence matrix, nhoodIndex(milo), and every results CSV's own "Nhood"
#   column, so everything is joined directly by Nhood rather than going
#   through milo_viz.R's index-cell-barcode indirection.
# - logFC/significance are kept raw here (not zeroed out for non-
#   significant neighborhoods the way milo_viz.R zeroes logFC for plot
#   coloring) -- this is a data table for further interactive analysis, not
#   a single fixed plot, so thresholding is left to the user. "Significant"
#   = SpatialFDR < 0.05, this project's standard threshold elsewhere
#   (DESeq2, GSEA).
# - UMAP embedding is still the index cell's own coordinates (reattached
#   from 19_subclustering3.R's harmony_umap.rds, same as milo_viz.R) --
#   that's just a plotting position for the neighborhood, not a value being
#   aggregated, so the index-cell convention is unaffected by the module-
#   score design choice above.
# - Every neighborhood is included, not just DA/significant ones -- the
#   user can filter interactively.
# - module_scores_ucell.csv is read with row.names = 1 (restoring real cell
#   barcode rownames), not read as a plain data frame with an "X" gene/cell
#   column the way this project's DESeq2 CSVs usually are -- barcodes are
#   needed here as actual rownames for the incidence-matrix join, not as a
#   data column.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(miloR)
  library(scater)
  library(Matrix)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15), collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_scrnaseq")

target_name <- "microglia"
wgcna_celltype <- "Microglia"

data_dir <- paste0("data/milo/", target_name, "/")
milo_results_dir <- paste0("results/milo/", target_name, "/")
wgcna_results_dir <- paste0("results/wgcna_consensus/", wgcna_celltype, "/")

out_dir <- paste0("results/milo_wgcna/", target_name, "/")
dir.create(out_dir, showWarnings = F, recursive = T)

# Load the Milo object and its UMAP embedding --------------------------------
# Same reattachment as milo_viz.R -- see its own header for why the saved
# Milo object has no UMAP reduction attached.

message2("Reading in Milo object and UMAP embedding")

milo <- readRDS(paste0(data_dir, "milo_obj.rds"))

umap <- readRDS(paste0("data/19_subclustering3/", target_name,
                       "/harmony_umap.rds"))
reducedDim(milo, "harmony_umap") <- Embeddings(umap)[colnames(milo), ]

nh_mat <- nhoods(milo)
if (is.null(rownames(nh_mat))) rownames(nh_mat) <- colnames(milo)

n_nhoods <- ncol(nh_mat)

# Neighborhood ID, index cell, embedding, and total size ---------------------
# One row per neighborhood, Nhood ID = column index into nhoods(milo) --
# see header note above.

message2("Assembling per-neighborhood ID, index cell, size, and embedding")

index_rows <- unlist(nhoodIndex(milo))

nhood_df <- data.frame(
  Nhood = seq_len(n_nhoods),
  cell = colnames(milo)[index_rows],
  size = Matrix::colSums(nh_mat)
)

umap_embed <- reducedDim(milo, "harmony_umap")[index_rows, ] %>%
  as.data.frame()

nhood_df <- bind_cols(nhood_df, umap_embed)

# Bring in WGCNA module scores, averaged per neighborhood ---------------------
# Mean across member cells, not index cell only -- see header note above.

message2("Reading in WGCNA module scores")

module_scores <- read.csv(paste0(wgcna_results_dir, "module_scores_ucell.csv"),
                          row.names = 1)

module_cols <- colnames(module_scores)[str_detect(colnames(module_scores),
                                                   "_UCell_kNN$")]

if (length(module_cols) == 0){
  stop("No '_UCell_kNN' columns found in ", wgcna_results_dir,
       "module_scores_ucell.csv -- check that wgcna_consensus_cns.R has ",
       "been run for ", wgcna_celltype, ".")
}

common_cells <- intersect(rownames(nh_mat), rownames(module_scores))
coverage <- length(common_cells) / nrow(nh_mat)

message2(paste0(round(coverage * 100, 1), "% of milo's ", target_name,
                " cells (", length(common_cells), " / ", nrow(nh_mat), ") ",
                "have a matching WGCNA module score"))

if (coverage < 0.5){
  stop("Fewer than 50% of milo's ", target_name, " cells could be matched ",
       "to a WGCNA module score by barcode -- the assumption that ",
       "wgcna_consensus_cns.R's ", wgcna_celltype, " population (from ",
       "data/18_full_integration/brain_sc) and milo.R's ", target_name,
       " population (from data/19_subclustering3/", target_name, ") are ",
       "the same set of cells (both traced back to the same ",
       "17_obj_reassembly.R cell_type3 annotations) appears to be wrong -- ",
       "check both scripts' source metadata before proceeding.")
}

nh_mat_common <- nh_mat[common_cells, , drop = F]
score_mat <- as.matrix(module_scores[common_cells, module_cols, drop = F])

n_scored_members <- Matrix::colSums(nh_mat_common)
score_sums <- as.matrix(Matrix::crossprod(nh_mat_common, score_mat))
mean_scores <- score_sums / n_scored_members
mean_scores[n_scored_members == 0, ] <- NA

colnames(mean_scores) <- str_remove_all(colnames(mean_scores), "_UCell_kNN$")

mean_scores_df <- as.data.frame(mean_scores) %>%
  mutate(Nhood = seq_len(n_nhoods),
         n_scored_members = as.numeric(n_scored_members)) %>%
  relocate(Nhood, n_scored_members)

nhood_df <- nhood_df %>%
  left_join(mean_scores_df, by = "Nhood")

# Bring in per-tissue, per-contrast logFC/significance ------------------------
# Raw (not zeroed) -- see header note above. Short tissue codes
# ("brain"/"sc") match deseq_viz2.R's naming convention.

message2("Reading in Milo differential abundance results")

combos <- tribble(
  ~tissue_file,           ~tissue_short, ~contrast,            ~contrast_short,
  "Motor_cortex",         "brain",       "sALS_vs_Control",    "sALS",
  "Motor_cortex",         "brain",       "C9orf72_vs_Control", "C9",
  "Cervical_spinal_cord", "sc",          "sALS_vs_Control",    "sALS",
  "Cervical_spinal_cord", "sc",          "C9orf72_vs_Control", "C9"
)

for (i in seq_len(nrow(combos))){

  results_csv <- paste0(milo_results_dir, combos$contrast[i], "_",
                        combos$tissue_file[i], "_nhood_results.csv")

  logfc_col <- paste0("logFC_", combos$tissue_short[i], "_",
                      combos$contrast_short[i])
  sig_col <- paste0("sig_", combos$tissue_short[i], "_",
                    combos$contrast_short[i])

  if (!file.exists(results_csv)){
    message2(paste0("Missing ", results_csv, " -- filling ", logfc_col, "/",
                    sig_col, " with NA (likely skipped by milo.R's ",
                    "min_donors_per_group check)"))
    nhood_df[[logfc_col]] <- NA_real_
    nhood_df[[sig_col]] <- NA
    next
  }

  res <- read.csv(results_csv) %>%
    dplyr::select(Nhood, logFC, SpatialFDR) %>%
    mutate(!!sig_col := !is.na(SpatialFDR) & SpatialFDR < 0.05) %>%
    dplyr::rename(!!logfc_col := logFC) %>%
    dplyr::select(-SpatialFDR)

  nhood_df <- nhood_df %>%
    left_join(res, by = "Nhood")
}

# Save -------------------------------------------------------------------

message2("Saving joint Milo/WGCNA data frame")

write.csv(nhood_df,
          file = paste0(out_dir, "milo_wgcna_joint_data.csv"),
          row.names = F)
