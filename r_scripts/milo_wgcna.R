# Integrates consensus hdWGCNA module scores (wgcna_consensus_cns.R) with
# Milo differential neighborhood abundance results (milo.R) for microglia
# into one long-format data frame: one row per (neighborhood, comparison)
# pair, with the neighborhood's UMAP embedding (for plotting), that
# comparison's logFC/significance, and each WGCNA module's mean score
# across just that comparison's member cells within the neighborhood.
# Assembly only, no plotting -- the user will build plots interactively
# from this data frame.
#
# Design notes:
# - Reworked from an earlier wide-format version of this script (one row
#   per neighborhood, module scores averaged across *all* of a
#   neighborhood's member cells regardless of tissue/group) after the user
#   caught a real bug in that design: a neighborhood's cell membership
#   doesn't change between comparisons (only abundance does), but pooling
#   every group's cells together to score a neighborhood meant every
#   comparison saw the exact same module score for a given neighborhood --
#   the "blue" module histograms for brain_C9/brain_sALS/sc_C9/sc_sALS came
#   out identically distributed, which is what surfaced the bug.
# - The fix: each neighborhood is now scored separately per comparison,
#   using only the member cells that actually belong to that comparison --
#   i.e. cells in the relevant tissue AND in either "Control" or that
#   comparison's disease group (excluding the third, irrelevant group).
#   Those per-comparison cell subsets genuinely differ (a neighborhood's
#   brain_C9 subset and brain_sALS subset share only their Control cells,
#   not their disease-group cells), so the resulting module scores can
#   legitimately differ across comparisons -- unlike logFC/size/embedding,
#   which are single neighborhood-level values miloR itself computes once
#   and that don't need this per-comparison treatment.
# - Module score aggregation within a comparison's cell subset is still a
#   mean across member cells (via nhoods(milo)'s cell x neighborhood
#   incidence matrix), not just the index cell's score -- confirmed with
#   the user previously, matches Milo's own convention for overlaying a
#   continuous per-cell covariate onto neighborhoods (e.g.
#   plotNhoodExpressionDA()).
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
#   runtime (not assumed) via a coverage message (both overall and per
#   comparison), with a hard stop() only if overall coverage is implausibly
#   low (<50%), since that would mean the assumption is actually wrong
#   rather than a handful of incidental stragglers.
# - `tissue`/`group` come straight from colData(milo) -- as.SingleCellExperiment()
#   (called in milo.R when building the Milo object) carries the Seurat
#   object's metadata columns over as colData, so these are the same
#   `tissue`/`group` columns milo.R itself filters/models on, not reloaded
#   from anywhere else.
# - Neighborhood ID is the same integer Nhood ID testNhoods() itself uses
#   (1:ncol(nhoods(milo))) -- the natural join key across nhoods(milo)'s
#   incidence matrix, nhoodIndex(milo), and every results CSV's own "Nhood"
#   column.
# - logFC/significance are kept raw (not zeroed out for non-significant
#   neighborhoods the way milo_viz.R zeroes logFC for plot coloring) --
#   this is a data table for further interactive analysis, not a single
#   fixed plot, so thresholding is left to the user. "Significant" =
#   SpatialFDR < 0.05, this project's standard threshold elsewhere (DESeq2,
#   GSEA).
# - UMAP embedding is still the index cell's own coordinates (reattached
#   from 19_subclustering3.R's harmony_umap.rds, same as milo_viz.R) and
#   doesn't vary by comparison -- it's a plotting position for the
#   neighborhood, not a value being aggregated, so it's joined back onto
#   every comparison's row for that neighborhood unchanged.
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
# One row per neighborhood -- comparison-independent, joined onto every
# comparison's rows below. Nhood ID = column index into nhoods(milo).

message2("Assembling per-neighborhood ID, index cell, size, and embedding")

index_rows <- unlist(nhoodIndex(milo))

nhood_meta <- data.frame(
  Nhood = seq_len(n_nhoods),
  cell = colnames(milo)[index_rows],
  size = Matrix::colSums(nh_mat)
)

umap_embed <- reducedDim(milo, "harmony_umap")[index_rows, ] %>%
  as.data.frame()

nhood_meta <- bind_cols(nhood_meta, umap_embed)

# Read WGCNA module scores ----------------------------------------------
# See header note above re: row.names = 1.

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
                "have a matching WGCNA module score overall"))

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

# Per-comparison neighborhood scoring ------------------------------------
# Each comparison gets its own cell subset (this tissue AND (Control OR
# that comparison's disease group)) -- see header note above for why this
# replaces the earlier "average over the whole neighborhood" design.

message2("Computing per-comparison module scores and DA results")

combos <- tribble(
  ~tissue_title,          ~tissue_short, ~contrast,            ~contrast_short, ~disease_group,
  "Motor cortex",         "brain",       "sALS_vs_Control",    "sALS",          "sALS",
  "Motor cortex",         "brain",       "C9orf72_vs_Control", "C9",            "C9orf72",
  "Cervical spinal cord", "sc",          "sALS_vs_Control",    "sALS",          "sALS",
  "Cervical spinal cord", "sc",          "C9orf72_vs_Control", "C9",            "C9orf72"
) %>%
  mutate(tissue_file = str_replace_all(tissue_title, " ", "_"),
         comparison = paste0(tissue_short, "_", contrast_short))

milo_tissue <- colData(milo)$tissue
milo_group <- colData(milo)$group

comparison_rows <- list()

for (i in seq_len(nrow(combos))){

  comparison <- combos$comparison[i]

  message2(paste0("Comparison: ", comparison))

  group_mask <- milo_tissue == combos$tissue_title[i] &
    milo_group %in% c(combos$disease_group[i])

  comparison_cells <- colnames(milo)[group_mask]
  scored_cells <- intersect(comparison_cells, common_cells)

  coverage_i <- length(scored_cells) / length(comparison_cells)
  message2(paste0("  ", length(comparison_cells), " cells in this ",
                  "tissue/group subset, ", length(scored_cells), " (",
                  round(coverage_i * 100, 1), "%) have a WGCNA module ",
                  "score"))

  # Module scores: mean over just this comparison's WGCNA-scored member
  # cells per neighborhood.
  nh_mat_scored <- nh_mat[scored_cells, , drop = F]
  score_mat_i <- as.matrix(module_scores[scored_cells, module_cols, drop = F])

  n_scored_members <- Matrix::colSums(nh_mat_scored)
  score_sums <- as.matrix(Matrix::crossprod(nh_mat_scored, score_mat_i))
  mean_scores <- score_sums / n_scored_members
  mean_scores[n_scored_members == 0, ] <- NA
  colnames(mean_scores) <- str_remove_all(colnames(mean_scores), "_UCell_kNN$")

  # Total member cells in this comparison's subset per neighborhood,
  # regardless of WGCNA-score availability -- for transparency alongside
  # n_scored_members.
  nh_mat_group <- nh_mat[comparison_cells, , drop = F]
  n_group_cells <- Matrix::colSums(nh_mat_group)

  scores_df_i <- as.data.frame(mean_scores) %>%
    mutate(Nhood = seq_len(n_nhoods),
           comparison = comparison,
           n_group_cells = as.numeric(n_group_cells),
           n_scored_members = as.numeric(n_scored_members)) %>%
    relocate(Nhood, comparison, n_group_cells, n_scored_members)

  results_csv <- paste0(milo_results_dir, combos$contrast[i], "_",
                        combos$tissue_file[i], "_nhood_results.csv")

  if (!file.exists(results_csv)){
    message2(paste0("  Missing ", results_csv, " -- filling logFC/sig with ",
                    "NA (likely skipped by milo.R's min_donors_per_group ",
                    "check)"))
    scores_df_i$logFC <- NA_real_
    scores_df_i$sig <- NA
  } else {
    res <- read.csv(results_csv) %>%
      dplyr::select(Nhood, logFC, SpatialFDR) %>%
      mutate(sig = !is.na(SpatialFDR) & SpatialFDR < 0.05) %>%
      dplyr::select(-SpatialFDR)

    scores_df_i <- scores_df_i %>%
      left_join(res, by = "Nhood")
  }

  comparison_rows[[comparison]] <- scores_df_i
}

comparison_df <- list_rbind(comparison_rows)

comparison_df %>% 
  # mutate(z = sin(turquoise) * cos(logFC)) %>%
  filter(sig == TRUE) %>%
  na.omit() %>%
  ggplot(aes(x = turquoise,
             y = logFC)) + 
  geom_density_2d() + 
  facet_wrap(. ~ comparison)

comparison_df %>% 
  filter(sig == TRUE & 
           logFC > 0) %>% 
  na.omit() %>%
  ggplot(aes(x = logFC,
             y = turquoise)) + 
  geom_point(aes(color = comparison)) +
  # geom_density_2d(aes(color = comparison)) +
  scale_color_manual(values = JCO_Four()) + 
  facet_wrap(. ~ comparison) + 
  theme_linedraw()

# Join comparison-independent neighborhood metadata onto every comparison's
# rows and save -----------------------------------------------------------

message2("Saving joint Milo/WGCNA data frame")

nhood_df <- comparison_df %>%
  left_join(nhood_meta, by = "Nhood") %>%
  relocate(Nhood, comparison, cell, size, harmonyumap_1, harmonyumap_2,
           n_group_cells, n_scored_members, logFC, sig)

write.csv(nhood_df,
          file = paste0(out_dir, "milo_wgcna_joint_data.csv"),
          row.names = F)
