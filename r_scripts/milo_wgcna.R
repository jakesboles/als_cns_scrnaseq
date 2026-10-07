# Integrates consensus hdWGCNA module scores (wgcna_consensus_cns.R) with
# Milo differential neighborhood abundance results (milo.R) for microglia
# into one long-format data frame -- one row per (neighborhood,
# comparison) pair, each module's mean score computed from only that
# comparison's own tissue/disease-group member cells (never Control, and
# never pooled across comparisons) -- plus a few exploratory plots of
# module score vs. neighborhood logFC. Interactive script, not a SLURM
# job.

suppressMessages({
  library(Seurat)
  library(tidyverse)
  library(miloR)
  library(scater)
  library(Matrix)
  library(ggExtra)
  library(patchwork)
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
# Each comparison gets its own cell subset (this tissue AND that
# comparison's disease group only, Control excluded) -- see header note
# above for why this replaces the earlier "average over the whole
# neighborhood" (and then "average over the whole neighborhood minus the
# other disease group") designs.

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

  # Total member cells from this tissue *regardless of group* (Control
  # included) -- diagnostic only, not used in any score. A neighborhood
  # with n_group_cells == 0 but a healthy n_tissue_cells just has no
  # cells of this specific disease group locally; one where
  # n_tissue_cells is *also* ~0 barely exists in this tissue at all (it's
  # dominated by the other tissue's cells), which is the more likely
  # explanation for the exact-zero logFC/untested pattern the user found
  # in some of these rows -- testNhoods()'s per-tissue GLM fit degenerates
  # for a neighborhood with essentially no counts across that tissue's
  # samples.
  tissue_mask <- milo_tissue == combos$tissue_title[i]
  nh_mat_tissue <- nh_mat[colnames(milo)[tissue_mask], , drop = F]
  n_tissue_cells <- Matrix::colSums(nh_mat_tissue)

  scores_df_i <- as.data.frame(mean_scores) %>%
    mutate(Nhood = seq_len(n_nhoods),
           comparison = comparison,
           n_tissue_cells = as.numeric(n_tissue_cells),
           n_group_cells = as.numeric(n_group_cells),
           n_scored_members = as.numeric(n_scored_members)) %>%
    relocate(Nhood, comparison, n_tissue_cells, n_group_cells, n_scored_members)

  results_csv <- paste0(milo_results_dir, combos$contrast[i], "_",
                        combos$tissue_file[i], "_nhood_results.csv")

  if (!file.exists(results_csv)){
    message2(paste0("  Missing ", results_csv, " -- filling logFC/sig with ",
                    "NA (likely skipped by milo.R's min_donors_per_group ",
                    "check)"))
    scores_df_i$logFC <- NA_real_
    scores_df_i$sig <- NA
  } else {
    # sig stays NA (not FALSE) when SpatialFDR itself is NA -- keeps
    # "genuinely tested and not significant" distinguishable from
    # "untested/degenerate fit for this neighborhood in this tissue" (see
    # n_tissue_cells note above), rather than collapsing both into FALSE.
    res <- read.csv(results_csv) %>%
      dplyr::select(Nhood, logFC, SpatialFDR) %>%
      mutate(sig = if_else(is.na(SpatialFDR), NA, SpatialFDR < 0.05)) %>%
      dplyr::select(-SpatialFDR)

    scores_df_i <- scores_df_i %>%
      left_join(res, by = "Nhood")
  }

  comparison_rows[[comparison]] <- scores_df_i
}

comparison_df <- list_rbind(comparison_rows)

# comparison_df %>% 
#   filter(sig == TRUE) %>% 
#   na.omit() %>% 
#   group_by(comparison) %>% 
#   summarize(n = n())
# 
# comparison_df %>%
#   filter(sig == TRUE) %>%
#   na.omit() %>%
#   group_by(comparison) %>%
#   summarise(
#     n = n(),
#     turquoise_min = min(turquoise), turquoise_max = max(turquoise), turquoise_sd = sd(turquoise),
#     logFC_min = min(logFC), logFC_max = max(logFC), logFC_sd = sd(logFC),
#     any_nonfinite = any(!is.finite(turquoise) | !is.finite(logFC))
#   )
# 
# comparison_df %>%
#   filter(sig == TRUE) %>%
#   na.omit() %>%
#   group_by(comparison) %>%
#   summarise(
#     n = n(),
#     n_distinct_turquoise = n_distinct(turquoise),
#     n_distinct_logFC = n_distinct(logFC),
#     turquoise_iqr = IQR(turquoise),
#     logFC_iqr = IQR(logFC),
#     turquoise_bw = MASS::bandwidth.nrd(turquoise),
#     logFC_bw = MASS::bandwidth.nrd(logFC)
#   )

module <- "blue"

p1 <- comparison_df %>%
  filter(sig == TRUE & logFC > 0 & 
           str_detect(comparison, "C9")) %>%
  mutate(tissue = str_split_i(comparison, "_", i = 1) %>% 
           factor(levels = c("brain", "sc"),
                  labels = c("Motor cortex", "Cervical spinal cord"))) %>%
  # mutate(group = case_when(str_detect(comparison, "C9") ~ "C9orf72-ALS",
  #                          str_detect(comparison, "sALS") ~ "sALS") %>% 
  #          factor(levels = c("sALS", "C9orf72-ALS"))) %>%
  na.omit() %>%
  ggplot(aes(x = !!sym(module), y = logFC)) +
  # geom_density_2d(aes(color = comparison), contour_var = "ndensity") +
  geom_point(aes(color = tissue), alpha = 0.7, size = 3) +
  scale_color_manual(values = c("#EFC000", "#0073C2")) +
  labs(y = "Nhood log2FC",
       x = paste0(str_to_title(module), " module score")) +
  # facet_wrap(. ~ comparison, scales = "fixed") +
  # scale_color_viridis_c() +
  ggtitle(paste0(str_to_title(module), " expression vs\nNhood fold change in C9orf72-ALS")) +
  theme_linedraw() + 
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))

p1h <- ggMarginal(p1, type = "density",
           groupFill = T)
ggsave(p1h,
       filename = paste0(out_dir, module, "_vs_milo_fc.png"),
       units = "in", dpi = 600,
       height = 5, width = 5)

# find a good way to statistically analyze the effect of comparison on the 
# relationship between logFC and module score

# Join comparison-independent neighborhood metadata onto every comparison's
# rows and save -----------------------------------------------------------

message2("Saving joint Milo/WGCNA data frame")

nhood_df <- comparison_df %>%
  left_join(nhood_meta, by = "Nhood") %>%
  relocate(Nhood, comparison, cell, size, harmonyumap_1, harmonyumap_2,
           n_tissue_cells, n_group_cells, n_scored_members, logFC, sig)

write.csv(nhood_df,
          file = paste0(out_dir, "milo_wgcna_joint_data.csv"),
          row.names = F)
