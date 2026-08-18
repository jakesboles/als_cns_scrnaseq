suppressMessages({
  library(tidyverse)
  library(ggplot2)
  library(paletteer)
  library(ggbeeswarm)
  library(scales)
  library(patchwork)
})

setwd("/projects/b1169/boles/als_cns_scrnaseq")

plots_dir <- "plots/05_qc3/"
dir.create(plots_dir, 
           showWarnings = F,
           recursive = T)

# in_dir <- "data/04_doubletfinder/metadata_persample/"
# 
# files <- list.files(in_dir)
# 
# for (i in seq_along(files)){
#   
#   message(paste0(i, " out of ", length(files)))
#   
#   if (i == 1) { 
#     meta <- readRDS(paste0(in_dir, files[i]))
#     
#     col_idx <- str_detect(colnames(meta), "pANN", negate = T)
#     
#     meta <- meta[, col_idx]
#   } else { 
#     meta2 <- readRDS(paste0(in_dir, files[i]))
#     
#     col_idx <- str_detect(colnames(meta2), "pANN", negate = T)
#     
#     meta2 <- meta2[, col_idx]
#     
#     meta <- rbind(meta, meta2)
#   }
# }

meta <- read.csv("tab_data/metadata.csv")


# Get pre-filter median stats & contamination ---------------------------------------------

pre_stats_list <- list()

pre_stats_files <- list.files("tab_data/02_qc1",
                              full.names = T)
stats_tissues <- c("Motor cortex", "Muscle", "Cervical spinal cord")

for (i in seq_along(pre_stats_files)){
  pre_stats_list[[i]] <- read.csv(pre_stats_files[i])[-31, ]
  
  pre_stats_list[[i]]$Tissue <- stats_tissues[i]
}

pre_stats <- list_rbind(pre_stats_list)

pre_stats <- pre_stats %>%
  dplyr::rename("ID" = "id") %>%
  # mutate(ID = str_replace_all(ID, "GWF-", "GWF")) %>%
  left_join(meta, by = c("ID", "Tissue")) %>%
  mutate(Batch = factor(Batch,
                        levels = c(1:6)),
         Group = factor(Group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9ORF72-ALS")),
         Tissue = factor(Tissue,
                         levels = c("Motor cortex", "Cervical spinal cord", "Muscle"),
                         labels = c("Motor cortex", "Cervical\nspinal cord", "Skeletal muscle")))

contam <- read.csv("tab_data/00_cellbender/cellbender_metrics.csv")[1:90, ] %>%
  dplyr::select(-c(pool, group)) %>% 
  # separate(id, into = c("ID", "Tissue"),
  #          sep = "_") %>%
  mutate(tissue = factor(tissue, 
                         levels = c("Motor cortex", "Cervical spinal cord", "Muscle"),
                         labels = c("Motor cortex", "Cervical\nspinal cord", "Skeletal muscle"))) %>%
    dplyr::rename("Tissue" = "tissue",
                  "ID" = "sample")

pre_stats <- pre_stats %>%
  left_join(contam,
            by = c("Tissue", "ID"))
# Median QC of counts, features, and mito before filtering ----------------


plot <- function(stats_df, grouping_var, dep_var){
  
  if (deparse(substitute(grouping_var)) == "Batch") {
    fill_pal <- c("dodgerblue1", "springgreen1", "gold", "orangered", "violetred", "slateblue1")
    color_pal <- c("dodgerblue4", "springgreen4", "darkgoldenrod", "orangered4", "violetred4", "slateblue4")
    xlab <- "Batch"
  } else {
    if (deparse(substitute(grouping_var)) == "Tissue") {
      fill_pal <- paletteer_d("ggthemes::few_Light")
      color_pal <- paletteer_d("ggthemes::few_Dark")
      xlab <- "Tissue"
    } else {
      if (deparse(substitute(grouping_var)) == "Group") {
        fill_pal <- c("orchid1", "chartreuse", "firebrick1")
        color_pal <- c("orchid4", "chartreuse4", "firebrick4")
        xlab <- "Group"
      } else {
        if (deparse(substitute(grouping_var)) == "ID") {
          fill_pal <- DiscretePalette_scCustomize(30, palette = "polychrome")
          color_pal <- DiscretePalette_scCustomize(30, palette = "polychrome")
          xlab <- "Patient identifier"
        }
      }
    }
  }
  
  if (deparse(substitute(dep_var)) == "Median_percent_mito") {
    ylab <- "% mitochondrial\ngene counts"
    title <- "Median mitochondrial gene content\n(pre-filtering)"
  } else {
    if (deparse(substitute(dep_var)) == "Median_nCount_RNA") {
      ylab <- "# of unique UMIs"
      title <- "Median number of UMIs per cell\n(pre-filtering)"
    } else {
      if (deparse(substitute(dep_var)) == "Median_nFeature_RNA") {
        ylab <- "# of unique genes"
        title <- "Median number of genes per cell\n(pre-filtering)"
      } else {
        if (deparse(substitute(dep_var)) == "Median_log10GenesPerUMI") {
          ylab <- "log10(Genes) / log10(UMIs)"
          title <- "Median cellular complexity\n(pre-filtering)"
        } else {
          if (deparse(substitute(dep_var)) == "retained_percent") {
            ylab <- "% of total cells"
            title <- "Percentage of cells kept\nafter QC"
          } else {
            if (deparse(substitute(dep_var)) == "retained_count") {
              ylab <- "Count"
              title <- "Cell counts after QC"
            } else {
              if(deparse(substitute(dep_var)) == "fraction_counts_removed_from_cells") {
                ylab <- "CellBender-inferred\ncontamination fraction"
                title <- "Estimated ambient RNA\ncontamination"
              }
            }
          }
        }
      }
    }
  }
  
  full_title <- paste0(title, " - by ", str_to_lower(xlab))
  
  # return(show_col(fill_pal))


  stats_df %>%
    ggplot(aes(x = {{grouping_var}},
               y = {{dep_var}})) +
    geom_quasirandom(shape = 21,
                     aes(fill = {{grouping_var}}),
                     show.legend = F,
                     size = 4,
                     alpha = 0.6) +
    stat_summary(aes(color = {{grouping_var}}),
                 fun = "mean",
                 geom = "crossbar",
                 width = 0.6,
                 show.legend = F) +
    stat_summary(aes(color = {{grouping_var}}),
                 fun.data = "mean_se",
                 geom = "errorbar",
                 width = 0.3,
                 show.legend = F,
                 linewidth = 1.2) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = color_pal) +
    labs(y = ylab,
         x = xlab) + 
    ggtitle(full_title) +
    theme_bw(base_size = 10) +
    theme(axis.text = element_text(color = "black"),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12))
}  

design <- "
AAABBB
CCCCCC
"

plot(pre_stats, Tissue, Median_percent_mito) + 
  plot(pre_stats, Group, Median_percent_mito) +
  plot(pre_stats, Batch, Median_percent_mito) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "prefilter_mito.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

plot(pre_stats, Tissue, Median_nCount_RNA) + 
  plot(pre_stats, Group, Median_nCount_RNA) +
  plot(pre_stats, Batch, Median_nCount_RNA) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "prefilter_numi.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

plot(pre_stats, Tissue, Median_nFeature_RNA) + 
  plot(pre_stats, Group, Median_nFeature_RNA) +
  plot(pre_stats, Batch, Median_nFeature_RNA) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "prefilter_nfeature.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

plot(pre_stats, Tissue, Median_log10GenesPerUMI) + 
  plot(pre_stats, Group, Median_log10GenesPerUMI) +
  plot(pre_stats, Batch, Median_log10GenesPerUMI) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "prefilter_complexity.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

plot(pre_stats, Tissue, fraction_counts_removed_from_cells) + 
  plot(pre_stats, Group, fraction_counts_removed_from_cells) +
  plot(pre_stats, Batch, fraction_counts_removed_from_cells) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "prefilter_contamination.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

# Get post-filter median stats --------------------------------------------

post_stats_list <- list()

post_stats_files <- list.files("tab_data/03_qc2/",
                               full.names = T)
post_stats_files <- post_stats_files[str_detect(post_stats_files, "median")]

for (i in seq_along(post_stats_files)){
  post_stats_list[[i]] <- read.csv(post_stats_files[i])[-31, ]
  
  post_stats_list[[i]]$Tissue <- stats_tissues[i]
}

post_stats <- list_rbind(post_stats_list)

post_stats <- post_stats %>%
  dplyr::rename("ID" = "id") %>%
  left_join(meta, by = c("ID", "Tissue")) %>%
  mutate(Batch = factor(Batch,
                        levels = c(1:6)),
         Group = factor(Group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9ORF72-ALS")),
         Tissue = factor(Tissue,
                         levels = c("Motor cortex", "Cervical spinal cord", "Muscle"),
                         labels = c("Motor cortex", "Cervical\nspinal cord", "Skeletal muscle")))

counts <- read.csv("tab_data/03_qc2/filtered_counts.csv") %>%
  dplyr::rename("ID" = "id", "Tissue" = "tissue") %>%
  mutate(Tissue = factor(Tissue,
                         levels = c("Motor cortex", "Cervical spinal cord", "Skeletal muscle"),
                         labels = c("Motor cortex", "Cervical\nspinal cord", "Skeletal muscle")))

post_stats <- post_stats %>%
  left_join(counts, by = c("ID", "Tissue")) %>%
  dplyr::select(-Cell_count)

# Median cell counts & percentage loss post-filtering ---------------------

plot(post_stats, Tissue, retained_count) + 
  plot(post_stats, Group, retained_count) +
  plot(post_stats, Batch, retained_count) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "postfilter_cell_count.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

plot(post_stats, Tissue, retained_percent) + 
  plot(post_stats, Group, retained_percent) +
  plot(post_stats, Batch, retained_percent) + 
  plot_layout(design = design)
ggsave(paste0(plots_dir, "postfilter_cell_percent.png"),
       units = "in", dpi = 600,
       height = 6, width = 8)

