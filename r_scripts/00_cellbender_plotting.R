library(tidyverse)
library(scCustomize)
library(ggplot2)
library(paletteer)

# Until cells begin getting called and objects start getting integrated,
# I will keep processing all 3 tissues together.
# It will be easier to start with a procured muscle object than have to redo 
# all this later

setwd("/projects/b1169/boles/als_cns_scrnaseq")

cellbender_outs <- "/projects/b1042/Gate_Lab/boles/als_multitissue/CellbenderOutput/" 

results_dir <- "results/00_cellbender/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

outs <- list()

cellbender_dirs <- list.dirs(cellbender_outs, 
                             full.names = F,
                             recursive = F)[2:7]

for (i in seq_along(cellbender_dirs)){
  files <- list.files(paste0(cellbender_outs, cellbender_dirs[i]),
                      full.names = T,
                      recursive = F)[str_detect(list.files(paste0(cellbender_outs, cellbender_dirs[i])), "_metrics.csv")]
  
  metrics <- map(files, read.csv,
                 header = F,
                 row.names = 1)
  
  metrics <- list_cbind(metrics)
  
  # metrics <- t(metrics)
  
  samples <- str_split_i(files, "[/]", i = 9)
  samples <- str_remove_all(samples, "_metrics.csv")
  
  colnames(metrics) <- samples
  
  outs[[i]] <- metrics
  
}

metrics <- list_cbind(outs)

metrics <- t(metrics) %>%
  as.data.frame() 

metrics <- metrics %>%
  rownames_to_column(var = "sample")

idx <- str_detect(metrics$sample, "cellbender",
                  negate = T)

metrics <- metrics[idx, ]

meta <- read.csv("tab_data/metadata.csv")
colnames(meta) <- c("sample", "tissue", "pool", "dummy_code", "barcode", "index", "group")

meta <- meta %>%
  dplyr::select(c(sample, tissue, pool, group))

metrics <- metrics %>%
  mutate(tissue = case_when(str_detect(sample, "_b") ~ "Motor cortex",
                          str_detect(sample, "_m") ~ "Muscle",
                          str_detect(sample, "_s") ~ "Cervical spinal cord",
                          .default = "Hypothalamus")) %>%
  mutate(sample = gsub("\\_.*", "", sample)) %>%
  mutate(sample = str_replace_all(sample, "GWF-", "GWF")) %>%
  mutate(sample = str_replace_all(sample, "GBB-19-47", "GWF19-47") %>%
           str_replace_all("GBB-20-54", "GWF20-54") %>%
           str_replace_all("GBB-21-56", "GWF21-56")) %>%
  left_join(meta, by = c("sample", "tissue")) %>%
  mutate(pool = if_else(str_detect(sample, "UWA"), 6, pool))


metrics <- metrics %>%
  mutate(pool = factor(pool,
                        levels = c(1:6),
                        labels = paste("Pool ", c(1:6))),
         tissue = factor(tissue,
                         levels = c("Motor cortex", "Cervical spinal cord",
                                    "Muscle", "Hypothalamus")))

write.csv(metrics,
          file = paste0(results_dir, "cellbender_metrics.csv"),
          row.names = F,
          quote = F)

# Left off here, can make plots later for publication with the above .csv
  

metrics %>% 
  ggplot(aes(x = sample, y = average_counts_removed_per_cell)) + 
  geom_col(aes(fill = tissue)) +
  scale_fill_paletteer_d(palette = "ggsci::default_jama") +
  facet_wrap(. ~ pool,
             ncol = 6,
             scales = "free_x") +
  labs(y = "Contamination (%)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 90,
                               hjust = 1, vjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_text(face = "bold",
                              size = 8),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    axis.text = element_text(color = "black")
  )
