library(Seurat)
library(tidyverse)

obj_dir <- "/projects/b1169/boles/als_multitissue_scfrp/data/03_qc2/"

meta_dir <- "/projects/b1169/boles/als_multitissue_scfrp/data/04_doubletfinder/"
meta_files <- list.files(meta_dir)
meta_list <- list()

data_out_dir <- "/projects/b1169/boles/als_multitissue_scfrp/data/06_obj_assembly/"
dir.create(data_out_dir,
           recursive = T,
           showWarnings = F)

# Read in object ----------------------------------------------------------

obj <- list.files(obj_dir)
obj <- readRDS(paste0(obj_dir, obj))

rownames(obj) <- str_replace_all(rownames(obj), "GWF-", "GWF")

# Read in metadata and collapse into one dataframe ------------------------

for (i in seq_along(meta_files)){
  meta_list[[i]] <- readRDS(paste0(meta_dir, meta_files[i]))
  
  meta_list[[i]] <- meta_list[[i]]@meta.data %>%
    dplyr::select(c(orig.ident:Phase, DF.unadj:DF.adj))
}

meta <- list_rbind(meta_list)

setdiff(colnames(obj), rownames(meta))

obj <- AddMetaData(obj, meta)

obj@meta.data %>% str()

names(which(sapply(obj@meta.data, anyNA)))

obj@meta.data %>%
  filter(is.na(Batch)) %>%
  pull(orig.ident) %>%
  unique()

tab <- read.csv("/projects/b1169/boles/als_multitissue_scfrp/tab_data/metadata.csv")
tab$ID

tab %>%
  filter(ID %in% c("GWF19-47", "GWF20-54", "GWF21-56"))

obj@meta.data <- obj@meta.data %>%
  mutate(Batch = if_else(orig.ident %in% c("GWF-19-47_b", "GWF-20-54_b", "GWF-21-56_b"),
                         1, Batch),
         Group = if_else(orig.ident %in% c("GWF-19-47_b", "GWF-20-54_b", "GWF-21-56_b"),
                         "sALS", Group))

saveRDS(obj,
        file = paste0(data_out_dir, "obj.rds"))


