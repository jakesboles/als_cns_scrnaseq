library(tidyverse)
library(DESeq2)
library(apeglm)
library(IHW)
library(Seurat)

setwd("/projects/b1169/boles/als_cns_scrnaseq")

results_dir <- "results/deseq2/"
dir.create(results_dir)

# Load BPCells object from 17 obj reassembly

# Factorize grouping variable

obj@meta.data$group <- factor(obj@meta.data$group,
                              levels = c("Control", "sALS", "C9orf72"))


# Extract cell type labels to establish for loop

cell_types <- unique(obj$cell_type3)

for (i in seq_along(cell_types)){
  
  message(paste0(celltypes[i], " in ", tissue))
  
  sub <- subset(obj,
                cell_type3 == cell_types[i])
  
  file <- str_replace_all(celltypes[i], " ", "_")
  
  bulk <- AggregateExpression(sub,
                              assays = "RNA",
                              return.seurat = F,
                              slot = "counts",
                              group.by = c("orig.ident"))
  
  exp <- bulk$RNA
  
  meta <- sub@meta.data %>%
    dplyr::select(c(orig.ident, group, sex, age)) %>% # make sure this is age at death
    distinct()
  
  idx <- match(colnames(exp), rownames(meta))
  meta <- meta[idx, ]
  
  meta$Genotype <- factor(meta$Genotype,
                          levels = c("WT", "HET", "HMM"))
  
  
  dds <- DESeqDataSetFromMatrix(countData = exp,
                                colData = meta,
                                design = ~ group) # change this as needed
  
  keep <- rowSums(counts(dds) >= 5) >= 4 # change these cutoffs as needed
  
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  saveRDS(dds,
          file = paste0(resuls_dir, file, "_dds.rds"))
  
  # resultsNames(dds)
  
  dir.create(paste0(results_dir, file))
  
  res <- results(dds,
                 contrast = c("group", "sALS", "Control"),
                 filterFun = ihw,
                 independentFilter = T)
  
  res <- as.data.frame(res)
  
  write.csv(res,
            file = paste0(results_dir, file, "/sALS_vs_Control.csv"))
  
  res_shrunk <- lfcShrink(dds,
                          coef = "group_sALS_vs_Control",
                          type = "apeglm")
  
  write.csv(res_shrunk,
            file = paste0(results_dir, file, "/sALS_vs_Control_lfc_shrunk.csv"))
  
  res <- results(dds,
                 contrast = c("group", "C9orf72", "ALS"),
                 filterFun = ihw,
                 independentFilter = T)
  
  res <- as.data.frame(res)
  
  write.csv(res,
            file = paste0(results, "/C9orf72_vs_Control.csv"))
  
  res_shrunk <- lfcShrink(dds,
                          coef = "group_C9orf72_vs_Control",
                          type = "apeglm")
  
  write.csv(res_shrunk,
            file = paste0(results, "/C9orf72_vs_Control_lfc_shrunk.csv"))
  
}
