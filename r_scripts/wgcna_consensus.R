suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(UCell)
  library(cowplot)
})

options(future.globals.maxSize=1048576000000)
load("/projects/b1169/projects/sea_ad_hypothalamus/results/preprocessing/qc/out_TW_05-04-2023/helperfunctions.RData")

dir <- "/gpfs/projects/b1169/thomas/als_multitissue/hdWGCNA/Consensus_jsb/"
dir.create(dir,
           showWarnings = F,
           recursive = T)
setwd(dir)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(256)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 16)

# Load Seurat object
s_all <- readRDS("/gpfs/projects/b1169/boles/als_multitissue_scfrp/data/23_full_projection/full_integrated.rds")

meta <- read.csv("/gpfs/projects/b1169/thomas/als_multitissue/hdWGCNA/AllMeta.csv")

# colnames(s_all@meta.data)

# colnames(meta)

changed <- rownames(s_all@meta.data[rownames(s_all@meta.data) %!in% meta$X,]) # sketch assay tags barcodes with an underscore

changed2 <- sub("_", "", changed)

sum(changed2 %in% meta$X)

meta$X[meta$X %in% changed2] <- paste0("_", meta$X[meta$X %in% changed2])

sum(meta$X %in% rownames(s_all@meta.data)) == nrow(meta)

rownames(meta) <- meta$X

meta <- meta[rownames(s_all@meta.data),]

all.equal(meta$X, rownames(s_all@meta.data))

s_all@meta.data <- s_all@meta.data[,c(1:8)]

s_all@meta.data <- cbind(s_all@meta.data, meta[,10:24])

# Join layers
s <- JoinLayers(s_all)

# Scale data for eigengene harmonization
s <- ScaleData(s)

s <- subset(s, tissue %in% c("Cervical spinal cord", "Motor cortex"))

cell = "Microglia"

dir.create(cell)

dir.create(paste0(cell, "/plots"))

dir.create(paste0(cell, "/csvs"))

setwd(cell)

s2 <- subset(s, PredictedCellType == cell)

genes = rownames(s2)

# s2 <- JoinLayers(s2)

pe <- as.data.frame(rowMeans(GetAssayData(s2, layer = "data", assay = "RNA") > 0))

colnames(pe) <- "Percent"

pe$gene <- rownames(pe)

pe <- as.data.frame(pe[pe$Percent > 0.05,])

fk <- pe$gene

gc()
# Set up for WGCNA
seurat_obj <- SetupForWGCNA(s, 
                            gene_select = "custom",
                            features = fk,
                            # fraction = 0.05,  
                            wgcna_name = cell)

# Group Metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("PredictedCellType", "orig.ident", "tissue"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'full_integrated_pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'PredictedCellType' # set the Idents of the metacell seurat object
)

# Normalize and scale metacells - must scale for harmonized eigengenes
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- ScaleMetacells(seurat_obj, 
                             features=VariableFeatures(seurat_obj))
seurat_obj <- RunPCAMetacells(seurat_obj, 
                              features=VariableFeatures(seurat_obj))
seurat_obj <- RunHarmonyMetacells(seurat_obj, 
                                  group.by.vars='orig.ident')

# Set up expression matrix
seurat_obj <- SetMultiExpr(
  seurat_obj,
  group_name = cell, # the name of the group of interest in the group.by column
  group.by='PredictedCellType',
  multi.group.by = "tissue",
  multi_groups = c("Cervical spinal cord", "Motor cortex"),
  assay = "RNA",
  slot = "data", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  use_metacells=T
)

# run soft power test
seurat_obj <- TestSoftPowersConsensus(seurat_obj)

# generate plots
plot_list <-  PlotSoftPowers(seurat_obj)

# saveRDS(plot_list, "plotlist.rds")

# saveRDS(seurat_obj, "seurat_obj.rds")

# seurat_obj <- readRDS("plotlist.rds")

# get just the scale-free topology fit plot for each group
consensus_groups <- unique(seurat_obj$tissue)
p_list <- lapply(1:length(consensus_groups), function(i){
  cur_group <- consensus_groups[[i]]
  plot_list[[i]][[1]] + ggtitle(paste0('Tissue: ', cur_group)) + theme(plot.title=element_text(hjust=0.5))
})

pdf(paste0("plots/SoftPower_", cell, ".pdf"), height = 8, width = 8)
wrap_plots(p_list, ncol=2)
dev.off()

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table,
          file = paste0("csvs/soft_powers.csv"))

# soft_power <- min(plot_list[[1]]$data$Power[plot_list[[1]]$data$SFT.R.sq > 0.8])

# build consensus network
seurat_obj <- ConstructNetwork(seurat_obj, 
                               # soft_power = c(8, 10), 
                               # tom_outdir = paste0(cell, "/"), 
                               tom_name = cell,
                               consensus = T,
                               overwrite_tom = TRUE)

# plot dendogram
pdf(paste0("plots/Dendrogram", "_", cell, ".pdf"), height = 8, width = 8)
PlotDendrogram(seurat_obj, main= paste0(cell, " Dendrogram"))
dev.off()

# change active hdWGCNA experiment to consensus
seurat_obj <- SetActiveWGCNA(seurat_obj, cell)

# compute eigengenes and connectivity
seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj, group_name = cell, group.by = "PredictedCellType") # Test to see what changing this from group to cell type does

# visualize network with semi-supervised UMAP
# seurat_obj <- RunModuleUMAP(
#   seurat_obj,
#   n_hubs = 5,
#   n_neighbors=5,
#   min_dist=0.1,
#   spread=2,
#   wgcna_name = cell,
#   target_weight=0.05,
#   supervised=TRUE
# )
# 
# # the ooooo aaahhhh plot
# pdf(paste0(cell, "/plots/UselessModuleUMAP", "_", cell, ".pdf"), height = 8, width = 8)
# print(ModuleUMAPPlot(
#   seurat_obj,
#   edge.alpha=0.5,
#   sample_edges=TRUE,
#   keep_grey_edges=FALSE,
#   edge_prop=0.075,
#   label_hubs=0
# ))
# dev.off()

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=4, text_size = 4)
pdf(paste0("plots/ModuleConnectivity", "_", cell, ".pdf"), height = 12, width = 12)
print(p)
dev.off()

# Compute module expression scores
library(UCell)
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# colnames(seurat_obj)

# Get harmonized module eigengenes
hMEs <- GetMEs(seurat_obj, harmonized = T)

# write module eigengenes CSV
write.csv(hMEs, paste0("csvs/ModuleEigengenes.csv"))

# UMAP of module eigengenes
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='MEs', # plot the hMEs
  reduction = "full_integrated_pca",
  order = TRUE,
  ucell = T
  # order so the points with highest hMEs are on top
)

# stitch together with patchwork
pdf(paste0("plots/ModuleEigengenesUMAP", "_", cell, ".pdf"), height = 8, width = 8)
print(wrap_plots(plot_list, ncol=4))
dev.off()

# UMAP of module expression scores
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='scores', # plot the hMEs
  reduction = "full_integrated_pca",
  order = TRUE,
  ucell = T
  # order so the points with highest scores are on top
)

# stitch together with patchwork
pdf(paste0("plots/ModuleScoresUMAP", "_", cell, ".pdf"), height = 8, width = 8)
print(wrap_plots(plot_list, ncol=4))
dev.off()

# Module eigengene correlogram
pdf(paste0("plots/ModuleEigengeneCorr", "_", cell, ".pdf"), height = 8, width = 8)
ModuleCorrelogram(seurat_obj, features = "MEs")
dev.off()

# Module score correlogram
pdf(paste0("plots/ModuleScoresCorr", "_", cell, ".pdf"), height = 8, width = 8)
ModuleCorrelogram(seurat_obj, features = "scores")
dev.off()

# Get data frame with gene-module memberships
mods <- seurat_obj@misc[[cell]][["wgcna_modules"]]

# Save for later
write.csv(mods, paste0("csvs/Modules.csv"))

# Save seurat object
saveRDS(seurat_obj, paste0(cell, ".rds"))

# # remove seurat_obj and s2
# rm(seurat_obj, s2)
# 
# # Be nice to the cluster
# gc()