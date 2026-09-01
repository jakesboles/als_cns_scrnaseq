suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(UCell)
  library(cowplot)
})

# Make directories, define variables --------------------------------------

cell <- commandArgs(trailingOnly = TRUE)[1]
tissue <- commandArgs(trailingOnly = T)[2]
file <- commandArgs(trailingOnly = T)[3]

cell_file <- str_replace_all(cell, "/", "-")

options(future.globals.maxSize=1048576000000)
load("/projects/b1169/projects/sea_ad_hypothalamus/results/preprocessing/qc/out_TW_05-04-2023/helperfunctions.RData")
setwd("/gpfs/projects/b1169/thomas/als_multitissue/hdWGCNA/")

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

dir.create(paste0("25_", file, "/", cell_file),
           recursive = T,
           showWarnings = F)

setwd(paste0("/gpfs/projects/b1169/thomas/als_multitissue/hdWGCNA/25_", file, "/", cell_file))

dir.create("plots",
           showWarnings = F)

dir.create("csvs",
           showWarnings = F)

# set random seed for reproducibility
set.seed(256)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 16)

# Load and prep full tissue object ----------------------------------------

s <- readRDS(paste0("/gpfs/projects/b1169/boles/als_multitissue_scfrp/data/24_reintegration/harmony/", file, "_reintegrated.rds"))
DefaultAssay(s) <- "RNA"
s[["sketch"]]<- NULL
s@graphs[c("sketch_nn", "sketch_snn")] <- NULL

s@meta.data$final_label2 <- str_replace_all(s@meta.data$final_label2, "/", "-")
# Using PredictedCellType throughout the script so it's easier to load final labels into it
s@meta.data$PredictedCellType <- s@meta.data$final_label2

# same for sample_ID
s@meta.data$sample_ID <- s@meta.data$id

# Make sure Batch is a factor
s@meta.data$Batch <- as.factor(s@meta.data$Batch)

# Create new group variable for when comparison is ALS vs Control
s@meta.data$Group2 <- s@meta.data$Group

# Combine ALS Groups
s@meta.data$Group2[s@meta.data$Group2 %in% c("C9orf72", "sALS")] <- "ALS"
# 
# # # Add new metadata
# newmeta <- read.csv("/gpfs/projects/b1169/thomas/als_multitissue/Metadata/SpinalCord_NewMeta.csv", row.names = 1)
# 
# newmeta <- newmeta[rownames(newmeta) %in% rownames(s_all@meta.data),]
# 
# # colnames(newmeta)[16] <- "sample_ID2"
# 
# # Make sure barcodes match
# all.equal(rownames(s_all@meta.data), rownames(newmeta))
# 
# # Get rid of extra columns
# s_all@meta.data <- s_all@meta.data[,c(1:8, 80:82)]
# 
# # Add new meta to meta
# s_all@meta.data <- cbind(s_all@meta.data, newmeta[,c(5:15)])
# 
# # Using PredictedCellType throughout the script so it's easier to load final labels into it
# s_all@meta.data$PredictedCellType <- s_all@meta.data$final_label2

# Join layers
s <- JoinLayers(s)

# Scale data for eigengene harmonization
s <- ScaleData(s)

# Identify genes expressed in at least 5% of the cells of interest --------

s2 <- subset(s, PredictedCellType == cell_file)

genes = rownames(s2)

# s2 <- JoinLayers(s2)

pe <- as.data.frame(rowMeans(GetAssayData(s2, layer = "data", assay = "RNA") > 0))

colnames(pe) <- "Percent"

pe$gene <- rownames(pe)

pe <- as.data.frame(pe[pe$Percent > 0.05,])

fk <- pe$gene

# Set up for WGCNA
seurat_obj <- SetupForWGCNA(s, 
                            gene_select = "custom",
                            features = fk,
                            # fraction = 0.05,  
                            wgcna_name = cell_file)

# construct metacells in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("PredictedCellType", "sample_ID"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'integrated_pca', # select the dimensionality reduction to perform KNN on
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
                                  group.by.vars='sample_ID')

# Set up expression matrix
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cell_file, # the name of the group of interest in the group.by column
  group.by='PredictedCellType',
  assay = "RNA",
  slot = "data", # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  use_metacells=T
)

# Find soft power ---------------------------------------------------------

# run soft power test
seurat_obj <- TestSoftPowers(seurat_obj)

# generate soft power plots
plot_list <-  PlotSoftPowers(seurat_obj)

png(paste0("plots/SoftPower_", cell_file, ".png"), 
    height = 8, width = 8,
    units = "in",
    res = 600)
print(wrap_plots(plot_list, ncol=2))
dev.off()

power_table <- GetPowerTable(seurat_obj)
write.csv(power_table,
          file = paste0("csvs/soft_powers.csv"))

# soft_power <- min(plot_list[[1]]$data$Power[plot_list[[1]]$data$SFT.R.sq > 0.8])

# Build TOM & cluster genes -----------------------------------------------------------

seurat_obj <- ConstructNetwork(seurat_obj, 
                               # soft_power = soft_power, 
                               # tom_outdir = paste0(cell, "/"), 
                               tom_name = cell_file,
                               overwrite_tom = TRUE)
# letting this select the soft power automatically

# plot dendogram
png(paste0("plots/Dendrogram", "_", cell_file, ".png"), 
    height = 8, width = 8,
    units = "in",
    res = 600)
PlotDendrogram(seurat_obj, main= paste0(cell_file, " Dendrogram"))
dev.off()

# Module stats ------------------------------------------------------------

seurat_obj <- SetActiveWGCNA(seurat_obj, cell_file)

# compute eigengenes and connectivity
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = "sample_ID")
seurat_obj <- ModuleConnectivity(seurat_obj, group_name = cell_file, group.by = "PredictedCellType")

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=4, text_size = 4)
png(paste0("plots/ModuleConnectivity", "_", cell_file, ".png"), 
    height = 12, width = 12,
    units = "in",
    res = 600)
print(p)
dev.off()

# Compute module expression scores
# Get data frame with gene-module memberships
mods <- seurat_obj@misc[[cell_file]][["wgcna_modules"]]

# Save for later
write.csv(mods, paste0("csvs/Modules.csv"))

gene_sets <- list()

module_names <- unique(mods$module)
module_names <- module_names[!(module_names %in% "grey")]

for (i in module_names){
  gene_sets[[i]] <- mods %>%
    filter(module == i) %>%
    dplyr::pull(gene_name)
}

maxrank <- max(unlist(lapply(gene_sets, length)))

sub <- subset(seurat_obj,
              PredictedCellType == cell_file)

sub <- AddModuleScore_UCell(sub,
                            features = gene_sets,
                            maxRank = maxrank)

sub <- SmoothKNN(sub,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "integrated_pca")

scores <- sub@meta.data %>%
  dplyr::select(matches("UCell_kNN|Group|orig.ident"))

write.csv(scores,
          paste0("csvs/ModuleScoresUCell.csv"))

# Get harmonized module eigengenes
hMEs <- GetMEs(seurat_obj, harmonized = T)

# write module eigengenes CSV
write.csv(hMEs, 
          paste0("csvs/ModuleEigengenes.csv"))

# Module eigengene correlogram
png(paste0("plots/ModuleEigengeneCorr", "_", cell_file, ".png"), 
    height = 8, width = 8,
    res = 600,
    units = "in")
ModuleCorrelogram(seurat_obj, features = "MEs")
dev.off()

saveRDS(seurat_obj,
        file = "wgcna_obj.rds")

