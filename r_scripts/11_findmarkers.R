library(Seurat)
library(scCustomize)
library(tidyverse)

proj_dir <- "/projects/b1169/boles/als_multitissue_scfrp/"

data_in_dir <- paste0(proj_dir, "data/11_clustering")
files_in <- list.files(data_in_dir,
                       full.names = T)

csv_dir <- paste0(proj_dir, "tab_data/12_findmarkers/")
dir.create(csv_dir,
           recursive = T,
           showWarnings = F)

obj_list <- list()
length(obj_list) <- 3

obj_list <- map(files_in, 
                readRDS)

# Already joined layers in script 11
# obj_list <- map(obj_list,
#                 JoinLayers)

res_choices <- c(1, 1.2, 1.4)

names(obj_list) <- c("brain", "muscle", "sc")

for (i in seq_along(obj_list)){
  Idents(obj_list[[i]]) <- paste0("res", res_choices[i], "_clusters")
  
  markers <- FindAllMarkers(obj_list[[i]])

  write.csv(markers,
            file = paste0(csv_dir, names(obj_list)[i],
                          "_res", res_choices[i], "_markers.csv"))
  
  Create_Cluster_Annotation_File(file_path = csv_dir,
                                 file_name = paste0(names(obj_list)[i], "_res", res_choices[i], "_annotations"))
                            
}
