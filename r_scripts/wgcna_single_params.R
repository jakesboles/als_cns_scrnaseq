# One-time helper that generates jobs/wgcna_single_params.txt: every
# (cell type, tissue) combination in 17_obj_reassembly.R's output, for
# wgcna_single.R's job array. Run by hand (not an array job) whenever
# cell_type3 annotations change.

setwd("/projects/b1169/boles/als_cns_scrnaseq")

tissues <- c("brain", "sc", "muscle")

params <- do.call(rbind, lapply(tissues, function(t){
  meta <- readRDS(paste0("data/17_obj_reassembly/", t, "/metadata.rds"))
  data.frame(cell_type = unique(meta$cell_type3), tissue_file = t)
}))

write.table(params,
            file = "jobs/wgcna_single_params.txt",
            sep = ",", row.names = F, col.names = F, quote = F)

message(paste0(nrow(params), " (cell_type, tissue) combinations written to ",
               "jobs/wgcna_single_params.txt -- set jobs/wgcna_single.sh's ",
               "--array range to 1-", nrow(params), "."))
