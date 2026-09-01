# One-time helper that generates jobs/wgcna_single_params.txt: every
# (cell_type, tissue) combination that exists in 17_obj_reassembly.R's
# output. wgcna_single.R runs hdWGCNA on literally every cell type in
# every tissue -- not a deliberately-curated subset like 13/15/19's
# targets -- so the params file is derived straight from the data rather
# than hand-written.
#
# Run this once as a plain Rscript call (not an array job) whenever
# 17_obj_reassembly.R's output or cell_type3 annotations change, then set
# jobs/wgcna_single.sh's --array range to match the row count printed at
# the end.

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
