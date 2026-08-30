#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name deseq2
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in deseq2.R -- currently 3 (brain,
# spinal cord, muscle). deseq2.R will fail fast with a clear error if
# these get out of sync.
#
# --mem/--time are an unmeasured estimate -- each task loads one tissue's
# full raw count matrix once, then loops over cell_type3 pseudobulking
# and running DESeq2 per cell type (much cheaper per-iteration than a
# PCA/Harmony/UMAP step, but the full-tissue BPCells load dominates
# memory). Sized well below 17_obj_reassembly.R's 200G/48:00:00 as a
# starting point, not a measurement. Check `seff <jobid>_<taskid>` after
# this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/deseq2.R
