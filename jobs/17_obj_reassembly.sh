#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 17_obj_reassembly
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 17_obj_reassembly.R -- currently 3
# (brain, spinal cord, muscle). 17_obj_reassembly.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time carried over from 07_norm_pca.R/09_*_integration.R's
# per-tissue sizing (this step does the same kind of work -- Normalize/
# FindVariableFeatures/ScaleData/RunPCA(npcs=100)/Harmony -- on the full
# tissue, minus JackStraw), scaled down slightly since the cleaned object
# has fewer cells than the original (cells removed via cell_type3 ==
# "Remove"). Unmeasured starting point -- check `seff <jobid>_<taskid>`
# after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/17_obj_reassembly.R
