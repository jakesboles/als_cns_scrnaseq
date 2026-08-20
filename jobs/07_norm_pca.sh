#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_norm_pca
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 07_norm_pca.R -- currently 3 (brain,
# spinal cord, muscle). 07_norm_pca.R will fail fast with a clear error if
# these get out of sync.
#
# --mem/--time are unchanged per-task from the single-job version (still a
# rough, unmeasured estimate) -- each task now handles one tissue's full
# cell count instead of the job looping over all 3 sequentially, so this
# should only reduce total wall-clock time (3 tissues in parallel instead
# of one after another), not change what any single tissue needs. JackStraw
# (100 replicates x 100 dims) is still the biggest unknown -- watch the
# first tasks closely and check `seff <jobid>_<taskid>` once they finish.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/07_norm_pca.R
