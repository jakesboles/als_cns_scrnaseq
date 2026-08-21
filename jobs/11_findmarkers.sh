#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 11_findmarkers
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 100G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 11_findmarkers.R -- currently 3
# (brain, spinal cord, muscle). 11_findmarkers.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- FindAllMarkers() (default
# Wilcoxon rank-sum test, one cluster vs. rest) is lighter than the
# neighbor graph/UMAP/clustering work in 10_clustering.R, so this is set
# lower than that job's 200G/48:00:00 as a starting point, not a
# measurement. Check `seff <jobid>_<taskid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/11_findmarkers.R
