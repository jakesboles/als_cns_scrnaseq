#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 15_annotation2
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 40G
#SBATCH --time 6:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 15_annotation2.R -- currently 3
# (brain, spinal cord, muscle). 15_annotation2.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- this only loads one tissue's
# expression data/UMAP plus a handful of lightweight metadata/CSV files
# per cell type and makes one plot, no PCA/clustering/marker-finding, so
# sized well below the other per-tissue jobs as a starting point, not a
# measurement. Check `seff <jobid>_<taskid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/15_annotation2.R
