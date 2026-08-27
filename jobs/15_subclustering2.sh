#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 15_subclustering2
#SBATCH --array 1-13
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match length(subclustering_targets) in 15_subclustering2.R
# -- currently 13 (5 brain + 4 muscle + 4 sc targets, some of which are
# multi-cell-type groups). 15_subclustering2.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate, carried over from
# 13_subclustering1.R's per-cell-type sizing (80G/24:00:00) -- this task
# does the same kind of work (PCA/Harmony/neighbors/UMAP/13-resolution
# clustering/FindAllMarkers) on one cell_type2 group's cells, but some of
# these groups combine multiple cell types and so may be larger than any
# single 13_subclustering1.R subset. Watch the first few tasks closely
# (especially the multi-cell-type groups) and check `seff <jobid>_<taskid>`
# once they finish.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/15_subclustering2.R
