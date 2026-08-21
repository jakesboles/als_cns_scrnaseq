#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 10_clustering
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 10_clustering.R -- currently 3
# (brain, spinal cord, muscle). 10_clustering.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate, carried over from
# 07_norm_pca.R's per-tissue sizing -- this script adds a neighbor
# graph/UMAP plus Leiden clustering at 19 resolutions, which 07 didn't
# need to do. The motor cortex task previously OOM-killed at 200G during
# silhouette scoring, because that step built an O(cells^2) distance
# matrix -- cluster scoring now uses centroid/graph-based approximations
# (see the top of 10_clustering.R) that scale with cells x clusters
# instead, so this should no longer be the memory bottleneck, but the
# actual ceiling is still unmeasured. Watch the first tasks closely
# (especially the largest tissue) and check `seff <jobid>_<taskid>` once
# they finish.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/10_clustering.R
