#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 13_subclustering1
#SBATCH --array 1-20
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match the number of rows in jobs/13_params.txt -- currently
# 20 (one per tissue/cell-type combination). 13_subclustering1.R will fail
# fast with a clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- each task re-runs
# PCA/Harmony/neighbors/UMAP/13-resolution Leiden clustering (the same
# kind of work 10_clustering.R does per tissue), but on one cell type's
# cells instead of a whole tissue, so this should need meaningfully less
# than 10_clustering.R's 200G/48:00:00. Sized down as a starting point,
# not a measurement -- watch the first few tasks (especially any large or
# very rare cell types) and check `seff <jobid>_<taskid>` once they
# finish.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/13_subclustering1.R
