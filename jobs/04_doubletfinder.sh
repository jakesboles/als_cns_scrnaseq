#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 04_doubletfinder
#SBATCH --array 1-90
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 32G
#SBATCH --time 3:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match the number of samples in
# sort(unique(readRDS("data/03_qc2/metadata.rds")$orig.ident)) -- currently
# 90. 04_doubletfinder.R will fail fast with a clear error if these get out
# of sync.
#
# --mem/--ntasks-per-node/--time are a starting-point guess for one sample's
# worth of clustering + DoubletFinder (much lighter than the whole-cohort
# jobs in 01/02) -- there's no profiling behind these numbers yet, so check
# `seff <jobid>_<taskid>` after the first handful of tasks finish and adjust.
#
# To cap how many tasks run at once (considerate on a shared cluster), use
# e.g. `#SBATCH --array 1-90%10` for at most 10 concurrent tasks.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/04_doubletfinder.R
