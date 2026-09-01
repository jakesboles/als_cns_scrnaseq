#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 19_subclustering3
#SBATCH --array 1-6
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match length(subclustering_targets) in 19_subclustering3.R
# -- currently 6 (microglia, astrocyte, oligodendrocyte, neurons,
# muscle_fiber, myeloid). 19_subclustering3.R will fail fast with a clear
# error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- same kind of work as
# 13_subclustering1.R/15_subclustering2.R (Normalize/FindVariableFeatures/
# ScaleData/RunPCA(npcs=50)/Harmony/UMAP) but some targets here span 2-3
# tissues instead of one, so sized up from those as a starting point, not
# a measurement. Targets vary a lot in cell count (e.g. "neurons" is
# likely much larger than "oligodendrocyte") -- check
# `seff <jobid>_<taskid>` per task once these finish and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/19_subclustering3.R
