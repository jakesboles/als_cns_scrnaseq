#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 18_full_integration
#SBATCH --array 1-2
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 250G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match length(integration_targets) in 18_full_integration.R
# -- currently 2 (all 3 tissues, and brain+sc only). 18_full_integration.R
# will fail fast with a clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- this does the same kind of
# work as 07_norm_pca.R/17_obj_reassembly.R (Normalize/
# FindVariableFeatures/ScaleData/RunPCA(npcs=100)/Harmony/neighbors/UMAP)
# but on a merged multi-tissue population, so sized up from those as a
# starting point, not a measurement. Task 1 (all_tissues) has meaningfully
# more cells than task 2 (brain_sc) -- watch it closely and check
# `seff <jobid>_<taskid>` once both finish.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/18_full_integration.R
