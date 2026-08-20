#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_norm_pca
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time are a rough starting-point estimate, not measured, and
# genuinely uncertain -- this now runs NormalizeData/ScaleData/RunPCA/
# JackStraw on each tissue's FULL cell count (tens to hundreds of thousands
# of cells per tissue), not the 2000-cells/sample sketch the old 07/08 used.
# JackStraw (100 replicates x 100 dims) is the biggest unknown: its cost
# scales with cell count, and it was only ever run at sketch scale before,
# so there's no prior data point for how long it takes at full scale. Watch
# the first run closely -- if JackStraw is the bottleneck, num.replicate is
# the parameter to reconsider (a scientific tradeoff, not changed here).
# --mem sized similarly to 06_obj_reassembly.sh's 200G (same order of
# magnitude of cells, just one tissue instead of the whole cohort at once).
# Check `seff <jobid>` after this runs and adjust both values.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/07_norm_pca.R
