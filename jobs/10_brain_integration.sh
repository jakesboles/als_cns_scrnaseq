#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 10_brain_integration
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 64G
#SBATCH --time 6:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time are a starting-point estimate, not measured -- the previous
# 150G/48:00:00 here was sized for CCA on the full (non-sketch) per-sample
# data. This now runs on the sketch subset from 08_sketch_pca.R (2000
# cells/sample) with variable features actually set (see that script's
# comments -- likely a major reason the old runs were so slow), so both
# should be considerably lower, but by how much isn't measured yet. Check
# `seff <jobid>` after this runs and adjust -- if it's still slow, k.anchor
# (currently 20, vs. Seurat's default of 5) is the next lever to try
# lowering.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/10_brain_integration.R
