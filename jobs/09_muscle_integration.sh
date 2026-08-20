#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 09_muscle_integration
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 150G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time are carried over unchanged from the previous (CCA, sketch
# data) version of this job -- genuinely uncertain now given two things
# changed at once: this runs on the full per-tissue data from 07_norm_pca.R
# instead of the 2000-cells/sample sketch, but Harmony integrates in the
# already-small PCA embedding space rather than needing dense access to
# expression data the way CCA did, so it may not need as much. Check
# `seff <jobid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/09_muscle_integration.R
