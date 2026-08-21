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

# --mem/--time are carried over unchanged from the Harmony version of this
# job, itself unmeasured -- back on CCA now, which is the heaviest
# combination tried yet: full per-tissue data (not the 2000-cells/sample
# sketch this was last measured on) plus CCA's dense expression access
# (not Harmony's PCA-embedding-only approach). Watch this closely and check
# `seff <jobid>` after it runs -- if it's too slow, k.anchor (currently 20,
# vs. Seurat's default of 5) is the next lever to try lowering.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/09_muscle_integration.R
