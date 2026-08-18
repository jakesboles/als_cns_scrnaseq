#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07_norm
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time are a starting-point estimate, not measured -- this reads
# through the full reassembled cohort (data/06_obj_reassembly/bpcells) for
# NormalizeData()/FindVariableFeatures() on each tissue subset in turn, so
# sized similarly to 06_obj_reassembly.sh. Check `seff <jobid>` after this
# runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/07_norm.R
