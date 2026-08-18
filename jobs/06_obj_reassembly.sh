#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 06_obj_reassembly
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 200G
#SBATCH --time 4:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time bumped up from the previous 64G/2:00:00 -- this job does the
# same kind of merge as 01_obj_creation.sh (which needed 240G for the raw
# ~1.3M-cell cohort), just on the somewhat smaller post-QC/post-doublet-
# filtering cell count. 200G/4:00:00 is a starting-point estimate, not
# measured -- check `seff <jobid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/06_obj_reassembly.R
