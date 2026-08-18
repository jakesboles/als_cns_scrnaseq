#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 06_obj_assembly
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 64G
#SBATCH --time 2:00:00
#SBATCH --output /projects/p31535/boles/als_multitissue_scfrp/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.4.0

Rscript /projects/p31535/boles/als_multitissue_scfrp/r_scripts/06_obj_assembly.R
