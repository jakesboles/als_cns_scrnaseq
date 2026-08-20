#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 10_brain_cca
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 150G
#SBATCH --time 48:00:00
#SBATCH --output /projects/p31535/boles/als_multitissue_scfrp/logs/%x.oe%j.log
#SBATCH --verbose

module load R/4.4.0

Rscript /projects/p31535/boles/als_multitissue_scfrp/r_scripts/10_brain_cca_integration.R
