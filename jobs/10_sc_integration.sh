#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 10_sc_integration
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 150G
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/10_sc_integration.R
