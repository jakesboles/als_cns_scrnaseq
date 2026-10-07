#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 13_subclustering1
#SBATCH --array 1-20
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 80G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/13_subclustering1.R
