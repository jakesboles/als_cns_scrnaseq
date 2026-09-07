#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name deseq2_gsea
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/deseq2_gsea.R
