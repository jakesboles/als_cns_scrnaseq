#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 12_annotation1
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 60G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 12_annotation1.R -- currently 3
# (brain, spinal cord, muscle). 12_annotation1.R will fail fast with a
# clear error if these get out of sync.
#
# --mem/--time are an unmeasured estimate -- this is mostly rendering a
# UMAP, one dot plot, and ~20 rasterized feature/violin plot pairs, which
# should be lighter than 10/11's compute (clustering, marker testing), so
# sized down from those as a starting point, not a measurement. Check
# `seff <jobid>_<taskid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/12_annotation1.R
