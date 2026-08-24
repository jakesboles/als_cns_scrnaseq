#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name 14_findmarkers2
#SBATCH --array 1-20
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 40G
#SBATCH --time 12:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match the number of rows in jobs/13_params.txt -- currently
# 20 (one per tissue/cell-type combination), same list 13_subclustering1.R
# uses. 14_findmarkers2.R will fail fast with a clear error if these get
# out of sync.
#
# --mem/--time are an unmeasured estimate -- each task runs FindAllMarkers()
# plus a handful of plots on one cell type's cells (a subset of a subset,
# smaller than anything 11_findmarkers.R processed), so sized well below
# that job's 100G/24:00:00 as a starting point, not a measurement. Check
# `seff <jobid>_<taskid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/14_findmarkers2.R
