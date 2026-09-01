#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name wgcna_single
#SBATCH --array 1-TODO
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 96G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array MUST be set before submitting -- run
# `Rscript r_scripts/wgcna_single_params.R` first (plain Rscript, not an
# array job) to (re)generate jobs/wgcna_single_params.txt from
# 17_obj_reassembly.R's current cell_type3 annotations, then replace
# 1-TODO above with 1-<row count printed by that script>.
# wgcna_single.R will fail fast with a clear error if this ever falls out
# of sync with jobs/wgcna_single_params.txt.
#
# --mem/--time are an unmeasured estimate -- MetacellsByGroups()/
# ConstructNetwork() run on the whole tissue's Harmony embedding (not a
# pre-subset object, see wgcna_single.R's header), so cost scales with
# tissue size more than with the specific cell type's own population.
# Sized above 13_subclustering1.R/15_subclustering2.R's per-cell-type
# scale as a starting point, not a measurement -- check
# `seff <jobid>_<taskid>` per task once these run and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/wgcna_single.R
