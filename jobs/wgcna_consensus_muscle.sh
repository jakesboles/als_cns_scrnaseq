#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name wgcna_consensus_muscle
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%j.log
#SBATCH --verbose

# Single task, no --array -- there's exactly one consensus comparison
# here (denervated vs. Type I vs. Type II fibers), matching
# 06_obj_reassembly.sh's precedent for genuinely single-task work.
#
# --mem/--time carried over from wgcna_consensus.sh's sizing -- same kind
# of work (consensus metacell/network construction per group), just 3
# fiber-type groups instead of 2 tissue groups. Check `seff <jobid>` once
# this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/wgcna_consensus_muscle.R
