#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name wgcna_consensus
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 24:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match length(subclustering_targets) in wgcna_consensus.R
# -- currently 3 (Microglia, Oligodendrocyte, Astrocyte).
# wgcna_consensus.R will fail fast with a clear error if these get out of
# sync.
#
# --mem/--time carried over from wgcna_single.sh's sizing (this project's
# own measured starting point from that script's completed run) --
# consensus WGCNA does the same kind of work per tissue group rather than
# once, so this could need more, not less. Check `seff <jobid>_<taskid>`
# per task once these run and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/wgcna_consensus.R
