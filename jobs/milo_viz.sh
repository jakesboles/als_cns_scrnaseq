#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name milo_viz
#SBATCH --array 1-6
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 32G
#SBATCH --time 2:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match the number of subdirectories in data/milo/ --
# currently 6, matching milo.R's own targets. milo_viz.R discovers
# targets from disk at runtime (see its header) rather than a hardcoded
# list, so this only needs updating if milo.R's own targets change --
# milo_viz.R will fail fast with a clear error if these fall out of sync.
#
# --mem/--time are an unmeasured estimate -- this only loads an already-
# built Milo object and a saved UMAP embedding and draws plots, no
# FindNeighbors()/graph-building/expression data involved, so sized well
# below milo.R's own job as a starting point, not a measurement. Check
# `seff <jobid>_<taskid>` once these run and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/milo_viz.R
