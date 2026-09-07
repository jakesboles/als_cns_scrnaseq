#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics
#SBATCH --job-name milo
#SBATCH --array 2,3,6
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 8
#SBATCH --mem 64G
#SBATCH --time 8:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match the number of subdirectories in data/19_subclustering3/
# -- currently 6 (microglia, astrocyte, oligodendrocyte, neurons,
# muscle_fiber, myeloid). milo.R discovers targets from disk at runtime
# (see its header) rather than a hardcoded list, so this only needs
# updating if a target is added to or removed from 19_subclustering3.R --
# milo.R will fail fast with a clear error if these fall out of sync.
#
# --mem/--time are an unmeasured estimate -- MiloR's neighborhood testing
# doesn't touch raw expression data, so sized well below the heavier
# integration/WGCNA jobs as a starting point, not a measurement. The
# neighbor graph and neighborhoods are now built once per task on the
# full cross-tissue object (not per-tissue subsets), so "myeloid" (3
# tissues combined into one object) is likely the slowest/highest-memory
# target here. Check `seff <jobid>_<taskid>` once these run and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/milo.R
