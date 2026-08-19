#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 08_sketch_pca
#SBATCH --array 1-3
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 10:00:00
#SBATCH --output /projects/b1169/boles/als_cns_scrnaseq/logs/%x_%A_%a.log
#SBATCH --verbose

# --array must match nrow(tissues) in 08_sketch_pca.R -- currently 3 (brain,
# spinal cord, muscle). 08_sketch_pca.R will fail fast with a clear error if
# these get out of sync.
#
# --mem/--time are a starting-point estimate, not measured -- the previous
# single-job version (all 3 tissues sequentially) timed out at 6:00:00
# partway through JackStraw() on the second tissue (spinal cord), so each
# tissue alone is somewhat under that, but not precisely known -- 10:00:00
# leaves margin above the observed partial run. --mem is unchanged from the
# sequential version since only one tissue was ever in memory at a time
# there too. Check `seff <jobid>_<taskid>` after the first tasks finish and
# adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_scrnaseq/r_scripts/08_sketch_pca.R
