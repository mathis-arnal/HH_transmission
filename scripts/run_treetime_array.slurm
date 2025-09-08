#!/bin/bash
#SBATCH --job-name=tt_array
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --nodelist= =nodes[20,24]
#SBATCH --time=1:00:00
#SBATCH --output=slurm_logs/tt_%A_%a.out
#SBATCH --error=slurm_logs/tt_%A_%a.err

module load singularity

# Parameters
WORKDIR="san:/users/arnal/VERDI_RECOVER/treetime_array"
BOOTSTRAP_DIR="san:/users/arnal/VERDI_RECOVER/bootstrap_iqtree"
ALIGNMENT="no_outliers_aligned_verdiGISAIDref_seq.fasta"
DATES="no_outliers_treetime_verdiGISAIDref_df.csv"
ID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)
TREE="$WORKDIR/res_2/trees/run_${ID}.nwk"


# Create task-specific working directory
SCRATCHDIR="/scratch/arnal-$SLURM_JOB_ID-run_${ID}"
mkdir -p "$SCRATCHDIR"

# Copy tree.nwk for this replicate
scp "$BOOTSTRAP_DIR/$DATES" "$SCRATCHDIR/"
scp "$BOOTSTRAP_DIR/$ALIGNMENT" "$SCRATCHDIR/"
scp $TREE "$SCRATCHDIR/"

cd "$SCRATCHDIR"

# Run TreeTime
singularity run https://depot.galaxyproject.org/singularity/treetime:0.11.4--pyhdfd78af_0 \
    treetime --tree tree.nwk \
    --aln "$SHAREDIR/$ALIGNMENT" \
    --dates "$SHAREDIR/$DATES" \
    --reroot NC_045512.2 \
    --stochastic-resolve \
    --clock-rate 0.001 \
    --clock-std-dev 0.0005
# Save output
scp -r "$SCRATCHDIR/treetime/timetree.nexus" "san:/users/arnal/VERDI_RECOVER/treetime_array/res_2/treetime_run_${ID}.nexus"


# Cleanup just task-specific files
rm -rf "$SCRATCHDIR"
