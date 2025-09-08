#!/bin/sh

########## Configuration Slurm #######
## Define a job name
#SBATCH --job-name=prep-array-tt
## Define the partition
#SBATCH -p normal
## Define the number of cores to use
#SBATCH --cpus-per-task=16
## Define our email address
#SBATCH --mail-user=mathis.arnal@ird.fr
## Define when to receive the email (BEGIN, ABORT, END or ALL)
#SBATCH --mail-type=ALL
#######################################

# Create a personal directory in the scratch
cd /scratch
mkdir arnal-$SLURM_JOB_ID

# Retrieve data from the san
scp -r san:/users/arnal/VERDI_RECOVER/bootstrap_iqtree /scratch/arnal-$SLURM_JOB_ID


ALIGNMENT="no_outliers_aligned_verdiGISAIDref_seq.fasta"
DATES="no_outliers_treetime_verdiGISAIDref_df.csv"
UFBOOT="iqtree.ufboot"
OUTDIR="treetime_array"

cd /scratch/arnal-$SLURM_JOB_ID/bootstrap_iqtree
mkdir -p "$OUTDIR"

echo "Splitting bootstrap trees..."
# Split .ufboot into individual tree files
split -l 1 --numeric-suffixes=1 --suffix-length=4 "$UFBOOT" "$OUTDIR/tree_"

cd "$OUTDIR"
mkdir -p res

# Loop through each tree file
for TREEFILE in tree_*; do
    BOOTNUM=$(basename "$TREEFILE" | sed 's/tree_//')
    RUNDIR="run_${BOOTNUM}"
    mkdir -p "$RUNDIR"
    mv "$TREEFILE" "$RUNDIR/tree.nwk"
done

# Copy the results back to the san
scp -r /scratch/arnal-$SLURM_JOB_ID/bootstrap_iqtree/treetime_array san:/users/arnal/VERDI_RECOVER/treetime_array

