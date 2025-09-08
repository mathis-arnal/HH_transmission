#!/bin/sh

########## configuration Slurm #######
## Définir un nom de job
#SBATCH --job-name=iq-tree
## Définir la partition
#SBATCH -p normal
## Definir le nombre de coeur a utiliser
#SBATCH --cpus-per-task=16
## définir notre adresse-mail
#SBATCH --mail-user=mathis.arnal@ird.fr
##  Définir qd recevoir le mail(BEGIN, ABORT, END ou ALL)
#SBATCH  --mail-type=ALL
#######################################

# Créer un répertoire perso dans le scratch
cd /scratch
mkdir arnal-$SLURM_JOB_ID

# récupérer les données depuis le san
scp -r san:/users/arnal/VERDI_RECOVER/bootstrap_iqtree /scratch/arnal-$SLURM_JOB_ID

module load singularity


ALIGNMENT="no_outliers_aligned_verdiGISAIDref_seq.fasta"
DATES="no_outliers_treetime_verdiGISAIDref_df.csv"
UFBOOT="iqtree.ufboot"
OUTDIR="treetime_bootstrap_runs"

cd /scratch/arnal-$SLURM_JOB_ID/bootstrap_iqtree 


mkdir -p "$OUTDIR"

echo "Splitting bootstrap trees..."

# Split .ufboot into individual tree files
split -l 1 --numeric-suffixes=1 --suffix-length=4 "$UFBOOT" "$OUTDIR/tree_"
cd "$OUTDIR"

# Loop through each tree file
for TREEFILE in tree_*; do
    BOOTNUM=$(basename "$TREEFILE" | sed 's/tree_//')

    RUNDIR="run_${BOOTNUM}"
    mkdir -p "$RUNDIR"
    mv "$TREEFILE" "$RUNDIR/tree.nwk"

    echo "Running TreeTime on bootstrap replicate $BOOTNUM"
    cd "$RUNDIR"

    singularity run https://depot.galaxyproject.org/singularity/treetime:0.11.4--pyhdfd78af_0 \
        treetime --tree tree.nwk \
                --aln ../../$ALIGNMENT \
                --dates ../../$DATES \
                --outdir output \
                --reroot NC_045512.2  \
                --stochastic-resolve
                

    cd ..
done

echo "All TreeTime runs completed."
scp -r /scratch/arnal-$SLURM_JOB_ID/treetime_bootstrap/treetime_bootstrap_runs  san:/users/arnal/VERDI_RECOVER/treetime_array



cd /scratch
rm -rf arnal-$SLURM_JOB_ID