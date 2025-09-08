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

set -e


# Créer un répertoire perso dans le scratch
cd /scratch
mkdir arnal-$SLURM_JOB_ID

# récupérer les données depuis le san
scp -r san:/users/arnal/VERDI_RECOVER/treetime_bootstrap /scratch/arnal-$SLURM_JOB_ID

module load iqtree
module load singularity


ALIGNMENT="no_outliers_aligned_verdiGISAIDref_seq.fasta"
DATES="no_outliers_treetime_verdiGISAIDref_df.csv"


cd /scratch/arnal-$SLURM_JOB_ID/treetime_bootstrap


# unfiltered sequences 
mkdir iqtree
iqtree2 -s $ALIGNMENT -m GTR+G -B 1000 -nt 16 -o "NC_045512.2"  -pre iqtree/iqtree --wbtl

echo "IQtree done"

scp -r /scratch/arnal-$SLURM_JOB_ID/treetime_bootstrap/iqtree  san:/users/arnal/VERDI_RECOVER/bootstrap_iqtree

echo "IQtree transferred"


# UFBOOT="iqtree/iqtree.ufboot"
# OUTDIR="treetime_bootstrap_runs"

# mkdir -p "$OUTDIR"

# echo "Splitting bootstrap trees..."

# # Split .ufboot into individual tree files
# split -l 1 "$UFBOOT" "$OUTDIR/tree_"
# cd "$OUTDIR"

# for TREEFILE in tree_*; do
#     BASENAME=$(basename "$TREEFILE")
#     BOOTNUM="${BASENAME#tree_}"
#     BOOTNUM="${BOOTNUM##*_}"

#     mkdir -p "run_${BOOTNUM}"
#     cp "../$ALIGNMENT" "run_${BOOTNUM}/aligned.fasta"
#     cp "../$DATES" "run_${BOOTNUM}/dates.tsv"
#     mv "$TREEFILE" "run_${BOOTNUM}/tree.nwk"

#     echo "Running TreeTime on bootstrap replicate $BOOTNUM"
#     cd "run_${BOOTNUM}"

#     singularity run https://depot.galaxyproject.org/singularity/treetime:0.11.4--pyhdfd78af_0 \
#         treetime --tree tree.nwk \
#                 --aln aligned.fasta \
#                 --dates dates.tsv \
#                 --outdir output \
#                 --reroot NC_045512.2  \
#                 --stochastic-resolve \
#                 --threads 16

#     cd ..
# done

# echo "All TreeTime runs completed."
# scp -r /scratch/arnal-$SLURM_JOB_ID/treetime_bootstrap/treetime_bootstrap_runs  san:/users/arnal/VERDI_RECOVER/treetime_bootstrap_runs



cd /scratch
rm -rf arnal-$SLURM_JOB_ID