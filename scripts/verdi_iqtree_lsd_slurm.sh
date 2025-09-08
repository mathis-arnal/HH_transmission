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
scp -r san:/users/arnal/VERDI_RECOVER/IQTREE_verdi /scratch/arnal-$SLURM_JOB_ID

module load iqtree
module load mafft
module load R/4.4.3
module load singularity

cd /scratch/arnal-$SLURM_JOB_ID/IQTREE_verdi

# ------------------ VERDI ONLY SEQUENCES -------------
mafft --keeplength --add verdi_seq.fasta --reorder NC_045512.2.fasta > aligned_verdiref_seq.fasta


# unfiltered sequences 
mkdir verdi_iqtree
iqtree2 -s aligned_verdiref_seq.fasta -m GTR+G -B 1000 -nt 16 -o "NC_045512.2"  -pre verdi_iqtree/lsd_iqtree --date verdiGISAIDref_df.txt --date-ci 100

mkdir verdi_treecluster
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.000084 -m max_clade  > verdi_treecluster/2.5SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.000034 -m max_clade  > verdi_treecluster/1SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.0002 -m max_clade  > verdi_treecluster/6SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.5 -m max_clade  > verdi_treecluster/0.5_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i verdi_iqtree/lsd_iqtree.timetree.nwk -t 1 -m max_clade  > verdi_treecluster/1_treecluster.txt

# filtered sequences
# filter the sequences
Rscript filtering_seq.R


mkdir filtered_verdi_iqtree

iqtree2 -s filtered_aligned_verdiref_seq.fasta -m GTR+G -B 1000 -nt 16 -o "NC_045512.2"  -pre filtered_verdi_iqtree/lsd_iqtree --date verdiGISAIDref_df.txt --date-ci 100

mkdir verdi_filtered_treecluster
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i filtered_verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.000084 -m max_clade  > verdi_filtered_treecluster/2.5SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i filtered_verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.000034 -m max_clade  > verdi_filtered_treecluster/1SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i filtered_verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.0002 -m max_clade  > verdi_filtered_treecluster/6SNPs_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i filtered_verdi_iqtree/lsd_iqtree.timetree.nwk -t 0.5 -m max_clade  > verdi_filtered_treecluster/0.5_treecluster.txt
singularity exec https://depot.galaxyproject.org/singularity/treecluster:1.0.4--pyh7cba7a3_0 TreeCluster.py -i filtered_verdi_iqtree/lsd_iqtree.timetree.nwk -t 1 -m max_clade  > verdi_filtered_treecluster/1_treecluster.txt

# Copy the result
scp -r /scratch/arnal-$SLURM_JOB_ID/IQTREE_verdi  san:/users/arnal/VERDI_RECOVER/IQTREE_res_verdi

# Supress the scratch   
cd /scratch
rm -rf arnal-$SLURM_JOB_ID
