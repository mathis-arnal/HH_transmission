#!/bin/sh

########## configuration Slurm #######
## Définir un nom de job
#SBATCH --job-name=beast
## Définir la partition
#SBATCH -p normal
## définir notre adresse-mail
#SBATCH --mail-user=mathis.arnal@ird.fr
##  Définir qd recevoir le mail(BEGIN, ABORT, END ou ALL)
#SBATCH  --mail-type=ALL
#######################################


# Créer un répertoir perso dans le scratch
cd /scratch
mkdir arnal-$SLURM_JOB_ID

# récupérer les données depuis le san
scp -r san:/users/arnal/VERDI_RECOVER/BEAST /scratch/arnal-$SLURM_JOB_ID

module load singularity

cd /scratch/arnal-$SLURM_JOB_ID/BEAST

singularity exec https://depot.galaxyproject.org/singularity/beast:1.10.4--hdfd78af_2 beast -overwrite 10.4.0_HKY_coalescent_exponential.xml


scp -r /scratch/arnal-$SLURM_JOB_ID/BEAST  san:/users/arnal/VERDI_RECOVER/beast_res

# supression du scratch
cd /scratch
rm -rf arnal-$SLURM_JOB_ID