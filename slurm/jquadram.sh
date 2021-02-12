#!/bin/sh
 
#SBATCH -p low
#SBATCH --output=jquadram.stdout
#SBATCH --error=jquadram.stderr
#SBATCH --job-name=jquadram
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8G
#SBATCH --exclude=gpunode001,gpunode002


now=$(date +"%T")
/home/nicola.vitulo/GIOELE/ovariancancer/slurm/tg_message.sh "START $SLURM_JOB_NAME @ $now"

#source /home/nicola.vitulo/GIOELE/miniconda3/bin/activate
#conda activate nf
cd /home/nicola.vitulo/GIOELE/MetaPhage
 
nextflow run main.nf -profile base,quadramvm

now=$(date +"%T")
/home/nicola.vitulo/GIOELE/ovariancancer/slurm/tg_message.sh "END $SLURM_JOB_NAME @ $now"