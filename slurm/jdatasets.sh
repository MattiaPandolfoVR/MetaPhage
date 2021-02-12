#!/bin/sh
 
#SBATCH -p low
#SBATCH --output=jdatasets.stdout
#SBATCH --error=jdatasets.stderr
#SBATCH --job-name=jdatasets
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --exclude=gpunode001,gpunode002


now=$(date +"%T")
/home/nicola.vitulo/GIOELE/ovariancancer/slurm/tg_message.sh "START $SLURM_JOB_NAME @ $now"

cd /home/nicola.vitulo/GIOELE/MetaPhage/datasets/base
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/003/SRR8299383/SRR8299383_1.fastq.gz && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/003/SRR8299383/SRR8299383_2.fastq.gz && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/004/SRR8299384/SRR8299384_1.fastq.gz && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/004/SRR8299384/SRR8299384_2.fastq.gz && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/006/SRR8299386/SRR8299386_1.fastq.gz && wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR829/006/SRR8299386/SRR8299386_2.fastq.gz

now=$(date +"%T")
/home/nicola.vitulo/GIOELE/ovariancancer/slurm/tg_message.sh "END $SLURM_JOB_NAME @ $now"