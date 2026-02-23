#!/bin/bash
#SBATCH --job-name=sra_download
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm-sra-%j.out
#SBATCH --error=slurm-sra-%j.err

module load StdEnv/2023
module load sra-toolkit

cd $HOME/Bulk-Transcriptomics-Assignment2/data/fastq

# List of SRA accessions
SRR_LIST=(
  SRR10551657
  SRR10551658
  SRR10551659
  SRR10551660
  SRR10551661
  SRR10551662
  SRR10551663
  SRR10551664
  SRR10551665
)

# Download FASTQ files
for SRR in "${SRR_LIST[@]}"; do
  echo "Downloading $SRR"
  fasterq-dump $SRR --split-files --threads 4
done
