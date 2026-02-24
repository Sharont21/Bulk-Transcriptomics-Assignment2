#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=slurm-fastqc-%j.out
#SBATCH --error=slurm-fastqc-%j.err

module load StdEnv/2023
module load fastqc

FASTQ_DIR=$HOME/Bulk-Transcriptomics-Assignment2/data/fastq
OUT_DIR=$HOME/Bulk-Transcriptomics-Assignment2/results/fastqc

fastqc \
  -t 4 \
  -o $OUT_DIR \
  $FASTQ_DIR/*.fastq
