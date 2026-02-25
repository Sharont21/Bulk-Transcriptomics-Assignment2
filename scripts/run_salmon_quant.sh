#!/bin/bash
#SBATCH --job-name=salmon_quant
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=slurm-salmon-quant-%j.out
#SBATCH --error=slurm-salmon-quant-%j.err

module load StdEnv/2023
module load salmon

FASTQ_DIR=$HOME/Bulk-Transcriptomics-Assignment2/data/fastq
INDEX_DIR=$HOME/Bulk-Transcriptomics-Assignment2/data/reference/salmon_index
OUT_DIR=$HOME/Bulk-Transcriptomics-Assignment2/data/salmon

# List of samples
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

for SRR in "${SRR_LIST[@]}"; do
  echo "Quantifying $SRR"
  salmon quant \
    -i $INDEX_DIR \
    -l U \
    -r $FASTQ_DIR/${SRR}_1.fastq \
    -p 4 \
    --validateMappings \
    -o $OUT_DIR/${SRR}_quant
done
