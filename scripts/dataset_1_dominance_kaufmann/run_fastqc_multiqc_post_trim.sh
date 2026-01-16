#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J run_fastqc_multiqc_post_trim
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules 

module load bioinfo-tools
module load FastQC/0.11.9
module load MultiQC/1.22.2

# Paths 

DATA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/post_trim_fastqc_multiqc"


# Run FastQC

echo "Running FastQC on all FASTQ files..."
fastqc -t 20 -o "$OUT_DIR" "$DATA_DIR"/*.fastq.gz

# Run MultiQC

echo "Running MultiQC to summarize results..."
multiqc "$OUT_DIR" -o "$OUT_DIR"

echo "All steps completed successfully!"

