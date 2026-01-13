#!/bin/bash
#SBATCH -A uppmax2025-2-148 
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 48:00:00
#SBATCH -J download_PRJEB70958
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load sratools/3.0.7
module load pigz

# Define directories
BASE=/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/
SRA=$BASE/raw_data/sra
FASTQ=$BASE/raw_data/fastq
META=$BASE/metadata

mkdir -p $SRA $FASTQ

# Download and convert all runs
for run in $(cat $META/SRR_Acc_List_PRJEB70958.txt); do
    prefetch $run --output-directory $SRA					# Download .sra file
    fasterq-dump $SRA/$run/$run.sra -O $FASTQ --split-files --threads 8		# Convert to FASTQ (paired-end, gzipped)
done

pigz -p 8 $FASTQ/*.fastq						# Compress
