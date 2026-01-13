#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 48:00:00
#SBATCH -J salmon_mapping_whole_genome_unfiltered
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load Salmon/1.10.1

# Paths
INDEX_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_index_whole_genome_unfiltered"
DATA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/salmon_quantification_whole_genome_unfiltered"

mkdir -p "$OUT_DIR"

# Loop over samples
for R1 in "$DATA_DIR"/*_1.fastq.gz; do
    R2="${R1/_1.fastq.gz/_2.fastq.gz}"
    SAMPLE_NAME=$(basename "$R1" "_1.fastq.gz")
    echo "Quantifying sample: $SAMPLE_NAME"

    salmon quant \
        -i "$INDEX_DIR" \
        -l A \
        -1 "$R1" \
        -2 "$R2" \
        -p 20 \
        --gcBias \
        --seqBias \
        --validateMappings \
        -o "$OUT_DIR/$SAMPLE_NAME"
done

echo "All quantifications completed successfully!"
