#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -t 48:00:00
#SBATCH -J run_salmon_mapping_based_consistent
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=32G

# Load modules
module load Salmon/1.10.3-GCC-13.3.0 #Pelle has different version than ive used before.

# Paths
INDEX_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_index_consistent"
DATA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/salmon_mapping_quant_consistent"

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
