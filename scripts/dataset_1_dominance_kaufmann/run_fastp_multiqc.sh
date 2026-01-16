#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J run_fastp_multiqc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules 

module load bioinfo-tools
module load fastp/0.23.4
module load MultiQC/1.22.2

# Paths 

DATA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/raw_data/fastq"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"


# Run fastp

echo "Running fastp trimming... "
for R1 in "$DATA_DIR"/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    R2="${DATA_DIR}/${SAMPLE}_2.fastq.gz"

    # Output files
    OUT_R1="${OUT_DIR}/${SAMPLE}_trimmed_1.fastq.gz"
    OUT_R2="${OUT_DIR}/${SAMPLE}_trimmed_2.fastq.gz"

    HTML="${OUT_DIR}/${SAMPLE}_fastp.html"
    JSON="${OUT_DIR}/${SAMPLE}_fastp.json"

    echo "Processing sample: $SAMPLE"

    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$OUT_R1" \
        -O "$OUT_R2" \
        --detect_adapter_for_pe \
        --trim_poly_g \
        --trim_poly_x \
        --thread 4 \
        --html "$HTML" \
        --json "$JSON"

done

echo "fastp trimming complete!"
echo "Running MultiQC..."

multiqc "$OUT_DIR" -o "$OUT_DIR"

echo "All steps completed successfully!"
