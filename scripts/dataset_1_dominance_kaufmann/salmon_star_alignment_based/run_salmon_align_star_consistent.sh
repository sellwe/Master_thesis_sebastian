#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -t 48:00:00
#SBATCH -J salmon_alignment_based_consistent
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=32G

module load Salmon/1.10.3-GCC-13.3.0 #Pelle has different version than ive used before.

# Paths
TRANSCRIPTS="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/transcripts_unfiltered_consistent.fa"
STAR_TRANS_BAM_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_transcriptome_bam"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/salmon_alignment_based_quant_consistent"

mkdir -p "$OUT_DIR"

for BAM in "$STAR_TRANS_BAM_DIR"/*_Aligned.toTranscriptome.out.bam; do
    SAMPLE_NAME=$(basename "$BAM" "_Aligned.toTranscriptome.out.bam")
    echo "=== Quantifying sample: $SAMPLE_NAME ==="

# validateMappings is not used in this version/for alignment based mode
    salmon quant \
        --libType A \
        --alignments "$BAM" \
        --targets "$TRANSCRIPTS" \
        -p 20 \
        --gcBias \
        --seqBias \
        -o "$OUT_DIR/$SAMPLE_NAME"

    if [[ $? -ne 0 ]]; then
        echo "ERROR: Salmon failed for $SAMPLE_NAME"
        exit 1
    fi
done

echo "=== All alignment-based quantifications completed successfully ==="
