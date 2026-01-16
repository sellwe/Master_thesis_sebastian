#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 4-00:00:00
#SBATCH -J star_alignment_continuation_1
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load star/2.7.11a

# No indexing, only alignment
GENOME_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/star_index"
RNA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_gene_counts"

# Process first batch: ERR12383258 to ERR12383287 (30 samples)
START_SAMPLE="ERR12383258_trimmed"
END_SAMPLE="ERR12383287_trimmed"

start_processing=false

for R1 in ${RNA_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" "_1.fastq.gz")
    R2="${RNA_DIR}/${SAMPLE}_2.fastq.gz"

    # Only start once we reach start sample
    if [[ "$SAMPLE" == "$START_SAMPLE" ]]; then
        start_processing=true
    fi

    if ! $start_processing ; then
        continue
    fi

    # Stop when we reach end sample
    if [[ "$SAMPLE" == "$END_SAMPLE" ]]; then
        echo "Reached end sample $END_SAMPLE"
        start_processing=false
    fi

    echo "Checking $SAMPLE ..."

    BAM="${OUT_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
    LOG="${OUT_DIR}/${SAMPLE}_Log.final.out"
    GENES="${OUT_DIR}/${SAMPLE}_ReadsPerGene.out.tab"

    # Skip samples that already completed successfully
    if [[ -s "$BAM" && -f "$LOG" && -f "$GENES" ]]; then
        echo "Skipping $SAMPLE (already complete)"
        continue
    fi

    echo "=== STARTING $SAMPLE ==="

    # Align reads and produce gene counts
    STAR --runThreadN 16 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --quantMode GeneCounts \
         --twopassMode Basic \
         --outFilterMultimapNmax 20 \
         --limitBAMsortRAM 20000000000 \
         --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_"

    echo "=== FINISHED $SAMPLE ==="
done

echo "Batch 1 completed processing from $START_SAMPLE to $END_SAMPLE"
