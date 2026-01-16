#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 2-12:00:00
#SBATCH -J star_transcriptome_3
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load star/2.7.11a

# Paths
GENOME_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/star_index"
RNA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_transcriptome_bam"

mkdir -p "$OUT_DIR"

# JOB 3: ERR12383294 to ERR12383317
START_SAMPLE="ERR12383294_trimmed"
END_SAMPLE="ERR12383317_trimmed"

echo "=== JOB 3: $START_SAMPLE to $END_SAMPLE (24 samples) ==="

start_processing=false

for R1 in ${RNA_DIR}/*_1.fastq.gz; do
    SAMPLE=$(basename "$R1" "_1.fastq.gz")
    R2="${RNA_DIR}/${SAMPLE}_2.fastq.gz"

    if [[ "$SAMPLE" == "$START_SAMPLE" ]]; then
        start_processing=true
    fi

    if ! $start_processing ; then
        continue
    fi

    # Check if already completed
    TRANSCRIPT_BAM="${OUT_DIR}/${SAMPLE}_Aligned.toTranscriptome.out.bam"
    LOG="${OUT_DIR}/${SAMPLE}_Log.final.out"
    
    if [[ -s "$TRANSCRIPT_BAM" && -f "$LOG" ]]; then
        echo "Skipping $SAMPLE (transcriptome BAM and log already exist)"
        if [[ "$SAMPLE" == "$END_SAMPLE" ]]; then
            break
        fi
        continue
    fi

    echo "=== PROCESSING $SAMPLE ==="

    STAR --runThreadN 16 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --readFilesCommand zcat \
         --twopassMode Basic \
         --quantMode TranscriptomeSAM \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes NH HI AS nM XS GX GN \
         --outSAMmapqUnique 60 \
         --alignEndsType EndToEnd \
         --outFilterMultimapNmax 20 \
         --winAnchorMultimapNmax 100 \
         --limitBAMsortRAM 20000000000 \
         --outFileNamePrefix "${OUT_DIR}/${SAMPLE}_"

    if [[ ! -f "$TRANSCRIPT_BAM" ]]; then
        echo "ERROR: Failed for $SAMPLE"
        exit 1
    fi

    if [[ "$SAMPLE" == "$END_SAMPLE" ]]; then
        break
    fi
done

echo "JOB 3 COMPLETED"
