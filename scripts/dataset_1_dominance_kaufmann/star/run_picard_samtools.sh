#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 4:00:00
#SBATCH -J picard_markdup
#SBATCH --output=%x.%A_%a.out  #adds the array number after the jobid
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --array=1-72

# Load modules
module load bioinfo-tools
module load picard/3.1.1
module load samtools/1.20

# Directories
INPUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_gene_counts"
OUTPUT_DIR="$INPUT_DIR/picard_marked_indexed"

mkdir -p "$OUTPUT_DIR"

#Make an array of all the .bam files 
BAM_FILES=("$INPUT_DIR"/*_Aligned.sortedByCoord.out.bam)

#pick out the specific bam file for each array job 
BAM="${BAM_FILES[$SLURM_ARRAY_TASK_ID - 1]}"  #arrays are 0-based

SAMPLE=$(basename "$BAM" "_Aligned.sortedByCoord.out.bam")
OUTPUT_BAM="$OUTPUT_DIR/${SAMPLE}_marked_duplicates.bam"
METRICS="$OUTPUT_DIR/${SAMPLE}_markdup_metrics.txt"

echo "Array task $SLURM_ARRAY_TASK_ID processing $SAMPLE..."

# Mark duplicates
java -jar $PICARD_HOME/picard.jar MarkDuplicates \
    INPUT="$BAM" \
    OUTPUT="$OUTPUT_BAM" \
    METRICS_FILE="$METRICS" \
    VALIDATION_STRINGENCY=LENIENT

# Index the marked BAM
samtools index "$OUTPUT_BAM"
