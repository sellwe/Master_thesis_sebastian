#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 1-00:00:00
#SBATCH -J featureCounts
#SBATCH --output=featureCounts.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se

module load bioinfo-tools
module load subread/2.0.3

#Directories
ANNOTATION="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/C_maculatus_annotation_unfiltered_fixed.gtf"
INPUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_gene_counts/picard_marked_indexed"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_gene_counts/featureCounts"

mkdir -p "$OUT_DIR"

# Mode 1: Standard counting (unique mappers only)
# No multimapping, at the gene level, used for differential expression analysis
echo "=== Mode 1: Standard (unique mappers only) ==="
featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/gene_counts_standard.txt" \
  -p -B -C \
  -g "gene_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam

# Mode 2: With multimappers, still on gene level (fractional counting, splitting reads between the multipe targets)
echo "=== Mode 2: With multimappers (fractional) ==="
featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/gene_counts_multimappers.txt" \
  -p -B -C \
  -M \
  --fraction \
  -g "gene_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam

# Mode 3: Transcript-level counting (seems like the paper used this)
echo "=== Mode 3: Transcript-level counting ==="
featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/transcript_counts.txt" \
  -p -B -C \
  -f \
  -g "transcript_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam

echo "=== All counting modes completed ==="
echo "Output files:"
echo "1. $OUT_DIR/gene_counts_standard.txt"
echo "2. $OUT_DIR/gene_counts_multimappers.txt"
echo "3. $OUT_DIR/transcript_counts.txt"
