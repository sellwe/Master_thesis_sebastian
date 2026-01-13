#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 08:00:00
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

# transcript level with multimapping
# -M and --fraction for multimapping
# -f for transcript-level 

featureCounts -T 16 \
  -a "$ANNOTATION" \
  -o "$OUT_DIR/transcript_counts_multimappers.txt" \
  -p -B -C \
  -M \
  --fraction \
  -f \
  -g "transcript_id" \
  -t "exon" \
  -s 2 \
  "$INPUT_DIR"/*_marked_duplicates.bam

