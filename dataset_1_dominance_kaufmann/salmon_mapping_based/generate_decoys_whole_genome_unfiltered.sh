#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 24:00:00
#SBATCH -J generate_decoys_whole_genome_unfiltered
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out
#SBATCH --mem=64G

# paths 

WORK_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered"

GENOME="${WORK_DIR}/C_maculatus_assembly.fna"
GTF="${WORK_DIR}/C_maculatus_annotation_unfiltered_fixed.gtf"
TRANSCRIPTS="${WORK_DIR}/transcripts_unfiltered.fa"

mkdir -p "$OUT_DIR"


grep "^>" "$GENOME" | cut -d " " -f 1 | sed 's/>//' > "$OUT_DIR/decoys_whole_genome.txt"

cat "$TRANSCRIPTS" "$GENOME" > "$OUT_DIR/gentrome_whole_genome_unfiltered.fa"

echo "decoys.txt and gentrome.fa created in $OUT_DIR"
