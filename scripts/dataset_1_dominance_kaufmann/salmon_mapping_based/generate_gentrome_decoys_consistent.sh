#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 06:00:00
#SBATCH -J generate_gentrome_decoys_consistent
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=32G

# on Pelle

# paths

WORK_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered"

GENOME="${WORK_DIR}/C_maculatus_assembly.fna"
GTF="${WORK_DIR}/C_maculatus_annotation_unfiltered_fixed.gtf"
TRANSCRIPTS="${WORK_DIR}/transcripts_unfiltered_consistent.fa"

mkdir -p "$OUT_DIR"

# extract the contig names as decoys (still the same)
grep "^>" "$GENOME" | cut -d " " -f 1 | sed 's/>//' > "$OUT_DIR/decoys_whole_genome_consistent.txt"

# build the gentrome (based on new transcriptome)
cat "$TRANSCRIPTS" "$GENOME" > "$OUT_DIR/gentrome_whole_genome_unfiltered_consistent.fa"

