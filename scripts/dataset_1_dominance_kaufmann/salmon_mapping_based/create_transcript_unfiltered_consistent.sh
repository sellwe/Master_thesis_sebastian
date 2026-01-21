#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 06:00:00
#SBATCH -J create_transcript_unfiltered_consistent
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=32G

# on Pelle

# Load required modules
module load gffread/0.12.7-GCCcore-13.3.0

# Define paths
WORKDIR=/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac
OUTDIR=$WORKDIR

cd $WORKDIR

echo "Running gffread on $(date)"
echo "Working directory: $WORKDIR"

# Input files
GTF=C_maculatus_annotation_unfiltered_fixed.gtf
FASTA=C_maculatus_assembly.fna #only use the hard-masked where it is required. Keep the generated assembly until then.
OUT=transcripts_unfiltered_consistent.fa

# Run gffread
gffread $GTF -g $FASTA -w $OUT

echo "Job finished on $(date)"

