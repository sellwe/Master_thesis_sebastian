#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -J create_transcript_unfiltered
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

# Load required modules
module load bioinfo-tools
module load gffread/0.12.7

# Define paths
WORKDIR=/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac
OUTDIR=$WORKDIR

cd $WORKDIR

echo "Running gffread on $(date)"
echo "Working directory: $WORKDIR"

# Input files
GFF=braker_unfiltered.gff3
FASTA=C_maculatus_assembly.fna #only use the hard-masked where it is required. Keep the generated assembly until then.
OUT=transcripts_unfiltered.fa

# Run gffread
gffread $GFF -g $FASTA -w $OUT

echo "Job finished on $(date)"



