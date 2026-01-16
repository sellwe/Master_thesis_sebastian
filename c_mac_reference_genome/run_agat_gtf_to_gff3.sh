#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J run_AGAT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

# Load modules
module load bioinfo-tools
module load AGAT

# Input / output
GTF_IN="C_maculatus_annotation_unfiltered_fixed.gtf"
GFF_OUT="C_maculatus_annotation_unfiltered.gff3"
GFF_FIXED="C_maculatus_annotation_unfiltered.fixed.gff3"

echo "Starting AGAT conversion at $(date)"

# Convert GTF to GFF3 
agat_convert_sp_gxf2gxf.pl \
   -g C_maculatus_annotation_unfiltered_fixed.gtf \
   -o C_maculatus_annotation_unfiltered.gff3

# Fix duplicated features / ordering 
agat_sp_fix_features_locations_duplicated.pl \
  -g C_maculatus_annotation_unfiltered.gff3 \
  -o C_maculatus_annotation_unfiltered_fixed.gff3

echo "AGAT conversion finished at $(date)"

