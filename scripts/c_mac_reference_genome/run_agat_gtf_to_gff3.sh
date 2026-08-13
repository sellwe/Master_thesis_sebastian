#!/bin/bash
#SBATCH -A uppmax2025-2-148 #change to the current UPPMAX settings (done before pelle)
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 12:00:00
#SBATCH -J run_AGAT
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se #change to your email 
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
   -g "$GTF_IN" \
   -o "$GFF_OUT"

agat_sp_fix_features_locations_duplicated.pl \
  -g "$GFF_OUT" \
  -o "$GFF_FIXED"

echo "AGAT conversion finished at $(date)"

