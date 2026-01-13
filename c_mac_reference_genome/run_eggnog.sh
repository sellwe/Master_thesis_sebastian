#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 48:00:00
#SBATCH -J run_eggnog
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se
#SBATCH --output=%x.%j.out

module load bioinfo-tools 
module load eggNOG-mapper/2.1.9

INPUT="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/braker_proteins.fa"
OUTDIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/eggnog_out"
mkdir -p "$OUTDIR"

emapper.py \
  -i "$INPUT" \
  --itype proteins \
  -o Cmaculatus_functional_annotation \
  --output_dir "$OUTDIR" \
  --cpu 16
