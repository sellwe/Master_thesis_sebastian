#!/bin/bash
#SBATCH -A uppmax2026-1-8
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH -t 12:00:00
#SBATCH -J salmon_index_consistent
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=32G

# on Pelle

#Load modules
module load Salmon/1.10.3-GCC-13.3.0 #Pelle has different version than ive used before.

GENTROME="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered/gentrome_whole_genome_unfiltered_consistent.fa"
DECOYS="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered/decoys_whole_genome_consistent.txt"
INDEXDIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_index_consistent"

mkdir -p "$INDEXDIR"

salmon index -t "$GENTROME" -i "$INDEXDIR" -d "$DECOYS" -p 20 --kmerLen 31 --gencode

