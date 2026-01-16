#!/bin/bash -l
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 20
#SBATCH -t 12:00:00
#SBATCH -J salmon_index_whole_genomeunfiltered
#SBATCH --output=%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sebastian.ellwe.1520@student.uu.se

# Load modules 
module load bioinfo-tools
module load Salmon/1.10.1

GENTROME="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered/gentrome_whole_genome_unfiltered.fa"
DECOYS="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_decoy_whole_genome_unfiltered/decoys_whole_genome.txt"
INDEXDIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/salmon_index_whole_genome_unfiltered"

mkdir -p "$INDEXDIR"

salmon index -t "$GENTROME" -i "$INDEXDIR" -d "$DECOYS" -p 20 --kmerLen 31 --gencode


