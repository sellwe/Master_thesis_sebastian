#!/bin/bash
#SBATCH -A uppmax2026-1-8 #change to the current UPPMAX settings (done before pelle)
#SBATCH -p pelle
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 06:00:00
#SBATCH -J parse_gene_tree_birth_nodes_2 
#SBATCH --mail-type=ALL
#SBATCH --output=%x.%j.out
#SBATCH --mem=16G

# Load Biopython (also loads Python 3 as a dependency)
module load Biopython/1.84-gfbf-2024a

# Paths
SCRIPT_DIR="/proj/naiss2023-6-65/Sebastian/Master_thesis_sebastian/scripts/c_mac_reference_genome"

echo "Starting gene tree birth node parsing: $(date)"
python3 ${SCRIPT_DIR}/parse_gene_trees_birth_nodes_2.py
echo "Finished: $(date)"
