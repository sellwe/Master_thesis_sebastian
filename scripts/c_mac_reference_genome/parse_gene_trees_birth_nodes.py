#!/usr/bin/env python3
"""
parse_gene_trees_birth_nodes.py

Assigns every gene in every OrthoFinder orthogroup a "birth node":
the point in evolutionary history where that specific gene copy first
became a distinct lineage.

The output is one row per gene per species. The output can be filtered
for the species of interest downstream. 

Birth node of a gene can be assigned in three ways:

1 duplication: 
    Using OrthoFinder's pre-reconciled, rooted gene trees.
    Walk from the leaf toward the root in the gene tree. The first
    duplication node encountered is the birth node. That node is mapped
    to a species tree node via Duplications.tsv, giving evolutionary age.

2 mrca_inferred:
   The gene was never duplicated in any lineage within this orthogroup,
   meaning it has been conserved as a single copy. Its age is estimated
   as the MRCA of all species that share this orthogroup: the oldest
   point in evolutionary history where the gene is documented to exist.
   Gene loss in intermediate lineages could make
   the gene appear older than it is, but it is biologically meaningful
   than assigning everything to rank 1. 

3 species_specific:
   No orthologs detected in any other species. The gene exists only in
   this species' orthogroup. Age is unknown (could be a young lineage-
   specific gene, a highly diverged gene that lost detectable homology,
   or an annotation artifact). Kept as NA for downstream handling.

In R: species_tree_birth_node will be merged with node_age_table to get age_rank.
birth_type shows how each age was determined.

Output columns:
  transcript_id          : e.g. g6090.t1 (matches annotation)
  gene_id                : e.g. g6090 (matches annotation)
  gene_id_raw            : full OrthoFinder gene ID
  species                : species prefix
  orthogroup             : OG identifier
  gene_tree_birth_node   : within gene tree duplication node (e.g. n7), or NA
  species_tree_birth_node: species tree node (e.g. N1, N5), or NA
  birth_type             : duplication | mrca_inferred | species_specific
  n_copies_in_og         : copies this species has in the OG
  n_species_in_og        : distinct species in the OG
"""

import os
import csv
import io
from collections import defaultdict
from Bio import Phylo

# -- Paths ----------------------------------------------------------------------
RESULTS_DIR  = "/proj/naiss2023-6-65/Sebastian/data/orthofinder/Results_May27"
TREE_DIR     = os.path.join(RESULTS_DIR, "Resolved_Gene_Trees")
DUP_FILE     = os.path.join(RESULTS_DIR, "Gene_Duplication_Events", "Duplications.tsv")
SP_TREE_FILE = os.path.join(RESULTS_DIR, "Species_Tree",
                             "SpeciesTree_rooted_node_labels.txt")
OUT_FILE     = ("/proj/naiss2023-6-65/Sebastian/Master_thesis_sebastian"
                "/scripts/c_mac_reference_genome/transcript_birth_nodes.tsv")

SEP = "_filtered_proteinfasta_TE_filtered_"


# -- Helpers ----------------------------------------------------------------------
def get_species(gene_id):
    """
    Extracts the full species prefix from an OrthoFinder gene ID.
    C_maculatus_filtered_proteinfasta_TE_filtered_C_maculatus_g6090.t1_1
    -> C_maculatus_filtered_proteinfasta_TE_filtered
    The species prefix is the same as the tip label in species_tree
    """
    idx = gene_id.find(SEP)
    if idx < 0:
        return "unknown"
    return gene_id[: idx + len(SEP) - 1]


def clean_transcript_id(gene_id):
    """
    C_maculatus_filtered_proteinfasta_TE_filtered_C_maculatus_g6090.t1_1
    -> g6090.t1
    Mirrors what i did in the R annotation scripts.
    """
    idx = gene_id.find(SEP)
    if idx < 0:
        return gene_id
    remainder = gene_id[idx + len(SEP):]          # C_maculatus_g6090.t1_1

    # Strip species name (Genus_species_ prefix)
    parts = remainder.split('_')
    if len(parts) >= 3:
        short_prefix = parts[0] + '_' + parts[1] + '_'
        if remainder.startswith(short_prefix):
            remainder = remainder[len(short_prefix):]  # g6090.t1_1

    # Strip the trailing OrthoFinder protein index (_1, _2, ...)
    tail = remainder.rsplit('_', 1)
    if len(tail) == 2 and tail[1].isdigit():
        return tail[0]    # g6090.t1
    return remainder


def get_gene_id(transcript_id):
    """g6090.t1 -> g6090"""
    return transcript_id.split('.')[0]


# -- Load species tree -------------------------------------------------------------
# The species tree is needed for the MRCA nodes for unduplicated genes
# For duplicated genens we get what we need from Duplications.tsv
print("Loading species tree ...", flush=True)
sp_tree      = Phylo.read(SP_TREE_FILE, 'newick')
sp_tree_tips = {leaf.name for leaf in sp_tree.get_terminals()}
print(f"  Species found: {sorted(sp_tree_tips)}", flush=True)


def get_mrca_label(species_set):
    """
    Returns the node label of the MRCA of all species in species_set
    that exist in the species tree.

    For unduplicated genes this represents the maximum age of that gene:
    the oldest point in evolutionary history where it is documented to
    have existed (based on the most distantly related species that share it).

    Returns None if fewer than 2 valid species are found. Single-species OG is handled 
    leter as species_specific
    """
    # keep only species in the species tree 
    valid = [s for s in species_set if s in sp_tree_tips]
    if len(valid) < 2:
        return None

    # get the Bio.Phylo clade object for each species 
    target_clades = []
    for sp in valid:
        clade = next(sp_tree.find_clades(sp), None)
        if clade is not None:
            target_clades.append(clade)

    if len(target_clades) < 2:
        return None

    # Bio.Phylo common_ancestor finds the mrca of all 
    # provided clades in the species tree
    mrca  = sp_tree.common_ancestor(target_clades)
    label = mrca.name
    if label and label.strip():
        return label.strip()
    return None


# -- Load Duplications.tsv ---------------------------------------------------------
# Duplications.tsv contains all duplicaiton events. Each row is one duplication event:
# it records the orthogroup, the internal node in the gene tree where the 
# duplication occurred (Gene Tree Node), and the node in the species tree where 
# OrthoFinder places that duplication in evolutionary time (Species Tree Node)

# Build a two level-dictionary: 
# og_dup_lookup[orthogroup][gene_tree_node] = species_tree_node
# to correlate each gene tree node for a orthogroup and which species tree node it correlates to 

print("Loading Duplications.tsv ...", flush=True)
og_dup_lookup = defaultdict(dict)

with open(DUP_FILE) as fh:
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        og  = row['Orthogroup'].strip()
        gtn = row['Gene Tree Node'].strip()
        stn = row['Species Tree Node'].strip()
        if gtn not in og_dup_lookup[og]:
            og_dup_lookup[og][gtn] = stn

print(f"  Loaded duplication events for {len(og_dup_lookup):,} orthogroups.",
      flush=True)


# -- Process gene trees ------------------------------------------------------------
# For each orthogroup we:
# 1 parse the resolved gene tree (Newick file, one per OG)
# 2 For each gene (leaf) in the tree, find its birth node
# 3 record results 

# The resolved gene trees are rooted against the species tree and have had
# DLC (Duplication-Loss-Coalescence) reconcilliation done. The duplications node 
# labels match Duplications.tsv 

print("Processing resolved gene trees ...", flush=True)

tree_files = sorted(f for f in os.listdir(TREE_DIR) if f.endswith('_tree.txt'))
n_files    = len(tree_files)
out_rows   = []

# Counters for final summary
n_duplication    = 0
n_mrca_inferred  = 0
n_species_specific = 0

for i, fname in enumerate(tree_files, 1):
    if i % 2000 == 0:
        print(f"  {i:,}/{n_files:,} trees processed ...", flush=True)

    # Parse this orthogroup
    # Orthogroup ID minus _tree.txt
    og    = fname[:-9]                              
    fpath = os.path.join(TREE_DIR, fname)

    with open(fpath) as fh:
        tree_str = fh.read().strip()
    #some OG files are empty and will be skipped
    if not tree_str:
        continue

    # Use Bio.Phylo to parse the newick string into a tree object
    # Internal nodes are labeled n1, n2 etc
    # These labels become the .name atrribute of each interal clade object
    # Leaves are labeled with the full Orthofinder gene ID
    tree      = Phylo.read(io.StringIO(tree_str), 'newick')

    # Look up all known duplication nodes of this OG.
    # dup_nodes is a dict: {gene_tree_node_label: species_tree_node_label}
    # ex. {"n7": "N5", "n12": "C_maculatus_filtered_proteinfasta_TE_filtered"
    # if the OG has no duplications, this is empty 
    dup_nodes = og_dup_lookup.get(og, {})

    # Get all leaf nodes here (the genes themselves)
    terminals = list(tree.get_terminals())

    # Count species representation within this OG
    # n_cpies_in_og: the paralog count per species 
    # n_species_in_og + get_mrca_label: to decide whether an unduplicated
    # gene should get an mrca_inferred age or be flagged species_specific
    sp_count = defaultdict(int)
    for leaf in terminals:
        sp_count[get_species(leaf.name)] += 1

    n_species_in_og = len(sp_count)

    # Precompute MRCA for unduplicated genes in this OG (once per OG)
    # og_mrca is is the species tree node label for the MRCA of all species in
    # this OG. For example, if the OG contains C_maculatus and T_castaneum,
    # og_mrca = "N4". If only C_maculatus is in the OG, og_mrca = None.
    og_mrca = get_mrca_label(sp_count.keys())

    # assign birth nodes for each gene
    for leaf in terminals:
        gene_id_raw   = leaf.name
        species       = get_species(gene_id_raw)
        transcript_id = clean_transcript_id(gene_id_raw)
        gene_id       = get_gene_id(transcript_id)

        # Walk from immediate parent toward root.
        # Every node in get_path() is a direct ancestor of this leaf. 
        # Stop at the first duplication, which created this particualr copy.
        # First duplication node hit = birth node.
        path      = tree.get_path(leaf)
        birth_gtn = None # gene tree node label (ex n7)
        birth_stn = None #species tree node label (ex N5)
        birth_type = None

        for clade in reversed(path[:-1]):
            if clade.name and clade.name in dup_nodes:
                # found birth duplication
                # clade.name is the gene tree node
                # dup_nodes[clade.name] is the specie tree node 
                birth_gtn  = clade.name
                birth_stn  = dup_nodes[birth_gtn]
                birth_type = "duplication"
                break

        # for genes with no duplication ancestor
        if birth_type is None:
            if og_mrca is not None:
                # Conserved single-copy gene: assign MRCA age as upper bound
                # assign latest MRCA of species in this OG as the age.
                birth_stn  = og_mrca
                birth_type = "mrca_inferred"
            else:
                # No orthologs outside this species: age unknown
                birth_stn  = "NA"
                birth_type = "species_specific"

        # Update summary counters 
        if birth_type == "duplication":
            n_duplication += 1
        elif birth_type == "mrca_inferred":
            n_mrca_inferred += 1
        else:
            n_species_specific += 1

        out_rows.append({
            'transcript_id':           transcript_id,
            'gene_id':                 gene_id,
            'gene_id_raw':             gene_id_raw,
            'species':                 species,
            'orthogroup':              og,
            'gene_tree_birth_node':    birth_gtn or 'NA',   # node within OG
            'species_tree_birth_node': birth_stn or 'NA',   # age rank 
            'birth_type':              birth_type,
            'n_copies_in_og':          sp_count[species],   # nr of genes in OG in this species
            'n_species_in_og':         n_species_in_og,     # nr of species in OG. Used for mrca, 
        })

print(f"  Total gene entries:        {len(out_rows):,}")
print(f"  birth_type=duplication:    {n_duplication:,}")
print(f"  birth_type=mrca_inferred:  {n_mrca_inferred:,}")
print(f"  birth_type=species_specific: {n_species_specific:,}")


# -- Write output ------------------------------------------------------------------
print(f"Writing output to:\n  {OUT_FILE}", flush=True)

fieldnames = [
    'transcript_id',
    'gene_id',
    'gene_id_raw',
    'species',
    'orthogroup',
    'gene_tree_birth_node',
    'species_tree_birth_node',
    'birth_type',
    'n_copies_in_og',
    'n_species_in_og',
]

with open(OUT_FILE, 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter='\t')
    w.writeheader()
    w.writerows(out_rows)

print("Done.", flush=True)
