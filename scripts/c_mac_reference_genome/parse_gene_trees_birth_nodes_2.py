#!/usr/bin/env python3
"""
parse_gene_trees_birth_nodes_2.py

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

    NOTE on root handling: Bio.Phylo's get_path() excludes the root node
    and includes the leaf. A gene sitting directly under the root would
    therefore miss the root duplication in the main loop. The root is
    checked explicitly AFTER the main loop, and only if the main loop
    found nothing. This avoids incorrectly overwriting a valid birth node
    found deeper in the tree with an older root duplication.

2 mrca_inferred:
    The gene was never duplicated in any lineage within this orthogroup,
    meaning it has been conserved as a single copy. Its age is estimated
    as the MRCA of all species that share this orthogroup: the oldest
    point in evolutionary history where the gene is documented to exist.
    Gene loss in intermediate lineages could make the gene appear older
    than it is, but it is biologically more meaningful than assigning
    everything to rank 1.

3 species_specific:
    No orthologs detected in any other species. The gene exists only in
    this species' orthogroup. Age is unknown (could be a young lineage-
    specific gene, a highly diverged gene that lost detectable homology,
    or an annotation artifact). Kept as NA for downstream handling.

In R: species_tree_birth_node will be merged with node_age_table to get
age_rank. birth_type shows how each age was determined.

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

# -- Paths ---------------------------------------------------------------------
RESULTS_DIR  = "/proj/naiss2023-6-65/Sebastian/data/orthofinder/Results_May27"
TREE_DIR     = os.path.join(RESULTS_DIR, "Resolved_Gene_Trees")
DUP_FILE     = os.path.join(RESULTS_DIR, "Gene_Duplication_Events", "Duplications.tsv")
SP_TREE_FILE = os.path.join(RESULTS_DIR, "Species_Tree",
                             "SpeciesTree_rooted_node_labels.txt")
OUT_FILE     = "/proj/naiss2023-6-65/Sebastian/data/orthofinder/transcript_birth_nodes_2.tsv"

SEP = "_filtered_proteinfasta_TE_filtered_"


# -- Helpers -------------------------------------------------------------------
def get_species(gene_id):
    """
    Extracts the full species prefix from an OrthoFinder gene ID.
    C_maculatus_filtered_proteinfasta_TE_filtered_C_maculatus_g6090.t1_1
    -> C_maculatus_filtered_proteinfasta_TE_filtered
    The species prefix matches the tip label in the species tree.
    """
    idx = gene_id.find(SEP)
    if idx < 0:
        return "unknown"
    return gene_id[: idx + len(SEP) - 1]


def clean_transcript_id(gene_id):
    """
    C_maculatus_filtered_proteinfasta_TE_filtered_C_maculatus_g6090.t1_1
    -> g6090.t1
    Same as the R annotation script.
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


# -- Load species tree ---------------------------------------------------------
# The species tree is needed for MRCA nodes for unduplicated genes.
# For duplicated genes we get what we need from Duplications.tsv.
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

    Returns None if fewer than 2 valid species are found. Single-species
    OGs are handled later as species_specific.
    """
    # Keep only species present in the species tree
    valid = [s for s in species_set if s in sp_tree_tips]
    if len(valid) < 2:
        return None

    # Get the Bio.Phylo clade object for each species
    target_clades = []
    for sp in valid:
        clade = next(sp_tree.find_clades(sp), None)
        if clade is not None:
            target_clades.append(clade)

    if len(target_clades) < 2:
        return None

    # Bio.Phylo common_ancestor finds the MRCA of all provided clades
    mrca  = sp_tree.common_ancestor(target_clades)
    label = mrca.name
    if label and label.strip():
        return label.strip()
    return None


# -- Load Duplications.tsv ----------------------------------------------------
# Duplications.tsv contains all duplication events. Each row is one event:
# it records the orthogroup, the internal node in the gene tree where the
# duplication occurred (Gene Tree Node), and the node in the species tree
# where OrthoFinder places that duplication in evolutionary time
# (Species Tree Node).

# Build a two-level dictionary:
# og_dup_lookup[orthogroup][gene_tree_node] = species_tree_node
# This lets us correlate each gene tree duplication node with its
# corresponding species tree placement.

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


# -- Process gene trees --------------------------------------------------------
# For each orthogroup:
# 1. Parse the resolved gene tree (Newick file, one per OG)
# 2. For each gene (leaf) in the tree, find its birth node
# 3. Record results

# The resolved gene trees are rooted against the species tree and have had
# DLC (Duplication-Loss-Coalescence) reconciliation done. The duplication
# node labels match Duplications.tsv.

print("Processing resolved gene trees ...", flush=True)

tree_files = sorted(f for f in os.listdir(TREE_DIR) if f.endswith('_tree.txt'))
n_files    = len(tree_files)
out_rows   = []

# Counters for final summary
n_duplication      = 0
n_mrca_inferred    = 0
n_species_specific = 0

for i, fname in enumerate(tree_files, 1):
    if i % 2000 == 0:
        print(f"  {i:,}/{n_files:,} trees processed ...", flush=True)

    # Orthogroup ID is the filename minus _tree.txt
    og    = fname[:-9]
    fpath = os.path.join(TREE_DIR, fname)

    with open(fpath) as fh:
        tree_str = fh.read().strip()
    # Some OG files are empty and will be skipped
    if not tree_str:
        continue

    # Use Bio.Phylo to parse the Newick string into a tree object.
    # Internal nodes are labeled n1, n2 etc. These labels become the
    # .name attribute of each internal clade object.
    # Leaves are labeled with the full OrthoFinder gene ID.
    tree = Phylo.read(io.StringIO(tree_str), 'newick')

    # Look up all known duplication nodes for this OG.
    # dup_nodes is a dict: {gene_tree_node_label: species_tree_node_label}
    # e.g. {"n7": "N5", "n12": "C_maculatus_filtered_proteinfasta_TE_filtered"}
    # If the OG has no duplications, this is empty.
    dup_nodes = og_dup_lookup.get(og, {})

    # Get all leaf nodes (the genes themselves)
    terminals = list(tree.get_terminals())

    # Count species representation within this OG.
    # sp_count is used for:
    #   n_copies_in_og: paralog count per species
    #   n_species_in_og + get_mrca_label: to decide whether an unduplicated
    #   gene gets mrca_inferred age or is flagged species_specific
    sp_count = defaultdict(int)
    for leaf in terminals:
        sp_count[get_species(leaf.name)] += 1

    n_species_in_og = len(sp_count)

    # Precompute MRCA for unduplicated genes in this OG (once per OG).
    # og_mrca is the species tree node label for the MRCA of all species
    # in this OG. For example, if the OG contains C_maculatus and
    # T_castaneum, og_mrca = "N4". If only C_maculatus, og_mrca = None.
    og_mrca = get_mrca_label(sp_count.keys())

    # Assign birth nodes for each gene
    for leaf in terminals:
        gene_id_raw   = leaf.name
        species       = get_species(gene_id_raw)
        transcript_id = clean_transcript_id(gene_id_raw)
        gene_id       = get_gene_id(transcript_id)

        # Build the ancestor list from immediate parent toward root.
        #
        # Bio.Phylo's get_path() returns [root ... leaf], including the
        # leaf but EXCLUDING the root. We drop the leaf with path[:-1]
        # and reverse so we walk from immediate parent toward the root.
        #
        # Stopping at the FIRST duplication node gives the birth event
        # that directly created this gene copy, not an older ancestor event.
        #
        # The root is NOT included in this loop. It is checked separately
        # below ONLY if the loop finds nothing. This prevents a root
        # duplication from overwriting a valid birth node found deeper
        # in the tree.

        path       = tree.get_path(leaf)
        ancestors  = list(path[:-1])
        birth_gtn  = None   # gene tree node label (e.g. n7)
        birth_stn  = None   # species tree node label (e.g. N5)
        birth_type = None

        for clade in reversed(ancestors):
            if clade.name and clade.name in dup_nodes:
                birth_gtn  = clade.name
                birth_stn  = dup_nodes[birth_gtn]
                birth_type = "duplication"
                break

        # get_path() excludes the root, so a gene sitting directly under
        # a root duplication node would miss it in the loop above.
        # Check the root explicitly, but only if the loop found nothing.
        if birth_type is None:
            root = tree.root
            if root.name and root.name in dup_nodes:
                birth_gtn  = root.name
                birth_stn  = dup_nodes[birth_gtn]
                birth_type = "duplication"

        # For genes with no duplication ancestor at all
        if birth_type is None:
            if og_mrca is not None:
                # Conserved single-copy gene: assign MRCA of species in
                # this OG as the upper bound age estimate.
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
            'gene_tree_birth_node':    birth_gtn or 'NA',
            'species_tree_birth_node': birth_stn or 'NA',
            'birth_type':              birth_type,
            'n_copies_in_og':          sp_count[species],
            'n_species_in_og':         n_species_in_og,
        })

print(f"  Total gene entries:          {len(out_rows):,}")
print(f"  birth_type=duplication:      {n_duplication:,}")
print(f"  birth_type=mrca_inferred:    {n_mrca_inferred:,}")
print(f"  birth_type=species_specific: {n_species_specific:,}")


# -- Write output --------------------------------------------------------------
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
