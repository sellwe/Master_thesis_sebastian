'''
Merges: 
- Structural annotation (C_maculatus_annotation_unfiltered_fixed.gff3)
- Functional annotation (Cmaculatus_functional_annotation.emapper.annotations),
- Raw OrthoFinder results from 
/proj/coleoptera-genomics-2025/snic2021-6-30/Milena/gene_family_analysis/orthofinder_orthoDB_TE_filtered/Results_May27
(SpeciesTree_rooted_node_labels_May27.txt, 
N0_May27.tsv, 
Duplications_May27.tsv)
- Age rank inferece (transcript_birth_nodes_2.tsv) created in parse_gene_trees_birth_nodes_2.py

Also generates some sanity checks and stats showcased in the GitHub README

'''


library(ape)
library(stringr)
library(tidyverse)
library(dplyr)


# =============================================
# GFF
# =============================================

gff <- read_tsv("C_maculatus_annotation_unfiltered_fixed.gff3",
                comment = "#",
                col_names = FALSE,
                show_col_types = FALSE,
                col_types = cols(.default = col_character()))

colnames(gff) <- c("seqname", "source", "type", "start", "end",
                   "score", "strand", "phase", "attributes")

gene_annotation <- gff %>%
  dplyr::filter(type == "transcript") %>%
  mutate(
    gene_id       = str_extract(attributes, 'gene_id=([^;]+)') %>% str_remove('gene_id='),
    transcript_id = str_extract(attributes, 'ID=([^;]+)')      %>% str_remove('ID=')
  ) %>%
  dplyr::select(gene_id, transcript_id, seqname, start, end, strand)

# sanity check for transcript id extraction from the GFF attributes
gff %>%
  filter(type == "transcript") %>%
  summarise(
    rows = n(),
    unique_ids = n_distinct(
      str_extract(attributes, 'ID=([^;]+)') %>%
        str_remove('ID=')
    )
  )

gff %>%
  filter(type == "transcript") %>%
  mutate(
    transcript_id = str_extract(attributes, 'ID=([^;]+)') %>%
      str_remove('ID=')
  ) %>%
  filter(is.na(transcript_id))

gff %>%
  filter(type == "transcript") %>%
  mutate(
    transcript_id = str_extract(attributes, 'ID=([^;]+)') %>%
      str_remove('ID=')
  ) %>%
  count(transcript_id) %>%
  filter(n > 1) %>%
  print()

# =============================================
# CHROMOSOMAL LOCATION
# =============================================
# Included in the full annotation, but never properly analyzed in this project

x_contigs <- c(
  'utg000057l', 'utg000114l', 'utg000139l', 'utg000191l',
  'utg000326l', 'utg000359l', 'utg000532l', 'utg000602l'
)

y_contigs <- c(
  'utg000322l', 'utg000312c', 'utg000610l', 'utg001235l'
)

confirmed_contigs <- read_tsv(
  "contigs_100kb_with_length.txt",
  col_names = c("contig", "length"),
  show_col_types = FALSE
)

autosomal_contigs <- confirmed_contigs %>%
  dplyr::filter(!contig %in% x_contigs, !contig %in% y_contigs) %>%
  pull(contig)

gene_annotation <- gene_annotation %>%
  mutate(
    chr_loc = case_when(
      seqname %in% x_contigs         ~ "X",
      seqname %in% y_contigs         ~ "Y",
      seqname %in% autosomal_contigs ~ "A",
      TRUE                           ~ "U"
    )
  )

cat("\nChromosomal location distribution:\n")
print(table(gene_annotation$chr_loc))


# =============================================
# EGGNOG FUNCTIONAL ANNOTATION
# =============================================

eggnog_annotation <- read_tsv("Cmaculatus_functional_annotation.emapper.annotations",
                              comment = "#",
                              col_names = FALSE,
                              show_col_types = FALSE)

colnames(eggnog_annotation) <- c(
  "query", "seed_ortholog", "evalue", "score",
  "eggNOG_OGs", "max_annot_lvl", "COG_category",
  "Description", "Preferred_name", "GOs", "EC",
  "KEGG_ko", "KEGG_Pathway", "KEGG_Module",
  "KEGG_Reaction", "KEGG_rclass", "BRITE",
  "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"
)


# =============================================
# ORTHOFINDER DATA (OGs and HOGs)
# =============================================

ortho <- read_tsv("N0_May27.tsv", show_col_types = FALSE)

ortho_c_mac <- ortho %>%
  dplyr::select(HOG, OG, `Gene Tree Parent Clade`,
                `C_maculatus_filtered_proteinfasta_TE_filtered`) %>%
  separate_rows(`C_maculatus_filtered_proteinfasta_TE_filtered`, sep = ",\\s*") %>%
  mutate(
    transcript_id = C_maculatus_filtered_proteinfasta_TE_filtered %>%
      str_remove("^C_maculatus_") %>%
      str_remove("_\\d+$")
  )


# =============================================
# MERGE ANNOTATIONS
# =============================================

full_annotation <- gene_annotation %>%
  left_join(eggnog_annotation, by = c("transcript_id" = "query")) %>%
  left_join(ortho_c_mac, by = "transcript_id") %>%
  dplyr::select(
    gene_id, transcript_id, seqname, chr_loc, start, end, strand,
    Description, Preferred_name, HOG, OG, `Gene Tree Parent Clade`,
    PFAMs, GOs, EC, KEGG_ko, KEGG_Pathway, COG_category, eggNOG_OGs
  )


# =============================================
# SPECIES TREE AGE TABLE (C_mac lineage only)
# =============================================

tree <- read.tree("SpeciesTree_rooted_node_labels_May27.txt")

c_mac_label <- "C_maculatus_filtered_proteinfasta_TE_filtered"

# Walk from C_mac tip to root, collecting all ancestor node indices
get_ancestors <- function(tree, tip_index) {
  ancestors    <- tip_index
  current_node <- tip_index
  repeat {
    parent_edge <- which(tree$edge[, 2] == current_node)
    if (length(parent_edge) == 0) break
    parent_node  <- tree$edge[parent_edge, 1]
    ancestors    <- c(ancestors, parent_node)
    current_node <- parent_node
  }
  return(ancestors)
}

c_mac_tip   <- which(tree$tip.label == c_mac_label)
ancestors   <- get_ancestors(tree, c_mac_tip)
node_depths <- node.depth.edgelength(tree)

node_labels <- sapply(ancestors, function(node) {
  if (node <= length(tree$tip.label)) {
    tree$tip.label[node]
  } else {
    tree$node.label[node - length(tree$tip.label)]
  }
})

branch_lengths <- sapply(ancestors, function(node) {
  edge_index <- which(tree$edge[, 2] == node)
  if (length(edge_index) == 0) return(0)
  tree$edge.length[edge_index]
})

# age_rank 1 = root (oldest), increasing toward C_mac tip (most recent)
node_age_table <- data.frame(
  SpeciesTreeNode      = node_labels,
  node_depth_from_root = node_depths[ancestors],
  branch_length        = branch_lengths,
  stringsAsFactors     = FALSE
) %>%
  arrange(node_depth_from_root) %>%
  mutate(age_rank = row_number())

print("Node age table for C_mac lineage:")
print(node_age_table)


# =============================================
# ALL DUPLICATION EVENTS FOR C_MAC
# Used for the duplication events distribution
# plot in VSCode. One transcript appears many
# times here as each ancestral duplication event
# is listed separately.
# =============================================

dup <- read_tsv("Duplications_May27.tsv", show_col_types = FALSE)

dup_long <- dup %>%
  mutate(all_genes = paste(`Genes 1`, `Genes 2`, sep = ",")) %>%
  separate_rows(all_genes, sep = ",") %>%
  mutate(all_genes = str_trim(all_genes)) %>%
  dplyr::filter(str_detect(all_genes, c_mac_label)) %>%
  dplyr::filter(`Species Tree Node` %in% node_age_table$SpeciesTreeNode) %>%
  left_join(node_age_table, by = c("Species Tree Node" = "SpeciesTreeNode")) %>%
  mutate(
    transcript_id = all_genes %>%
      str_remove("^C_maculatus_filtered_proteinfasta_TE_filtered_") %>%
      str_remove("^C_maculatus_") %>%
      str_remove("_\\d+$")
  )

cat("\nTotal duplication events for each age rank:\n")
print(table(dup_long$age_rank))

# print a .csv if needed

#write.csv(
#  dup_long,
#  "C_mac_dup_long_all_events_May27.csv",
#  row.names = FALSE
#)


# =============================================
# GENE TREE BIRTH NODES
#
# transcript_birth_nodes.tsv was produced by
# parse_gene_trees_birth_nodes.py on UPPMAX.
# It assigns every gene a species_tree_birth_node
# by traversing the resolved gene trees from
# OrthoFinder. Three birth types:
#
#   duplication:      birth node from gene tree traversal.
#                     First duplication node from leaf to root.
#                     Most reliable.
#   mrca_inferred:    gene never duplicated; age is MRCA of all
#                     species in OG. Upper bound, but informative.
#   species_specific: no orthologs in any other species and no
#                     duplication recorded at root. Age unknown (NA).
#                     For C_maculatus in this dataset: 0 genes.
#                     Kept in code for generalizability to other species.
#
# OrthoFinder runs on one protein per gene, keeping the longest isoform.
# This is most often the .t1 transcript but can be .t2 or higher.
# The representative isoform receives age and HOG data. All other isoforms
# for that gene retain NA for all age columns. The full annotation retains
# all isoforms from the structural annotation for completeness.
# =============================================

birth_nodes_raw <- read_tsv(
  "transcript_birth_nodes_2.tsv",
  show_col_types = FALSE
) %>%
  dplyr::filter(species == c_mac_label)

cat("\nBirth type distribution from gene trees:\n")
print(table(birth_nodes_raw$birth_type))

# Join species_tree_birth_node with node_age_table to get age_rank.
# species_specific genes have species_tree_birth_node = "NA" and will
# get no match here, leaving age_rank as NA.
birth_nodes <- birth_nodes_raw %>%
  left_join(
    node_age_table %>%
      dplyr::select(SpeciesTreeNode, age_rank,
                    node_depth_from_root, branch_length),
    by = c("species_tree_birth_node" = "SpeciesTreeNode")
  ) %>%
  dplyr::select(
    transcript_id,
    gene_id,
    orthogroup,
    species_tree_birth_node,
    gene_tree_birth_node,
    birth_type,
    age_rank,
    node_depth_from_root,
    branch_length,
    n_copies_in_og,
    n_species_in_og
  )

cat("\nAge rank distribution (representative isoforms only):\n")
print(table(birth_nodes$age_rank, useNA = "ifany"))


# =============================================
# ADD AGE TO FULL ANNOTATION
#
# Join by transcript_id only. OrthoFinder used one protein per gene
# (longest isoform, usually .t1 but not always). Only the representative
# isoform receives age and HOG data. All other isoforms retain NA for
# all age columns, and are excluded from paralog analyses by the
# existing filter(!is.na(HOG)) and filter(!is.na(age_rank)) steps.
# =============================================

age_cols <- c("species_tree_birth_node", "gene_tree_birth_node",
              "birth_type", "age_rank",
              "node_depth_from_root", "branch_length",
              "n_copies_in_og", "n_species_in_og")

full_annotation_with_age <- full_annotation %>%
  left_join(
    birth_nodes %>% dplyr::select(transcript_id, all_of(age_cols)),
    by = "transcript_id"
  )

cat(paste0("\nTranscripts matched by transcript_id: ",
           sum(!is.na(full_annotation_with_age$age_rank)), "\n"))
cat(paste0("Transcripts without age (non-representative isoforms or absent from OrthoFinder): ",
           sum(is.na(full_annotation_with_age$age_rank)), "\n"))

cat("\nFinal age rank distribution in full annotation:\n")
print(table(full_annotation_with_age$age_rank, useNA = "ifany"))

cat("\nBirth type distribution in full annotation:\n")
print(table(full_annotation_with_age$birth_type, useNA = "ifany"))


# =================================================
# ISOFORM SANITY CHECK
# Are there genes with more than one isoform that
# both received a HOG assignment from OrthoFinder?
# If yes, that means two proteins from the same gene
# entered OrthoFinder, which violates the one-protein-
# per-gene design and creates non-independence.
# =================================================

isoform_check <- full_annotation_with_age %>%
  filter(!is.na(HOG)) %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms    = n(),
    transcripts   = paste(sort(transcript_id), collapse = ", "),
    .groups       = "drop"
  ) %>%
  filter(n_isoforms > 1)

cat("\nGenes with multiple HOG-assigned isoforms (bad case):", nrow(isoform_check), "\n")

if (nrow(isoform_check) > 0) {
  cat("These genes have more than one isoform in OrthoFinder:\n")
  print(isoform_check)
} else {
  cat("All genes have exactly one HOG-assigned isoform. Clean.\n")
}

# =================================================
# CHECK 2: Any gene with multiple isoforms in the
# full structural annotation (regardless of HOG)?
# This would mean the unfiltered GTF contains both
# .t1 and .t2 for the same gene, which is expected
# since we deliberately used the non-isoform-filtered
# annotation. 
# =================================================

all_isoforms <- full_annotation_with_age %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms  = n(),
    transcripts = paste(sort(transcript_id), collapse = ", "),
    .groups     = "drop"
  ) %>%
  filter(n_isoforms > 1)

cat("\nGenes with multiple isoforms in full structural annotation:", nrow(all_isoforms), "\n")
cat("Total extra isoforms (beyond one per gene):", sum(all_isoforms$n_isoforms - 1), "\n")

cat("\nIsoform count distribution:\n")
isoform_dist <- all_isoforms %>%
  group_by(n_isoforms) %>%
  summarise(n_genes = n(), .groups = "drop") %>%
  arrange(n_isoforms)

print(isoform_dist)

# =================================================
# CHECK 3: Any gene with multiple isoforms that both
# received an age rank?
# Same logic as the HOG check but for age_rank.
# Should also be 0 since age rank comes from the
# same OrthoFinder one-protein-per-gene run.
# =================================================

age_isoform_check <- full_annotation_with_age %>%
  filter(!is.na(age_rank)) %>%
  group_by(gene_id) %>%
  summarise(
    n_isoforms  = n(),
    transcripts = paste(sort(transcript_id), collapse = ", "),
    .groups     = "drop"
  ) %>%
  filter(n_isoforms > 1)

cat("\nGenes with multiple age-rank-assigned isoforms (bad case):", nrow(age_isoform_check), "\n")

if (nrow(age_isoform_check) > 0) {
  cat("These genes have more than one isoform with an age rank:\n")
  print(age_isoform_check)
} else {
  cat("All genes have exactly one age-rank-assigned isoform. Clean.\n")
}


# =============================================
# SAVE
# =============================================

write.csv(
  full_annotation_with_age,
  "C_mac_full_annotation_with_age_fixed_2.csv",
  row.names = FALSE
)

cat("\nDone. Saved to C_mac_full_annotation_with_age_fixed_2.csv\n")

