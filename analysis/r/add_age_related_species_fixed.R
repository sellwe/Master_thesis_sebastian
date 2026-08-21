
library(dplyr)
library(tidyverse)
library(ape)
library(tidyr)
library(stringr)


# =============================================
# LOAD GFFs
# =============================================

c_mac_gff <- read_tsv("C_maculatus_annotation_unfiltered_fixed.gff3",
                      comment = "#", col_names = FALSE, show_col_types = FALSE,
                      col_types = cols(.default = col_character()))
colnames(c_mac_gff) <- c("seqname", "source", "type", "start", "end",
                         "score", "strand", "phase", "attributes")

c_mac_annotation <- c_mac_gff %>%
  dplyr::filter(type == "transcript") %>%
  mutate(
    gene_id       = str_extract(attributes, 'gene_id=([^;]+)') %>% str_remove('gene_id='),
    transcript_id = str_extract(attributes, 'ID=([^;]+)')      %>% str_remove('ID=')
  ) %>%
  dplyr::select(gene_id, transcript_id, seqname, start, end, strand)

# All related species use isoform-filtered GFFs with BRAKER attribute formatting.
# transcript_id is parsed from ID=, gene_id from Parent= (differs from C_mac GFF).
read_braker_gff <- function(path) {
  gff <- read_tsv(path, comment = "#", col_names = FALSE,
                  show_col_types = FALSE, col_types = cols(.default = col_character()))
  colnames(gff) <- c("seqname", "source", "type", "start", "end",
                     "score", "strand", "phase", "attributes")
  gff %>%
    dplyr::filter(type == "transcript") %>%
    mutate(
      transcript_id = str_extract(attributes, 'ID=([^;]+)')     %>% str_remove('ID='),
      gene_id       = str_extract(attributes, 'Parent=([^;]+)') %>% str_remove('Parent=')
    ) %>%
    dplyr::select(gene_id, transcript_id, seqname, start, end, strand)
}

a_obtectus_annotation       <- read_braker_gff("A_obtectus_annotation_braker_isoform_filtered.gff")
b_siliquastri_annotation    <- read_braker_gff("B_siliquastri_annotation_braker_isoform_filtered.gff")
c_chinensis_annotation      <- read_braker_gff("C_chinensis_annotation_braker_isoform_filtered.gff")
c_septempunctata_annotation <- read_braker_gff("C_septempunctata_annotation_braker_isoform_filtered.gff")
t_castaneum_annotation      <- read_braker_gff("T_castaneum_annotation_braker_isoform_filtered.gff")


# =============================================
# SHARED RESOURCES
# =============================================

tree  <- read.tree("SpeciesTree_rooted_node_labels_May27.txt")
ortho <- read_tsv("N0_May27.tsv", show_col_types = FALSE)

# transcript_birth_nodes.tsv was produced by parse_gene_trees_birth_nodes.py
# on UPPMAX. Contains one row per gene per species. All species are included
# so we filter per species inside the function below.
birth_nodes_all <- read_tsv("transcript_birth_nodes_2.tsv", show_col_types = FALSE)

cat("Loaded birth nodes for", n_distinct(birth_nodes_all$species), "species\n")
cat("Total entries:", nrow(birth_nodes_all), "\n")

# Walk from a tip to the root, collecting all ancestor node indices.
# Used to build the lineage-specific age table for each species.
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


# =============================================
# CORE FUNCTION
# =============================================

# Builds a full annotation table for one species, including HOG membership
# and gene age ranks from gene tree birth nodes. Each species gets its own
# lineage-specific age table since the number of ancestors and node depths
# differ across the phylogeny.
build_species_annotation <- function(gene_annotation, tip_label, output_path) {
  
  species_prefix <- paste0(tip_label, "_")
  short_prefix   <- paste0(str_extract(tip_label, "^[^_]+_[^_]+"), "_")
  
  cat("\n=== Processing:", tip_label, "===\n")
  
  # ------------------------------------------
  # 1. Lineage-specific age table
  #
  # Maps species_tree_birth_node labels to age_rank
  # for this species. Age rank 1 = root (oldest),
  # increasing toward the species tip (most recent).
  # ------------------------------------------
  tip_index   <- which(tree$tip.label == tip_label)
  ancestors   <- get_ancestors(tree, tip_index)
  node_depths <- node.depth.edgelength(tree)
  
  node_labels <- sapply(ancestors, function(node) {
    if (node <= length(tree$tip.label)) tree$tip.label[node]
    else tree$node.label[node - length(tree$tip.label)]
  })
  
  branch_lengths <- sapply(ancestors, function(node) {
    edge_index <- which(tree$edge[, 2] == node)
    if (length(edge_index) == 0) return(0)
    tree$edge.length[edge_index]
  })
  
  node_age_table <- data.frame(
    SpeciesTreeNode      = node_labels,
    node_depth_from_root = node_depths[ancestors],
    branch_length        = branch_lengths,
    stringsAsFactors     = FALSE
  ) %>%
    arrange(node_depth_from_root) %>%
    mutate(age_rank = row_number())
  
  cat("Lineage nodes found:", nrow(node_age_table), "\n")
  print(node_age_table)
  
  # ------------------------------------------
  # 2. HOG memberships from N0 table
  # ------------------------------------------
  ortho_species <- ortho %>%
    dplyr::select(HOG, OG, `Gene Tree Parent Clade`, all_of(tip_label)) %>%
    separate_rows(all_of(tip_label), sep = ",\\s*") %>%
    mutate(
      transcript_id = .data[[tip_label]] %>%
        str_remove(fixed(species_prefix)) %>%
        str_remove(fixed(short_prefix)) %>%
        str_remove("_\\d+$")
    ) %>%
    dplyr::select(-all_of(tip_label))
  
  # ------------------------------------------
  # 3. Gene tree birth ranks
  #
  # Filter transcript_birth_nodes.tsv to this species,
  # then join species_tree_birth_node with the lineage-
  # specific node_age_table to get age_rank.
  #
  # Join by transcript_id only. OrthoFinder used one
  # isoform per gene, so only representative isoforms
  # receive age data. Non-representative isoforms
  # retain NA for all age columns.
  # ------------------------------------------
  birth_nodes_sp <- birth_nodes_all %>%
    dplyr::filter(species == tip_label) %>%
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
  
  cat("\nBirth type distribution:\n")
  print(table(birth_nodes_sp$birth_type))
  
  cat("\nAge rank distribution (representative isoforms only):\n")
  print(table(birth_nodes_sp$age_rank, useNA = "ifany"))
  
  age_cols <- c("species_tree_birth_node", "gene_tree_birth_node",
                "birth_type", "age_rank",
                "node_depth_from_root", "branch_length",
                "n_copies_in_og", "n_species_in_og")
  
  # ------------------------------------------
  # 4. Build full annotation
  # ------------------------------------------
  
  # Structural annotation + HOG membership
  full_annotation <- gene_annotation %>%
    left_join(
      ortho_species %>%
        dplyr::select(transcript_id, HOG, OG, `Gene Tree Parent Clade`),
      by = "transcript_id"
    )
  
  # Join age by transcript_id only
  full_annotation <- full_annotation %>%
    left_join(
      birth_nodes_sp %>% dplyr::select(transcript_id, all_of(age_cols)),
      by = "transcript_id"
    )
  
  cat(paste0("\nTranscripts matched by transcript_id: ",
             sum(!is.na(full_annotation$age_rank)), "\n"))
  cat(paste0("Transcripts without age (non-representative isoforms or absent from OrthoFinder): ",
             sum(is.na(full_annotation$age_rank)), "\n"))
  
  full_annotation <- full_annotation %>%
    dplyr::select(
      gene_id, transcript_id, seqname, start, end, strand,
      HOG, OG, `Gene Tree Parent Clade`,
      age_rank,
      species_tree_birth_node, gene_tree_birth_node,
      birth_type,
      node_depth_from_root, branch_length,
      n_copies_in_og, n_species_in_og
    )
  
  cat("\nFull annotation rows:", nrow(full_annotation), "\n")
  cat("Transcripts with age_rank assigned:",
      sum(!is.na(full_annotation$age_rank)), "\n")
  
  cat("\nFinal age rank distribution:\n")
  print(table(full_annotation$age_rank, useNA = "ifany"))
  
  write.csv(full_annotation, output_path, row.names = FALSE)
  cat("Saved to:", output_path, "\n")
  
  return(invisible(full_annotation))
}


# =============================================
# RUN ALL SPECIES
# =============================================

build_species_annotation(
  gene_annotation = a_obtectus_annotation,
  tip_label       = "A_obtectus_filtered_proteinfasta_TE_filtered",
  output_path     = "A_obtectus_full_annotation_with_age_fixed_2.csv"
)

build_species_annotation(
  gene_annotation = b_siliquastri_annotation,
  tip_label       = "B_siliquastri_filtered_proteinfasta_TE_filtered",
  output_path     = "B_siliquastri_full_annotation_with_age_fixed_2.csv"
)

build_species_annotation(
  gene_annotation = c_chinensis_annotation,
  tip_label       = "C_chinensis_filtered_proteinfasta_TE_filtered",
  output_path     = "C_chinensis_full_annotation_with_age_fixed_2.csv"
)

build_species_annotation(
  gene_annotation = c_septempunctata_annotation,
  tip_label       = "C_septempunctata_filtered_proteinfasta_TE_filtered",
  output_path     = "C_septempunctata_full_annotation_with_age_fixed_2.csv"
)

build_species_annotation(
  gene_annotation = t_castaneum_annotation,
  tip_label       = "T_castaneum_filtered_proteinfasta_TE_filtered",
  output_path     = "T_castaneum_full_annotation_with_age_fixed_2.csv"
)

