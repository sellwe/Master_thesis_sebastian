# Master_thesis_sebastian

## Repository Structure

Master_thesis_sebastian/  
├── scripts/ # Shell scripts from UPPMAX   
│ ├── c_mac_reference_genome/ # Reference genome, Genome annotation and preparation  
│ └── dataset_1_dominance_kaufmann/ # RNA-seq analysis pipeline  
├── analysis/ # Statistical analysis and visualization  
│ ├── r/ # R scripts   
│ └── python/ # Python scripts for plotting and additional analysis  
├── metadata/ # Sample metadata and information  
├── .gitignore # Files to exclude from version control  
└── README.md # This file  

Here I am investigating how paralog age might predict sex-biased gene expression in C. maculatus, with a subtheme of how gene duplications might resolve intralocus sexual conflicts. 

Some research questions to be answered: 
 
* Is sex-biased expression consistent among paralogs within a gene family?

* Does paralog age predict whether a copy evolves sex-biased expression?

* Are older copies within a family more or less biased than younger ones? 

* Does this depend on when the gene family itself originated? 

# Datasets

## Reference Genome

The reference genome comes from the paper "Y-Linked Copy Number Polymorphism of Target of Rapamycin Is Associated with Sexual Size Dimorphism in Seed Beetles" by Kaufmann et. al 2023. Male virgin C_maculatus abdominal tissue samples from the Lomé population. I am using the small male Y haplotype assembly from the paper, as it is more continuous, and the small Y haplotype is the most abundant haplotype in the population. 

**Summary:**
| Metric        | Value   |
|--------------|--------|
| Genome size  | 1.2 Gbp |
| Genes        | 35,865  |
| Repeats      | 72%     |
| Y assembled  | 10 Mbp  |

## Annotations  
Scripts titled _unfiltered are based on the non-isoform filtered annotation, to preserve all isoform information for downstream use. 

**Source files:**  
| Softlinked File Name | Source path |
|------|-------------|
| C_maculatus_annotation_nonfiltered.gtf | /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/braker/braker.gtf |
| C_maculatus_annotation_isoform_filtered.gff | /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/braker/braker_isoform_filtered.gff |
| braker_proteins.fa | /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/braker/braker.aa |
| C_maculatus_assembly.fna | /proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/C_maculatus/assembly_genomic.fna.masked |
| OrthoFinder results | /proj/naiss2023-6-65/Milena/gene_family_analysis/orthofinder_orthoDB_TE_filtered/Results_May27/ |

### Structural annotation:  
C_maculatus_annotation_nonfiltered.gtf was converted to a gff3 file using:   
agat_convert_sp_gxf2gxf.pl \  
    --gxf C_maculatus_annotation_nonfiltered.gtf \  
    -o braker_unfiltered.gff3  

Then, since STAR requires a gtf file to be run, this gff3 file was converted back to a .gtf file using:  
gffread braker_unfiltered.gff3 -T -o C_maculatus_annotation_unfiltered_fixed.gtf 

This standardizes BRAKER attribute formatting and ensures consistent transcript_id fields across all features. The resulting C_maculatus_annotation_unfiltered_fixed.gtf was used consistently for the STAR index and Salmon transcriptome.(why some scripts are named _consistent).  

After the mapping softwares were finished a final gff3 file was created for merging with the functional annotation and downstream comparative genomic analyses (gff is easier to handle) (**run_agat_gtf_to_gff3.sh**). 

### Functional annotation:  
eggNOG-mapper was run on the BRAKER protein sequences to assign functional annotation including GO terms, KEGG pathways, COG categories and PFAM domains (**run_eggnog.sh**). 

### Orthology Inference:    
The structural annotation, eggNOG results and OrthoFinder output (N0.tsv) were merged in R to produce the full annotation table (**create_full_annotation_fixed.R**). Each transcript was assigned a chromosomal location based on known X, Y and autosomal contigs.  
Sex chromosome contigs had already been identified to:  
x_contigs <- c(
  'utg000057l', 'utg000114l', 'utg000139l', 'utg000191l',
  'utg000326l', 'utg000359l', 'utg000532l', 'utg000602l'
)

y_contigs <- c(
  'utg000322l', 'utg000312c', 'utg000610l', 'utg001235l'
)  

Contigs shorter than 100 kb were excluded, consistent with the paper. Unassigned contigs are labeled U. 

| Location   | Transcripts |
|------------|------------|
| Autosomal  | 28,651     |
| X-linked   | 1,827      |
| Y-linked   | 334        |
| Unassigned | 7,176      |

Gene age ranks were assigned by traversing the OrthoFinder resolved gene trees on UPPMAX (**parse_gene_trees_birth_nodes.py**, run via **run_parse_gene_trees.sh**). The script runs on the full OrthoFinder dataset and assigns a birth node to every gene across all species, which can be used for future studies of this phylogeny. Each gene gets one of three birth types:  
| Birth type        | Description                                                                 |
|------------------|-----------------------------------------------------------------------------|
| duplication      | Birth node from gene tree traversal. Most reliable.                         |
| mrca_inferred    | Gene never duplicated; age is the MRCA of all species in the orthogroup. Upper bound. |
| species_specific | No orthologs detected in any other species. Age unknown (NA).              |

Age ranks were mapped to nodes along the C_maculatus lineage in the species tree, from root (rank 1, oldest) to the C_maculatus tip (rank 9, most recent). 
| SpeciesTreeNode        | node_depth_from_root | branch_length | age_rank |
|------------------------|----------------------|---------------|----------|
| N0                     | 0.000                | 0.000         | 1        |
| N1                     | 0.351                | 0.351         | 2        |
| N3                     | 0.398                | 0.047         | 3        |
| N4                     | 0.421                | 0.023         | 4        |
| N6                     | 0.458                | 0.037         | 5        |
| N8                     | 0.629                | 0.171         | 6        |
| N11                    | 0.661                | 0.032         | 7        |
| N12                    | 0.697                | 0.036         | 8        |
| C_maculatus (tip)      | 0.759                | 0.062         | 9        | 

**Phylogeny with age ranks:**
![alt text](image-3.png)

**All duplication events in the C.mac lineage:**  
![alt text](image-1.png)

Birth type distribution across all C_maculatus transcripts in OrthoFinder:  
| Birth type        | Transcripts |
|------------------|------------|
| duplication      | 20,026     |
| mrca_inferred    | 6,733      |
| species_specific | 559        |

Age rank distribution across representative isoforms:
 
| Age rank | 1     | 2     | 3    | 4    | 5    | 6     | 7    | 8    | 9      | NA   |
|----------|-------|-------|------|------|------|-------|------|------|--------|------|
| Transcripts    | 4,447 | 3,382 | 561  | 624  | 280  | 1,511 | 282  | 667  | 15,005 | 559  |

Rank 9 contains the most transcripts as it captures all C_maculatus-specific and recently diverged genes. Makes sense in light of C. macs large and highly repetitive genome. Many of these might be TE-induced duplications.   
Age ranks were propagated to all isoforms via a two-step join: first by transcript_id (26,759 transcripts matched directly), then by gene_id as fallback for non-representative isoforms of genes present in OrthoFinder (2,207 additional transcripts rescued). 9,022 transcripts remain without an age rank; these are genes fully absent from OrthoFinder or species-specific singletons.  

Final age rank and birth type distribution in the full annotation:

| Age rank | 1     | 2     | 3    | 4    | 5    | 6     | 7    | 8    | 9      | NA   |
|----------|-------|-------|------|------|------|-------|------|------|--------|------|
| Transcripts    | 4,995 | 3,695 | 625  | 692  | 301  | 1,609 | 297  | 710  | 16,042 | 9,022 | 

| Birth type        | Transcripts |
|------------------|------------|
| duplication      | 21,572     |
| mrca_inferred    | 7,394      |
| species_specific | 572        |
| NA               | 8,450      |

### Related species
Out of curiosity, the same rank approach was applied to five related species using (**add_age_related_species.R**):  
* Coccinella septempunctata
* Tribolium castaneum
* Acanthoscelides obtectus
* Bruchidius siliquastri
* Callosobruchus chinensis 

All species annotation files (braker.gtf) from their respective folder in:  
/proj/naiss2023-6-65/Milena/annotation_pipeline/only_orthodb_annotation/

Script function: builds a lineage specific age table for each species by walking its path through the shared species tree, then joins gene tree birth nodes from transcript_birth_nodes.tsv to assign age ranks. The same two-step transcript_id / gene_id join handles isoforms. Since transcript_birth_ndoes.tsv covers all species in the phylogeny, no additional UPPMAX jobs were needed, and deeper analyses are available in the fututre. 

![alt text](image-4.png)

The same pattern is seen for C. septempunctata, A. obtectus and C. chinensis. T.castaneum and B. siliquastri does not share the pattern. Hypethesis being that large and repetitive genomes lead to gene fsmaily expansion through TEs. 

## RNA-Seq Data

The dataset comes from the paper "Sex-Specific Dominance of Gene Expression in Seed Beetles" by Kaufmann et.al 2024. It provides sex-specific expression data across multple genotypes (will be accounted for). 

- 10 gen full sib inbred lines, Lomé Population  
- Three pairwise crosses of six isogenic lines.  
- Six homozygous lines (no signs of inbreeding depression), three heterozygote.  
- Each isogenic line used as both maternal and paternal (reciprocal crosses), giving 4 genotypes per cross  
- 3 replicates per genotypes per sexes, 24 samples per cross, 72 samples in total.  
- Each sample consists of RNA extracted from abdominal tissues from a pool of 6 virgin beetle abdomens.  

The data was downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB70958&o=acc_s%3Aa with the script (**download_PRJEB70968.sh**).   

Run FastQC and MultiQC (**run-fastqc_multiqc.sh**)  
Run fastp for trimming (**run_fastp_multiqc.sh**)  
Run fastqc and multiqc again to confirm improvements (**run_fastqc_multiqc_post_trim.sh**).

### Metadata

The metadata was found from the same SRA website with accesion-nr PREJB70958.  
The original metadata from the study is found in **Philipp_dominance_metafile.xlsx** and **Philipp_dominance_notes_meta.pdf**.   

Label-corrected metadata is in **dominance_meta_corrected.xlxs** (all 72 entries), and **dominance_meta_corrected_outlier_corrected.xlsx** and **dominance_meta_corrected_outlier_corrected.csv** (where the outliers are removed).   

**Sex-correction:** Sample sexes were determined by linking original FASTQ filenames submitted to SRA (ex. TF-2581-3_S3_L001_R1_001.fastq.gz) to the TF ID in the original Excel file, which has the correct sexes in the Sex column.  

**Genotype-translation:** Cross IDs from the original study were translated to letter codes similar to those in the published study as follows:  
| Cross ID | Code |
|----------|------|
| 13:20    | A    |
| 42:13    | B    |
| 21:5     | C    |
| 1:11     | D    |
| 47:1     | E    |
| 4:18     | F    | 

| Cross | Samples   | Genotypes          |
|-------|-----------|--------------------|
| A x B | 1-24      | AA, AB, BA, BB     |
| C x D | 25-48     | CC, CD, DC, DD     |
| E x F | 49-72     | EE, EF, FE, FF     |

Reciprocal crosses were collapsed into a single Genotype column (AB + BA = AB), as cross direction did not affect results in the original study. Genotype is treated as a background variable in the differential expression analysis.

**Outlier correction:** After PCA visualization one sample (ERR12383283, DD) was changed from male to female due to clustering. Two samples were removed due to ambiguous sex:   
ERR12383297 (male, FE) and   ERR12383303 (male, EF). 
Three remaining samples are suspected of ambigous sex as they stray from the respective clusters in the PCA, but are kept (ERR12383254 (female, BA), ERR12383278 (male, AA), ERR12383310 (male, FF)). 

Final sample count: 70 samples, of whom 37 are female, and 33 are male. 

**Early PCA plot indicating outliers:**
![alt text](<Screenshot 2025-12-12 140934.png>) 

# Mapping methods

Three methods were tested to assess how each handles multi-mapping, which is a big concern when studying paralogs. All three were run on the transcript level using RNA-seq data from the Kaufmann 2024 dominance dataset.

- Salmon mapping-based mode (quasi-mapping / selective alignment)
- STAR + featureCounts
- STAR + Salmon alignment-based mode

After comparing mapping rates, DE signal and correlation between methods, salmon-map
was selected for the main analysis. 

For details on the two Salmon modes, see:
https://salmon.readthedocs.io/en/latest/salmon.html#

## Salmon-mapping based mode (Quasi-mapping) 

Salmon maps directly to the transcriptome, on the fragment level (one fragment = read pairs = one RNA-molecule). Quasi-mapping uses a k-mer index rather than full base-by-base alignment. Salmon runs a sliding window across each read and identifies candidate transcript matches based on consistent k-mer hits. Reads compatible withthe same set of transcripts are grouped into equivalence classes, and a probabilistic EM algorithm resolves which transcript each fragment most likely originated from. This makes it well suited for paralog studies since ambiguous reads are modelled
rather than discarded.

**Step 1: Build transcriptome**    
A transcriptome was created from `C_maculatus_assembly.fna` and
`C_maculatus_annotation_unfiltered_fixed.gtf` using gffread
(**create_transcript_unfiltered_consistent.sh**).

**Step 2: Build gentrome and decoy file**  
The whole-genome decoy approach was used, where genome sequences serve as decoys. Fragments that map best to a decoy are discarded. A gentrome (transcripts followed by the full genome) and a decoy names file were generated (**generate_gentrome_decoys_consistent.sh**).

**Step 3: Build Salmon index**  
A k-mer index was built from the gentrome and decoy files
(**create_salmon_index_unfiltered_consistent.sh**).

**Step 4: Quantify**
Salmon was run with the following flags (**run_salmon_map_consistent.sh**):

- `--gcBias` corrects for GC-content bias during quantification
- `--seqBias` corrects for sequence-specific bias at fragment starts
- `--validateMappings` enables selective alignment mode (now default)

**Step 5: Differential expression in R**
(**salmon_map_dominance_consistent_script.R**)  

- tximport on the transcript level
- Pre-filter transcript list exported before expression filtering as
  `prefilter_transcripts_annotated_May27.csv`, annotated with HOG and gene_id.
  Used downstream to get gene family sizes at the mapped genome level.
- Expression filter: ≥5 counts in at least 5 samples. 
- DESeq2: `~ Genotype + Sex`, contrast Male vs. Female. Genotypes from the cross-design in the paper taken into account here. 
- VST normalization for PCA
- Results annotated with `C_mac_full_annotation_with_age_fixed.csv`, which includes
structural annotation, eggNOG functional annotation, OrthoFinder HOG membership
and gene age ranks.

**Step 6: Visualization in Python**
Results were combined with the full annotation and imported into VS Code for PCA and volcano plots (**salmon_map_unfiltered_plotting_transcript_new_filtering.ipynb**).

## STAR (with featureCounts)

STAR aligns RNA reads to the genome and is splice-junction aware. The alignment was split across three scripts due to UPPMAX time limits, with the continuation scripts checking for already completed samples before restarting. 
(**star_alignment_dominance.sh**, **star_alignment_dominance_continuation.sh**, **star_alignment_dominance_continuation_2.sh**).

**Step 1: Build STAR index and align**

Index was built with:
- `--sjdbGTFfile` and `--sjdbOverhang 149` (read length 150 - 1) for splice-junction awareness

Alignment flags:
- `--outSAMtype BAM SortedByCoordinate` sorts output by genomic location
- `--quantMode GeneCounts` performs gene-level quantification
- `--twopassMode Basic` enables novel junction discovery
- `--outFilterMultimapNmax 20` discards reads mapping to more than 20 locations

**Step 2: Mark duplicates and index BAM files**

Picard was used to flag duplicate reads based on identical start positions. samtools indexed the resulting BAM files. Run as a SLURM array job, one task per sample (**run_picard_samtools.sh**).

featureCounts was run in four modes on the exon level using the Picard-marked BAM files. Only mode 4 was used for the downstream comparison with Salmon, as it is the closest equivalent.

| Mode | Level | Multimapping | Script |
|------|-------|-------------|--------|
| 1 | Gene | No | run_subread_featurecounts.sh |
| 2 | Gene | Yes, fractional | run_subread_featurecounts.sh |
| 3 | Transcript | No | run_subread_featurecounts.sh |
| 4 | Transcript | Yes, fractional | run_subread_featureCounts_transcript_multi.sh

Multimapped reads are split fractionally across all targets. This does not use sequence uniqueness or transcript abundance models the way Salmon does, which risks inflating counts or missassigning reads for paralogous transcripts. 

**Step 4: Differential expression in R** (**star_DE_analysis.R**)

- Exon-level counts aggregated to transcript level by summing
- Only mode 4 (transcript-level multimappers) used for DE analysis
- Fractional counts rounded to integers before DESeq2 input
- Expression filter: ≥5 counts in at least 5 samples
- DESeq2: `~ Genotype + Sex`, contrast Male vs. Female
- VST normalization for PCA
- Results annotated with `C_mac_full_annotation_with_age_fixed.csv`, which includes
structural annotation, eggNOG functional annotation, OrthoFinder HOG membership
and gene age ranks.

**Step 5: Visualization in Python**
(**STAR_plotting_new_filtering.ipynb**)

## Salmon-alingment based mode 

In this mode, STAR first aligns reads to the genome and outputs transcript-coordinate BAM files. Salmon then quantifies from those alignments using the same transcriptome as salmon-map. This is the most conservative of the three methods: reads that are too ambiguous for STAR alignment are discarded before Salmon ever sees them.

**Step 1: STAR alignment to transcriptome coordinates**

The same STAR index from the featureCounts pipeline was reused. STAR was run with `--quantMode TranscriptomeSAM` to produce BAM files in transcript coordinates for Salmon. Run was split across three scripts due to UPPMAX time limits, processing
24 samples each (**star_transcriptome_for_salmon_1.sh**, **star_transcriptome_for_salmon_2.sh**, **star_transcriptome_for_salmon_3.sh**).

Key flags differing from the featureCounts STAR run:
- `--quantMode TranscriptomeSAM` outputs BAM aligned to transcript coordinates
- `--outSAMattributes NH HI AS nM XS GX GN` adds tags required by Salmon
- `--alignEndsType EndToEnd` disables soft-clipping
- `--outSAMmapqUnique 60` sets mapping quality for uniquely mapped reads
- `--winAnchorMultimapNmax 100` allows more anchor points per window
- `--outFilterMultimapNmax 20` discards reads mapping to more than 20 locations

These flags differ from the featureCounts STAR run because the output goes to Salmon rather than featureCounts. The additional SAM attributes are required by Salmon to process the BAM. `--winAnchorMultimapNmax 100` increases the number
of candidate anchor positions, giving STAR more opportunity to find valid alignments for reads near repetitive regions before passing them to Salmon. `--alignEndsType EndToEnd` disables soft-clipping to ensure only full-match reads enter Salmon's probabilistic model, though this also makes the alignment step slightly more conservative than the featureCounts run.

**Step 2: Salmon quantification from BAM files**

Salmon used the transcript-coordinate BAM files from STAR and the same transcriptome as salmon-map (**run_salmon_align_star_consistent.sh**).

Flags used:
- `--gcBias` corrects for GC-content bias
- `--seqBias` corrects for sequence-specific bias at fragment starts

**Step 3: Differential expression in R**
(**salmon_align_dominance_consistent_script.R**)

- tximport on the transcript level
- Expression filter: ≥5 counts in at least 5 samples
- DESeq2: `~ Genotype + Sex`, contrast Male vs. Female
- VST normalization for PCA
- Results annotated with `C_mac_full_annotation_with_age_fixed.csv`

**Step 4: Visualization in Python** (**salmon_align_plotting_new_filtering.ipynb**)

# Mapping software comparison  

The three methods were compared on DE signal, expression agreement and mapping statistics. All results are filtered on (≥5 counts in ≥5 samples) and use the DESeq2 design `~ Genotype + Sex`. Log file summaries were computed by parsing the output logs from each software across all 70 samples (**mapping_software_comparison.R**).

## Differential Expression Analysis (male vs. female)

## PCA Plots 
### Salmon-Map  
![alt text](image-10.png) 
### Salmon-Align  
![alt text](image.png)
### STAR  
![alt text](image-8.png)

## Volcano Plots  
### Salmon-Map  
![alt text](image-11.png)
### Salmon-Align  
![alt text](image-7.png) 
### STAR  
![alt text](image-9.png)



## Post-DE method comparison summary table  

| Method              | Total Transcripts | Transcripts Retained | Significant DE (padj < 0.05) | Significant DE (padj < 0.05 & abs(log2FC) > 1) | Higher in Males | Higher in Females | PCA1 Variance Explained (%) |
|---------------------|------------------|----------------------|------------------------------|------------------------------------------------|-----------------|-------------------|-----------------------------|
| Salmon_mapping      | 36,382           | 17,574               | 13,923                       | 7,270                                          | 5,070           | 2,200             | 65.4                        |
| STAR_featureCounts  | 37,989           | 17,566               | 14,101                       | 6,964                                          | 4,843           | 2,121             | 69.1                        |
| Salmon_alignment    | 37,989           | 16,563               | 12,645                       | 6,382                                          | 4,277           | 2,105             | 66.0                        |


## Expression correlation   
Pairwise comparisons of log10(baseMean + 1) per transcript between methods, computed on transcripts present in all three DE results (**mapping_software_comparison_correlation.ipynb**). 

### STAR vs. Salmon map
![alt text](image-6.png)

### STAR vs. Salmon align
![alt text](image-5.png)  

### Salmon map vs Salmon align

![alt text](image-2.png)

Pearson, Spearman and R² values are high across all three pairs, indicating strong agreement on overall expression levels. Some divergence is visible between STAR and the Salmon methods, likely reflecting differences in how multimappers are handled.

## Log File Sumamries  

### Salmon-Map Log files 

| Method | n samples | Avg total fragments | Avg reads in eq. classes | Avg disc. (align score) | Avg disc. (fragment score) | Avg disc. (decoy) | Avg mapping rate (%) | Targets | Decoys |
|---|---|---|---|---|---|---|---|---|---|
| Salmon-map | 70 | 43,044,599 | 15,698,474 | 24,226,429 | 9,643,825 | 6,633,015 | 36.2 | 37,320 | 938 |

The mapping rate of 36% is low. Out of 43M fragments, 9.6M failed alignment scoring, 6.6M mapped best to a genome decoy and were discarded, and 15.7M were assigned to equivalence classes for probabilistic quantification. The low rate could  reflects the fragmented and repetitive nature of the C. maculatus genome rather than a methodological issue, since featureCounts also discards a large fraction of reads.  
Mappings discarded due to alignment score = represents individual alignments that fail the scoring threshold.  
Fragments discarded due to alignment score = fragments that failed all their mappings.  
Fragments discarded due to decoy matching = they mapped best to a decoy.   

### STAR Log files 

| Method | n samples | Avg input reads | Avg uniquely mapped | Avg multi-mapped | Unique (%) | Multi-mapped (%) | Too many loci | Too many loci (%) |
|---|---|---|---|---|---|---|---|---|
| STAR | 70 | 43,044,599 | 26,447,066 | 4,964,529 | 61.1 | 11.6 | 7,722,021 | 18.2 |

61.1% of reads mapped uniquely. 11.6% mapped to multiple loci and are reported in the BAM. 18.2% mapped to too many loci and were discarded. 

### FeatureCounts Log files 

| Method | Avg assigned | Avg unassigned (no feature) | Avg unassigned (ambiguous) |
|---|---|---|---|
| featureCounts | 22,497,598 | 42,087,519 | 17,832,473 |

22.5M reads assigned, 42.1M discarded as they did not overlap an annotated exon, and 17.8M overlapping multiple features resolved fractionally.

### STAR (transcriptome BAM for Salmon-align)

| Method | n samples | Avg input reads | Avg uniquely mapped | Avg multi-mapped | Unique (%) | Multi-mapped (%) | Too many loci | Too many loci (%) |
|---|---|---|---|---|---|---|---|---|
| STAR transcriptBAM | 70 | 43,044,599 | 23,277,654 | 4,459,149 | 53.8 | 10.5 | 7,606,532 | 17.9 |

53.8% mapped uniquely, 17.9% discarded as too multi-mapped. The 10.5% that multi-mapped were passed to Salmon for probabilistic assignment.

### Salmon-Align

| Method | n samples | Avg total mapped | Avg uniquely mapped | Avg multi-mapped | Unique (%) | Multi (%) | Avg reads in eq. classes |
|---|---|---|---|---|---|---|---|
| Salmon-align | 70 | 11,685,928 | 10,533,184 | 1,152,744 | 89.9 | 10.1 | 11,558,587 |

Of the reads passed by STAR, 89.9% were assigned uniquely by Salmon and 10.1% were ambiguous and modelled probabilistically.

### Method Selection 






For raw statistical power STAR + featureCounts has the most kept and significant transcripts. However, some of this apparent power can be inflated due to the fractional assignment of the paralogous transcripts.

Salmon-Align has the lowest DE counts, lowest PC1 variance and is the most conservative approach as the most ambiguous reads were already filtered out during genome alignment.

Salmon-mapping considers all fragments, uses probabilistic modelling for assigning ambiguous reads, doesn't discard ambiguous reads before the probibalistic modelling and still has a relatively high DE signal despite being more conservative. The low mapping rate is still a bit concerning, but seeing as featureCounts also discards a majority of reads, this might be a genome issue? 

In our case, since we are interested in paralogs downstream, i think salmon map might be the best here. Its not as conservative as salmon align. It still has a high biological signal while using probibalistic modelling and decoy filterining to handle the ambigous reads instead of just fractal counting. We might trade some statistical raw power for more confidence in the origin on the reads.

Later downstream the data from STAR could be used as a comparison and see if the results are comparable. 

# Paralog analyses  
I used preexisting data from Orthofinder that include 14 beetle species + Drosophyla melanogaster as an outgroup. This assigns each gene to an Orthogroup (OG) and a Hierarchical Orthogroup (HOG). OG groups all genes across all species that descend from a single gene in the root of the tree.  
HOGs provide more context by breaking down OG at each node of the species tree, creating subfamilies within the same OG. A single OG might split into several HOGs depending on when speciation events occurred.  

As we are only looking at genes/transcripts in C.mac for these analyses, using Hierarchical Orthologous Groups (HOGs) will be synonomous with gene family, as each HOG represents all transcripts in C_mac that share a common ancestor at a given node in the tree.  

In the OrthoFinder analysis, the proteome was isoform filtered to retain only the canonical isoform per gene. Consequently, in my following analyses based on this data on things like HOG and age rank annotations, i will only have one transcript per gene represented. And no isoform-level random effets are necessary in the mixed models. 

# Gene family size and sex bias (**HOG_size_sex_bias.ipynb**)
First I wanted to investigate if the larger the gene family, the more sex baised paralogs it will have, as having more copies would lead me to suspect that there are more opportunities for sex-specific neofunctionalization to occur.  
First i removed all the transcript that doest belong to a HOG, resulting in 15.648 transcripts.  
I then defined three different sets of gene family sizes:  
1. Genome wide  
The size of the gene families are defined by all of the transcripts in the annotation.  
2. Mapped / not expressed   
These are all transcripts that managed to be mapped to the transcriptome by Salmon. These stem from the dataset prefilter_df which comes after tximport in the script **salmon_map_dominance_consistent_script.R**, but before the expression filtering (≥5 counts in ≥5 samples) and the differential expression analysis.  
3. Expressed  
These are the sizes defined by the transcripts in the result dataset, post filtering and DE-analysis.  

I did this since only looking at the size definitions of expressed transcript would not reflect reality. For example, if we only looked at the Expressed set, the largest gene family size appears to be 25. But that is based only on the "visible" transcripts. But If we look at the entire genome from full_annotation, the largest gene family has 115 members (if all of those were expressed). This indicates that a majority of members in the large gene families are not expressed or lowly expressed in abdominal tissues.  
In the dataset each expressed transcript has a column for each size definition. 

Regardless of size definition however, a majority of the expressed transcripts belong to very small gene families, and the median size is consistently 1. For most plots later we filter out size 1.  

| Level | Total transcripts | Transcripts with Gene family | Transcripts without Gene family | # Gene families | Max family size | Median family size |
|---|---|---|---|---|---|---|
| Genome | 37,988 | 31,507 | 6,481 | 14,563 | 115 | 1 |
| Mapped | 36,382 | 30,041 | 6,341 | 14,563 | 114 | 1 |
| Expressed | 17,574 | 15,648 | 1,926 | 10,850 | 25 | 1 |

### Density plot of size definition comparison 

![alt text](image-55.png)  

log10 transformed:  
![alt text](image-56.png)  

### Size category expression differences 

If we look at the full genome definitions, what is the proportion of mapped and expressed transcripts?  
![alt text](image-59.png)  

Here we can see that consistently aross sizes, <50% of transcripts are not expressed, but a majority was mapped by salmon (but too lowly expressed). 

### lFC within gene families vs. gene family size. 
To investigate sex bias and gene family size i started by looking at the transcript logFoldChange (male vs female) within each gene family, against the size of each gene family.  
First with the expressed size definition:   
![alt text](image-11.png)  

And with the full genome size definition, with :  
![alt text](image-58.png)  

This indicate however that the strength of sex-bias decreases the larger the gene family size is. Does indicate that male bias is more common than female bias, especially for larger sizes where female-bias is rarer. By looking at the full genome plot its made even more clear that sex-bias decreases with the size of the gene family.   

### Bias direction among biased transcripts  
To investigare further; within each gene family, which proportion is male biased at each gene family size? Due to the difference in sample sizes I added wilson confidence intervals (preferred over standard CI as proportions reaches 0 at some points).  
![alt text](image-9.png)   
Out of the transcripts that are sex-biased, more are male-biased. 

What is the proportion of male biased, female biased and unbiased transcripts at each gene family size?  
![alt text](image-10.png)

A general trend can be observed with most transcripts being unbiased, followed by male-biased and female-biased.  
The same two plots but with the full genome size definitions (a bit cluttered):  
![alt text](image-46.png)  
![alt text](image-47.png)

### Variance vs. gene family size  
Variance for all transcripts that belong to a gene family of a certain size. Here only looking at the expressed size definitions.   
![alt text](image-48.png)  
Expression differene decreases as the family size increases. But for HOG sizes larger than 15 we have very small sample sizes so variance is noisy and unreliable as small sample size makes variance unstable.

### Variance within each gene family 
If we instead look at each gene family individually:  
![alt text](image-49.png)  
Size 1 filtered out. This plot again shows a trend that the direction of expression in the larger gene families are more "in agreement" in direction than the smaller sizes. The outlier in red is N0.HOG0000055. 

### Within gene family directional bias  
Next I categorized the gene families based on the expression direction of the transcripts they contain. These being all unbiased, all male-biased, all female-biased and the combinations of the three:  
![alt text](image-50.png)  

If we divide the multiple category bars into the overall proportion transcript bias within the categories, we see that its around 50% bias and unbiased in the double categories, and in the all three category its 40.2% unbiased, 36.3% male-biased and23.5% female-biased.    
![alt text](image-68.png)  

We can visualize this further with the categories and the gene family sizes they include:  
Expressed:    
![alt text](image-52.png) 

Full genome:  
![alt text](image-51.png)

The largest families are found in the two sex biased + unbiased categories.  

### Gene family size sex-bias category proportions  

Expressed:  
![alt text](image-53.png)  

Full genome:  
![alt text](image-60.png)  

# Paralog ancestry - Age rank analyses
I will also look at the relative branch lengths within each HOG as a indicator of the age of each paralog using mixed models.  
I downloaded SpeciesTree_rooted.txt, SpeciesTree_rooted_node_labels.txt, Duplications.tsv and SpeciesTree_Gene_Duplications_0.5_Support.txt from the OrthoFinder output.  

By uploading this file into Phylo.io I could visually inspect the phylogeny, see screenshot below:  

![alt text](image-69.png)    

In order to add the age rank i updated the script (**create_full_annotation.R**).   
I loaded SpeciesTree_rooted_node_labels.txt in order to create a node age tree with the depths of the nodes to the roots and the branch lengths of the nodes that lead to C_maculatus in the phylogeny (N0, N1, N2, N5, N8, N10, N12, N13 and the tip node C_maculatus_proteinfaste_TE_filtered) and i ranked the nodes from oldest (N0 = 1) to youngest (C_maculatus_proteinfaste_TE_filtered = 9) in decending order, based on the node's depth from the root.  

![alt text](image-71.png) 

Then i load Duplications.tsv to get information about when duplication events occured, and which genes resulted from each duplicaiton, and i filtered for only C_mac transcripts and duplications that occured on the branches leading to C_Mac in the phylogeny. 
I merged this with the age_rank table and exported this as dup_long which contains all the duplication events in the C_mac lineage.  

But this counts every single duplication event, so one transcript can appear many times (nested in the family tree history). 

I therefore assigned each transcript its most recent duplication as we want the latest point where the copy became independant, into dup_most_recent. I created a transcript_id column and joined this with the rest of the full_annotation. 

I imported dup_long.csv and full_annotation_with_age.csv into VSCode for plotting: 

Age diversity of all duplication events:  

![alt text](image-36.png)  

A majority of all duplication events occurred in C_mac.

I did a combined plot with 4 dataset levels:  

1. Full annotation: All of the C-mac transcripts most recent copy.  

2. Prefilter transcripts: I then merged the annotation with the mapping results from Salmon-map (**salmon_map_dominance_consistent_script.R**).  
This was done after tximport but before the expression filter, exporting a new file "prefilter_transcripts_annotated.csv". Here we get the age_rank information about all of the transcripts mapped by Salmon.  

3. DE results: Post-DE-analysis. Includes all expressed transcripts (5 counts in 5 samples).  

4. Only the significnatly DE transcripts (M v F)

| Dataset         | Info                                                        | Transcript count | Transcript count with age info |
|-----------------|-------------------------------------------------------------|------------------|-------------------------------|
| dup_long        | Total duplication events (multiple per transcript)          | 92,785           | 92,785 (100%)                 |
| full_annotation | Most recent copy. All annotated transcripts                 | 37,988           | 21,484 (56.6%)                |
| prefilter_df    | Mapped transcripts by Salmon                                | 36,382           | 20,467 (56.3%)                |
| de_results      | Expressed transcripts (≥5 counts in ≥5 samples)             | 17,574           | 9,796 (55.7%)                 |
| significant_de  | Significantly DE transcripts (padj < 0.05, \|LFC\| > 1)     | 7,270            | 4,319 (59.4%)                 |

dup_long includes all of the duplication events in C_mac only, hence 100%. Full_annotation is defined from the .gff. In each level only around 50-60% of transcripts from the GFF have had a duplication event in the C_mac lineage. The remaining ~40% of transcripts could still belong to a HOG, but have no duplication events in the C_mac lineage. Or they could not belong to a HOG at all, and are novel C_Mac specific genes with no homologs, TE-derived or very fragmented in the assembly. 

![alt text](image-37.png)  

The common thread in all of these levels is that the majority of transcripts most recent copy originate in C.mac, consistent with the high number of predicted genes in the species, and the repetetive nature of the genome. There is a noteable drop in the amount of age 9 transcripts from mapped to expressed. Do these belong to the outstadningly large amount of annotated "genes" in the C_mac genome?  

I wanted to see if we see the same pattern in related species. So i tried to replicate the grey bar for the full genome duplication ages for the species C_chienensis (sister species), A_obtectus (Bruchnae) and T_castaneum (Tenebrionidae).

I softlinked their isoform_filtered annotation files and cleaned them up with the script **clean_gff_related_species.sh**. 

![alt text](image-72.png)

## Age rank and sex bias

### Box plots:  
Here i tried to replicate Milenas plot of the significantly DE expressed / sex-biased transcripts in each age rank, split between male and female biased (padj < 0.05, |log2FC| > 1, where log2FC > 1 = male-biased, log2FC <1 = female-biased).

4319 transcripts with age rank were plotted. (2951 significantlly DE transcripts are missing age rank. Out of these 2158 belong to a gene family, and 793 does not).

![alt text](image-38.png)  

Higer male bias across all ranks. 

Age ranks 4(N5), 5(N8), 6(N10) interesting here? Big variance, high mean and high median?  
N5 is where R_ferrugineus and D_ponderosae splits off  
N8 is where A_obtectus splits off  
N10 is where B_siliquastri splits off  

### Proportion sex-bias in age ranks  
Here i tried to replicate Milenas plot of the proportions of unbiased, male-biased or female-biased transcripts within each age rank.  
9796 expressed transcripts were plotted, significance was defined by padj < 0.05, |log2FC| > 1, where log2FC > 1 = male-biased, log2FC <1 = female-biased. 

![alt text](image-39.png)  

In all ages, most are unbiased, followed by male-biased and last female-biased. There is an increase in proportionally male-biased transcripts in ages 4(N5), 5(N8) and 6(N10) once again, hinting at intermediate age increase in male-bias. 

## Age rank and gene family size  
 
| Dataset         | Total transcripts | With age_rank        | With HOG             | With both (plotted)     |
|-----------------|------------------|-----------------------|-----------------------|---------------------------|
| full_annotation | 37,988           | 21,484 (56.6%)        | 31,507 (82.9%)        | 21,467 (56.5%)            |
| prefilter_df    | 36,382           | 20,467 (56.3%)        | 30,041 (82.6%)        | 20,450 (56.2%)            |
| de_results      | 17,574           | 9,796 (55.7%)         | 15,648 (89.0%)        | 9,785 (55.7%)             |
| sig             | 7,270            | 4,319 (59.4%)         | 6,472 (89.0%)         | 4,314 (59.3%)             |

A higher proportion of transcripts belong to a HOG than have age rank information, making age rank the limiting factor. 

### Within gene family transcript age rank diversity:   
Looking at the number of different ages within each gene family:

![alt text](image-40.png)

Most gene families only have transcripts of the same age, but there are still candidates of gene families with 2 or 3 different ages in the significantlly DE dataset.  

### Top 10 most age diverse 
Top 10 most age rank diverse gene families. This is all gene families with more than 3 DE transcripts in them. The age range might indicate an continously expanding gene family (darker), with duplication events occuring throughout time rather in shorter bursts (lighter) (max age rank - min age rank +1). 

![alt text](image-42.png) 

I might come back to this later and investigate those families.  

### Age rank distribution by gene family size  
Plotted the 9785 transcripts with both HOG and age rank info. 
![alt text](image-66.png)   

### Age rank vs. Gene family size with sex-bias info.  

Same sex bias definition. With the gene family sizes as the y-axis.    
Genome sizes:  
![alt text](image-61.png)
Expressed sizes:   
![alt text](image-62.png) 

Again indicates that most gene family expansions are happening within C_mac. 

I also tried to do it with box plots separated into unbiased, male-biased and female-biased:    
Genome sizes:  
![alt text](image-63.png)
Expressed sizes:  
![alt text](image-64.png)   
These are all of the expressed trasncripts, not just the significant ones. 

And as proporitonal bars, where i grouped together the genome-defined sizes into singeltons, small, medium and large families to limit this to 4 plots.  
![alt text](image-67.png) 

# Chromosomal location analyses 

Sex Chromosomal contigs: 

{ X : ['utg000057l_1','utg000114l_1','utg000139l_1','utg000191l_1','utg000326l_1','utg000359l_1','utg000532l_1','utg000602l_1'],
  Y : ['utg000322l_1','utg 000312c_1','utg 000610l_1','utg 001235l_1']}


From the paper i filtered contigs longer than 100kbp to identify "true" autosomes. The locations X, Y, A (autosome) and U (unassigned/belongs to the short contigs). 

#transcripts on each before expression filtering:  
    A     U     X     Y 
28651  7176  1827   334

After expression filtering:  
    A     U     X     Y 
16630   148   719    76

# Mixed model analyses 
I use lme4 on the post-DESeq2 data as a two-stage approach. The DESeq2 step identifies which transcripts are sex-biased while controlling for underlying genotypes. The lme4 step then asks whether evolutionary age predicts the pattern and magnitude of that sex bias across transcripts, with gene families as a random effect to account for the non-independence of transcripts within the same family.
Two types of models are used throughout:

lmer for continuous responses (log2FoldChange, |log2FoldChange|)  
glmer for binary responses (is this gene sex-biased: yes/no)  

Mixed models have fixed and random effects. Fixed effects are the variable we wish to estimate and interpret. Such as age rank or the sex-bias class. 

Random effects are variables we want to control for but not interpret directly. Here the random effect is always (1 | HOG) which is a random intercept that gives each gene family its own baseline expression bias. This will account for the fact that paralogs within the same family are not independent observations. 

I use REML (Restricted Maximum Likelihood) to fit the variance components in lmer models. This gives unbiased estimates of the random effect variance compared to regular maximum likelihood. For likelihood ratio tests (LRT) comparing two models, both models must be fitted with regular maximum likelihood so they are on the same scale.

Two tools are used to interpret results beyond the coefficient tables:

drop1 performs a likelihood ratio test by dropping one term at a time from the full model and testing whether removing it significantly worsens the fit. This is a Type III test, meaning it tests each term after accounting for all others. This lets me see if an interaction or main effect is justified in the model. 

emmeans (estimated marginal means) computes the predicted value of the response for each combination of factor levels, averaged over all other predictors in the model. This is what the model thinks the mean response would be for a typical gene in each group. Contrasts then compare pairs of these predicted means directly.

## Dataset Summary

| | Count |
|---|---|
| Total transcripts loaded | 17,574 |
| Transcripts with age rank | 9,796 |
| Transcripts with HOG | 15,648 |
| **Model data** (age rank + HOG + logFC) | **9,785** |
| HOGs represented | 5,048 |
| Unique gene_ids | 9,785 |

OrthoFinder was run on an isoform-filtered proteome (one protein per gene) so eac htranscript maps to only one gene and one HOG. There is no isofrm nesting to wory about, so (1 | HOG) should be sufficient for the random  effects. 

### Dataset Issues
Two structural properties of the data shape the analysis throughout:  
24.8% of genes are singletons (gene families of size 1). These genes have no family members to compare against, so they have no relative age within a family. They are excluded from any within-family analysis.  

89% of gene families have all paralogs at the same age rank. This means most gene family expansions happened in a single burst duplication event at one node. Relative age is only biologically meaningful in the 554 families where paralogs arose at different nodes. This is why i did the two-dataset design used in Parts 2 and 3.

I divided the analyses in three parts. The raw model outputs can be found in the corresponding mixed_model_partX_results.txt files.


# Part 1: Absolute Age Models 
Does the absolute evolutionary age of a gene predict how sex-biased it is? Age is deined using 9 discrete ranks (1 = oldest, shared with Diptera; 9 = youngest, C. mac specific).


## Model 1 - Direction of sex bias  
**Formula:** `log2FoldChange ~ age_rank + (1 | HOG)`  
**Response:** log2FoldChange (M vs F). Positive = male-biased, Negative = female-biased.  

lmer treats all log2FC estimates as equally certain regarless of expression level. To adress this, each model is run twice. Once unweighted (all transcripts equal) and once weighted by `1/lfcSE²`, (inverse variance from DESeq2), so that transcripts with precise logFC estimates contribute more.

Three age parameterisations are compared:

* M1a: Continuous scaled age rank. One unit = 1 SD increase toward younger (SD of age rank = 2.9 units). Assumes equal spacing between ranks.  

* M1b: Age as 9 discrete factor levels. Rank 1 (N0) is the reference. Each coefficient shows how much a given rank differs from the oldest. Captures non-linearity but costs 7 extra degrees of freedom.  

* M1c: Real cumulative branch lengths from the tree root, scaled. More biologically honest than arbitrary rank spacing.  

AIC is used to compare the three parameterisations within the same weighting scheme (the samaller AIC the better the model). ΔAIC > 4 is considered meaningful.

### AIC Comparison - Unweighted:

| Model | Parameterisation  | df | AIC   |
| ----- | ----------------- | -- | ----- |
| M1a   | Continuous rank   | 4  | 42514 |
| M1b   | Factor (9 levels) | 11 | 42491 |
| M1c   | Node depth        | 4  | 42516 |

M1b wins by ΔAIC = 23. The non-linear node-specific pattern is real.

### M1b Results (Unweighted, best model)
Each coefficient shows the estimated log2FC at that rank relative to rank 1 (N0).

| Age rank | Node  | Estimate vs rank 1 | p-value  | Sig. |
| -------- | ----- | ------------------ | -------- | ---- |
| 1 (ref)  | N0    | +0.697 (intercept) | 1.87e-07 | ***  |
| 2        | N1    | -0.279             | 0.061    | .    |
| 3        | N2    | -0.219             | 0.177    |      |
| 4        | N5    | +0.261             | 0.238    |      |
| 5        | N8    | +0.302             | 0.048    | *    |
| 6        | N10   | +0.478             | 0.015    | *    |
| 7        | N12   | -0.064             | 0.686    |      |
| 8        | N13   | -0.441             | 0.061    | .    |
| 9        | C_mac | +0.208             | 0.140    |      |

Ranks 5 (N8) and 6 (N10) are significantly more male-biased than rank 1. These are bruchid-specific nodes (similar to what we have seen in earlier plots).

### AIC Comparison - Weighted:
| Model | Parameterisation  | df | AIC   |
| ----- | ----------------- | -- | ----- |
| M1a   | Continuous rank   | 4  | 39145 |
| M1b   | Factor (9 levels) | 11 | 39151 |
| M1c   | Node depth        | 4  | 39146 |

When expression precision is accounted for, the non-linearity disappears and the simpler continuous model wins.  

The weighted M1b tells the opposite story: ranks 2, 3, 4 and 7 are significantly less male-biased than rank 1, meaning the oldest genes are most male-biased when expression precision is accounted for. The discrepancy shows that the unweighted male-bias signal at young ranks is partly driven by lowly expressed genes with noisy fold-change estimates.

![alt text](image-73.png)

## Model M1_prob: Probability of Sex Bias by Node
**Formula**: `glmer(sex_biased_bin ~ age_rank_factor + (1 | HOG), family = binomial)`
**Response:** Binary (1 = sex-biased by DESeq2 thresholds |logFC| > 1 and padj < 0.05, 0 = unbiased).

As a complement to M1b, asking wheter genes are biased at all, rather than how biased they are at each node. drop1 tests wheter the factor as a whole significantly improves the fit. EMMs give the predicted probability of sex bias at each node (using type = "response"). 

![alt text](image-74.png)
The plot indicates that genes at the bruchid-specific nodes, are not only more male-biased in direction (from M1b) but also significantly more likely to be classified as sex-biased generally. Ranks 7-9 show some decline back toward baseline, suggesting the elevated sex bias is specific to these intermediate bruchid nodes rather than a continous increase with younger age.


## Model 2: Magnitude of Sex Bias
**Formula:** `lmer(abs(log2FoldChange) ~ age + (1 | HOG))`    
**Response:** Absolute log2FoldChange - measure of how sex-biased a gene is regardless of which sex is higher.

### M2a: Continuous Age, Random Intercept
**Formula:** `lmer(abs_logFC ~ age_rank_scaled + (1 | HOG)`
Run on the full dataset, and on HOG >= 2 (multi-member families only) for comparison.

| Dataset           | Slope (age_rank_scaled) | p-value  |
| ----------------- | ----------------------- | -------- |
| Full (n=9,785)    | +0.162                  | 2.45e-13 |
| HOG ≥ 2 (n=7,361) | +0.094                  | 0.001    |

Younger genes show significantly hihger sex bias magnitude in both datasets. The effect is smaller once singletons are excluded, suggesting singletons in young families inflate the signal somewhat.

### M2a_factor: Factor Age, Random Intercept
**Formula:** `lmer(abs(log2FoldChange) ~ age_rank_factor + (1 | HOG))`

The factor version shows the predicted magnitude at each specific node rather than an average slope. This is the magnitude equivalent of M1b to see which nodes stand out. EMMs give the predicted |log2FC| at each rank. drop1 tests whether the factor as a whole significantly predicts magnitude.

![alt text](image-75.png)

The predicted magnitude of sex bias follow a similar pattern. Dramatic increase at rank 4. Interestingly, this decreases after rank 4, indicating that intermediate-age families at the bruchid node drive the sex-biased signal.

### M2b: Random Slope Model
**Formula:** `lmer(abs_logFC ~ age_rank_scaled + (1 + age_rank_scaled | HOG))`

This allows each gene family to have its own age-magnitude relationship, instead of assuming the slope is the same for all families. A random slope requires at least 3 observations per gene family to estimate reliably (only 2 points will give zero residual degrees of freedom within the HOG), so this is run on a HOG >= 3 dataset (3,725 transcripts, 806 HOGs).

I used the bobyqa optimizer because random slope models are harder to fit and the default optimizer can have trouble near convergence boundaries.

To test whether the random slope is justified, M2a and M2b are compared using a likelihood ratio test on the same HOG >= 3 dataset. Both models must be fitted with maximum likelihood (REML = FALSE) for a valid LRT comparison.

**LRT result:** p < 2.2e-16, so the random slope model fits significantly better.

However, the average slope across all families is not significant (p = 0.14). Within-family age trajectories vary alot across families, and those trajectories cancel out at the population level.  
The intercept-slope correlation is +0.30, meaning gene families with higher sex bias also show a steeper age gradient within themselves. The age-sex bias connection is concentrated in families that already have sex-biased expression.

| Model    | n transcripts | n HOGs | Fixed slope (age) | p    | Slope variance | Intercept–slope r |
| -------- | ------------- | ------ | ----------------- | ---- | -------------- | ----------------- |
| M2a_filt | 3,372         | 5,806  | 0.077             | 0.14 | -              | -                 |
| M2b_filt | 3,372         | 5,806  | 0.077             | 0.14 | 0.537          | +0.30             |

Model 3: Does the Age-Magnitude Relationship Differ by Sex Bias Class?
Formula: lmer(abs(log2FoldChange) ~ sex_bias * age_rank_scaled + (1 | HOG))

### Model 3: Does the Age-Magnitude Relationship Differ by Sex Bias Class?
**Formula:** `lmer(abs_logFC ~ sex_bias * age_rank_scaled + (1 | HOG)`
**Response:** |log2FoldChange|  
**Predictors:** sex_bias (three levels: unbiased / male-biased / female-biased) crossed with continuous age rank. Unbiased is the reference group. "Unbiased" means the gene did not pass |logFC| > 1 and padj < 0.05

A significant interaction means the age-magnitude slope differs between sex bias classes. drop1 tests the interaction term. EMMs predict the magnitude for each class at low (-1 SD), average (0) and high (+1 SD) age, allowing direct comparison of the slopes.

| Dataset           | Interaction F | p-value |
| ----------------- | ------------- | ------- |
| Full (n=9,785)    | 4.04          | 0.018   |
| HOG ≥ 2 (n=7,361) | 8.47          | 0.0002  |

The interaction term is significant in both datasets, stronger in the multi-member family. 

**Range from EMM table:** 
| Sex bias class | At -1 SD (older) | At average age | At +1 SD (younger) | Change | Interpretation |
|----------------|-----------------|----------------|--------------------|--------|----------------|
| Unbiased       | 0.58            | 0.73           | 0.88               | +0.30  | Steep positive slope |
| Male-biased    | 2.81            | 2.86           | 2.91               | +0.10  | Nearly flat |
| Female-biased  | 2.11            | 2.18           | 2.26               | +0.15  | Nearly flat |

Predicted values are |log2FoldChange|. Age is on a standardised scale where -1 SD correspond to roughly ranks 3-4 and +1 SD roughly ranks 6-7.  
This indicate that the age-magnitude gradient from M2a is driven by unbiased genes. Once a gene is sex-biased, its magnitude stays roughly constant regardless of evolutionary age. Younger unbiased genes show higher fold-change variation, but it doesnt cross the significance threshold. 

### Model 3b Factor Age Version (Supplementary)
**Formula:** `lmer(abs(log2FoldChange) ~ sex_bias * age_rank_factor + (1 | HOG))`

Using age_rank_factor shows which nodes drives interaction, and drop1 gives an interation p=2.8e-16, much stronger than the continous version. 

EMMs gives the predicted magnitude for each sex bias class at each node. Two nodes stand out, a spike in male magnitude at N5 a reversal at N13. All other nodes show male-biased genes exceeding female-biased genes in magnitude. 

| Node         | Male abs(log2FC) | Female abs(log2FC) | Notable                                         |
| ------------ | ------------- | --------------- | ----------------------------------------------- |
| Rank 4 (N5)  | 4.12          | 2.24            | Male magnitude nearly double any other node     |
| Rank 8 (N13) | 2.32          | 3.56            | Only node where female exceeds male (p = 0.014) |

![alt text](image-76.png)

## Covariate controlled models 

To test wether the age effects are confounded by gene length, mean expression level or gene family size. Young genes tend to be shorter and lower expressed, and larger families have more observed paralogs per bin. If the age signals simply reflect these properties, they are not informative.

Covariates added:

* log_gene_length_scaled:  
log-transformed transcript length from GFF coordinates (end - start), scaled to mean=0 SD=1.
* log_baseMean_scaled:  
log-transformed DESeq2 mean normalized count across all samples, scaled.
* hog_size_genome_scaled:  
gene family size at the genome level, scaled.

### M2a_cov
**Formula:** `lmer(abs_logFC ~ age_rank_scaled + log_gene_length + log_baseMean + hog_size + (1 | HOG))`

Table from drop1: 
| Term                   | F value | p-value |
| ---------------------- | ------- | ------- |
| age_rank_scaled        | 47.5    | 5.8e-12 |
| log_gene_length_scaled | 35.1    | 3.3e-09 |
| log_baseMean_scaled    | 49.1    | 2.5e-12 |
| hog_size_genome_scaled | 2.5     | 0.115   |

The age effect strenghtens after controlling for these factors. Gene length and expression are predictors of magnitude, but they do not explain the age factor. 

### M3_cov
**Formula:** `lmer(abs_logFC ~ sex_bias * age_rank_scaled + log_baseMean + (1 | HOG))`

Table from drop1:
| Term                     | F value | p-value   |
| ------------------------ | ------- | --------- |
| log_baseMean_scaled      | 108.6   | < 2.2e-16 |
| sex_bias:age_rank_scaled | 6.15    | 0.002     |

The sex bias x age interaction holds after controlling for expression level. 

# Part 2: Relative age models 
Does a paralogs relative position within its gene familys duplication history predict sex bias and does this depend on the family itself? 

**New definitions:**  
**Family age:** The age rank of the oldest paralog in the HOG. This defines when the gene family first appeared in the phylogeny. Binned into three categories using tertiles of the distribution across all gene families:

| Category     | Age rank range                                  |
| ------------ | ----------------------------------------------- |
| Old          | Ranks 1–3 (shared with Coleoptera or older)     |
| Intermediate | Ranks 4–7 (bruchid-specific nodes)              |
| Young        | Ranks 8–9 (Callosobruchus or *C. mac* specific) |

**Relative age:** The percentile rank of each paralog's age_rank within its gene family. The oldest copy gets 0, the youngest gets 1. Copies that share the same age_rank receive the same percentile rank. If all copies in a family share the same age_rank (a burst duplication at one node), all of them receive 0, because none is older or younger than any other.

Binned into two levels: old (percentile rank <= 0.5) and young (> 0.5). Three bins were structurally impossible: a young family cannot have a "youngest paralog" relative to an older within-family copy, because paralogs in a young family most of the time originate at the same recent node. The young-family x young-paralog interaction cell was always empty.

**Two analysis datasets**  
**Analysis A:** Full dataset excluding singletons (n = 7,361, 2,624 gene families). Singletons have no relative age within a family.  
**Analysis B:** Restricted to the 554 gene families where paralogs arose at different nodes (sequential duplication). Only in these families is relative age biologically meaningful. 89% of gene families have all paralogs at the same age rank and are excluded here.

### Model MS1: Probability of Sex Bias
**Formula:** `glmer(sex_biased_bin ~ rel_age_cat * family_age_2lev + (1 | HOG), family = binomial)`

Family age is collapsed to two levels (old vs non-old) for the glmer, because the young-family x young-paralog cell is more or less empty and causes the model to fail with 3 levels.

drop1 tests the interaction term. EMMs give predicted probabilities per group. Contrasts compare young vs old paralogs within each family age class.

MS1_A showed convergence warnings due to large HOG random effect variance relative to the small fixed effects. I used allFit to confirm these are numerical artifacts and not a real problem: all five optimizers gave essentially the same fixed effect estimates (maximum range across optimizers: 0.028) and log-likelihood (range: 0.23).

| Model                     | Interaction p-value          |
| ------------------------- | ---------------------------- |
| MS1_A (full)              | 0.44                         |
| MS1_B (age-variable gene families) | 0.44                |

The overall interaction is not significant. However, the marginal contrast in non-old families is significant:

| Family age       | Contrast (old vs young paralog) | Odds ratio | p-value |
| ---------------- | ------------------------------- | ---------- | ------- |
| Old families     | old vs young                    | 0.82       | 0.35    |
| Non-old families | old vs young                    | 0.57       | 0.043   |

In non-old families, young paralogs have significantly higher odds of being sex-biased than old paralogs (OR=0.57 means old paralogs have 43% lower odds). This does not hold for the old families, but as mentioned the overall interaction is non-significant. 

### Model MS2: Magnitude of Sex Bias
**Formula:** `lmer(abs_logFC ~ rel_age_cat * family_age_cat + (1 | HOG))`

Family age uses 3 levels here. The lmer is more tolerant of sparse cells than the glmer, and the young-family x young-paralog cell, while very sparse, is still estimable for a continuous response.

| Model                     | Interaction F | p-value |
| ------------------------- | ------------- | ------- |
| MS2_A (full)              | 5.11          | 0.006   |
| MS2_B (age-variable HOGs) | 3.75          | 0.024   |


The interaction is significant in both datasets. 

Table from EMM(MS2_A):
| Relative age | Family age | Predicted abs(log2FC) |
|-------------|------------|----------------------|
| Old paralog | Old family | 1.38 |
| Young paralog | Old family | 1.55 (+0.17, p=0.15) |
| Old paralog | Intermediate family | 2.04 |
| Young paralog | Intermediate family | 1.76 (-0.29, p=0.05) |
| Old paralog | Young family | 1.65 |
| Young paralog | Young family | excluded (near-empty cell) |

The highest sex bias magnitude is in old paralogs from intermediate-age families. Within those families, young paralogs are actually less biased than old ones - the opposite of our prediction.  In old families the direction fits the prediction (young slightly higher) but the difference is not significant.

This connects to the M1b and M3b findings: intermediate-age families correspond to the bruchid-specific nodes (ranks 4-7) where sex bias is elevated overall. Within those families, it is the established older copies driving the bias, not the most recently derived ones.

![alt text](image-77.png)

### Model MS2_A: Covariate Control
**Formula:** `lmer(abs_logFC ~ rel_age_cat * family_age_cat + hog_size_genome_scaled + log_baseMean_scaled + (1 | HOG)`

Controlling for family size and expression.  

| Term                       | F value | p-value |
| -------------------------- | ------- | ------- |
| hog_size_genome_scaled     | 2.58    | 0.109   |
| log_baseMean_scaled        | 34.7    | 4.1e-09 |
| rel_age_cat:family_age_cat | 5.17    | 0.006   |

The interaction survives covariate control. Expression level matters independently (higher expressed genes are more sex-biased), but it does not explain why old paralogs in intermediate-age families show elevated sex bias magnitude.

# Part 3: Duplication Offset Models
Instead of relative age i use a more direct approach. Since with relative age, the copy labeled as "youngest" within a family doesnt capture how long since the gamilies origin it originated. With offsets I wanted to capture different the founder copies are to the derived copies, and use this evolutionary distance instead. 

Duplication offset = paralog age_rank − min(age_rank in HOG)

* Offset = 0: founder copy — arose when the gene family was first established
* Offset > 0: derived copy — arose through a later duplication event at a more recent node

This uses the phylogenetic distance from the families node of origin, rather than relative position among siblings. All founder copies get offset 0 regardless of how many there are.

This has some other limitations though.
Three offset bins (founder / early / late) were structurally impossible. High offsets can only occur in old families, since a family that originated at node 8 cannot produce a copy five nodes later. The interaction table would have the same empty cell problem as Part 2. Binary (founder vs derived) is the only approach that worked. 

Both models MO1 and MO2 use Analysis B (554 HOGs with sequential duplication), since the offset concept is only meaningful where paralogs arose at different nodes.

Relative offset rescales the offset to 0-1 within each family's own evolutionary window: rel_offset = dup_offset / (9 - family_origin). This corrects for the fact that old families have had more time to accumulate high-offset copies. Used in the sensitivity model MO2, restricted to old and intermediate families because young families can only have rel_offset values of 0 or 1 (no continuous gradient to fit).

### Model MO1: Binary Offset
**Formula (magnitude):** `lmer(abs_logFC ~ offset_bin1 + family_age_2lev + (1 | HOG)`  
**Formula (magnitude):** `glmer(sex_biased_bin ~ offset_bin1 + family_age_2lev + (1 | HOG)`

No interaction term - the interaction between offset and family age is not estimable (empty cells, same structural reason as the three-bin problem).

**EMMs - magnitude:**
| Copy type | Family age | Predicted abs(log2FC) |
|-----------|------------|----------------------|
| Founder | Old | 1.50 |
| Derived | Old | 1.60 (+0.11) |
| Founder | Non-old | 1.98 |
| Derived | Non-old | 2.08 (+0.11) |

**EMMs - probability**
| Copy type | Family age | Predicted P(sex-biased) |
| --------- | ---------- | ----------------------- |
| Founder   | Old        | 0.462                   |
| Derived   | Old        | 0.509                   |
| Founder   | Non-old    | 0.510                   |
| Derived   | Non-old    | 0.557                   |

**drop1 results:**
| Term            | Model    | F / LRT    | p-value |
| --------------- | -------- | ---------- | ------- |
| offset_bin1     | MO1_mag  | F = 2.34   | 0.126   |
| family_age_2lev | MO1_mag  | F = 9.71   | 0.002   |
| offset_bin1     | MO1_prob | LRT = 2.93 | 0.087   |
| family_age_2lev | MO1_prob | LRT = 1.39 | 0.239   |

Family age is the dominant predictor of magnitude. Derived copies are 0.11 |log2FC| higher than founders on average, in the predicted direction, but this does not reach significance. For probability, offset is borderline (p = 0.087): founders have 17% lower odds of being sex-biased than derived copies (OR = 0.83).

![alt text](image-78.png)

### Model MO2: Continuous Relative Offset (Sensitivity)
**Formula (magnitude):** `lmer(abs_logFC ~ rel_offset + family_age_2lev + (1 | HOG)`  
**Formula (magnitude):** `glmer(sex_biased_bin ~ rel_offset + family_age_2lev + (1 | HOG)`

Restricted to old and intermediate families (n = 2,406, 543 HOGs). Young families are excluded because their rel_offset is either 0 or 1 with nothing in between, making a continuous slope uninterpretable.

| Term            | Model    | F / LRT    | p-value |
| --------------- | -------- | ---------- | ------- |
| rel_offset      | MO2_mag  | F = 2.15   | 0.143   |
| family_age_2lev | MO2_mag  | F = 8.40   | 0.004   |
| rel_offset      | MO2_prob | LRT = 3.31 | 0.069   |
| family_age_2lev | MO2_prob | LRT = 0.80 | 0.371   |

The continuous and binary approaches agree: moving from founder (rel_offset = 0) to the latest possible derived copy (rel_offset = 1) adds ~ 0.11 to |log2FC| (p = 0.14). The fact that the effect size and direction are identical in MO1 and MO2 rules out the binning choice as the reason for non-significance.

### MO1_mag_cov: Covariate Control
**Formula (magnitude):** `lmer(abs_logFC ~ offset_bin1 + family_age_2lev + hog_size_genome_scaled + (1 | HOG)`  
| Term                   | F value | p-value |
| ---------------------- | ------- | ------- |
| offset_bin1            | 2.33    | 0.127   |
| family_age_2lev        | 9.69    | 0.002   |
| hog_size_genome_scaled | 0.0004  | 0.983   |

Family size has essentially zero effect (F < 0.001). The offset result is not a confounder artifact - the effect is simply small.

### Summary of results
| Question                                                  | Model              | Key result                                                                                                                                       |
| --------------------------------------------------------- | ------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| Does age predict direction of sex bias?                   | M1b unweighted     | N8 (rank 5) and N10 (rank 6) significantly more male-biased. Signal disappears when weighted.                                                    |
| Does age predict magnitude of sex bias?                   | M2a_full           | Younger genes more sex-biased in magnitude (p = 2.5e-13). Effect smaller in multi-member families.                                               |
| Does the magnitude–age relationship vary across families? | M2b                | Random slope justified (p < 2.2e-16). Average slope flat across families. High-bias families show steeper within-family age gradients.           |
| Does the age–magnitude slope differ by sex bias class?    | M3                 | Yes (p = 0.002–0.018). Unbiased genes drive the age gradient. Sex-biased genes are flat across ages.                                             |
| Which nodes show extreme sex bias magnitude?              | M3b                | Rank 4 (N5): extreme male magnitude (4.12). Rank 8 (N13): only node where female exceeds male.                                                   |
| Does relative paralog age predict P(sex-biased)?          | MS1                | No significant interaction. Young paralogs in non-old families have borderline higher odds (OR = 0.57, p = 0.043).                               |
| Does relative paralog age × family age predict magnitude? | MS2                | Yes (p = 0.006 and 0.024 in the two datasets). Old paralogs in intermediate-age families show highest magnitude.                      |
| Does founder vs derived status predict sex bias?          | MO1                | Directional trend (derived +0.11).                                                                                                               |
| Do covariates explain the age effects?                    | M2a_cov, MS2_A_cov | No. Age effects survive and strengthen after controlling for gene length, expression, and family size.                                           |

**Overall conclusion**  
Absolute family age is the strongest and most consistent predictor of sex bias. Gene families that originated at bruchid-specific nodes (particularly N5 and N10 for direction, and the intermediate age category for magnitude) show the highest sex bias. The within-family age effects (relative age, duplication offset) show directional trends consistent with the prediction that derived copies evolve more sex bias, but these effects are smaller and less consistent than the family age effect. Once a gene is sex-biased, its magnitude does not continue to increase with younger evolutionary age. It is unbiased genes that show the strongest age-magnitude gradient, possibly suggesting that evolutionary timing determines whether a gene crosses the threshold into sex-biased expression, but not how extreme that bias becomes once it does.