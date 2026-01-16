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

Analyses to be run:  
Run salmon-align on Uppmax, continue with pipeline in R and python

Orthofinder analysis in python 

Below follows the pipeline used. 

# Datasets

## Reference Genome and Annotaions 

The reference genome comes from the paper "Y-Linked Copy Number Polymorphism of Target of Rapamycin Is Associated with Sexual Size Dimorphism in Seed Beetles" by Kaufmann et. al 2023. Male virgin C_maculatus abdominal tissue samples from the Lomé population. I am using the the small male Y haplotype assembly from the paper, as it is more continuous, and the small Y haplotype is the most abundant haplotype in the population.

Most scripts are titeled "_unfiltered" as they are based the non-isoform filtered annotation files, as i want to conserv that information. 

I created a symbolic link to braker_proteins.aa in order to run eggNog to get functional annotation, which was combined with the structural annotation to create the "full annotation" in R (**run_eggnog.sh**).  

In R, I first converted the non-isoform filtered gtf file to a gff3 file. This structural annotation file was merged with the results from eggnog as well as the results from OrthoFinder (N0.tsv file) to create a more comprehensive structural and functional annotation.

## RNA Dataset 1

Dataset 1 is based on "Sex-Specific Dominance of Gene Expression in Seed Beetles" by Kaufmann et.al 2024. 

10 gen full sib inbred lines.  
Lomè population of C_maculatus.  
Three pairwise crosses of six isogenic lines.  
Six homozygous (no signs of inbreeding depression), three heterozygote.  
Used each isogenic line as both maternal and paternal = reciprocal crosses (4 genotypes per cross).  
RNA seq: They sequenced 3 replicates for each of the 4 genotypes and sexes, giving 24 samples per cross and 72 samples in total.  
The data was downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB70958&o=acc_s%3Aa with the script (**download_PRJEB70968.sh**).   

Run FastQC and MultiQC (**run-fastqc_multiqc.sh**).
Run fastp for trimming (**run_fastp_multiqc.sh**).
Ran fastqc and multiqc again to confirm improvements (**run_fastqc_multiqc_post_trim.sh**).

### Metadata

Label corrected metadata for dataset 1 is found in **dominance_meta_corrected_outlier_corrected.xlsx** and **dominance_meta_corrected_outlier_corrected.csv**.  

After PCA visualization one sample (ERR12383283) was changed from male to female due to clustering and suspected misidentification. Two samples (ERR12383297 male and ERR12383303 male) were removed due to ambiguous sexes. Three remaining samples are suspected of ambigous sex as they stray from the respective clusters, but are kept (ERR12383254 female, ERR12383278 male, ERR12383310 male). 


# Mapping methods

Three different mapping methods are used and will be compared. Salmons mapping based mode/quasi mapping, STAR with featureCounts and the combined method of Salmons alignment based mode + STARs .bam files. For the main part of the project RNA-Seq data from dataset 1 was used. Analyses were run on the transcript level rather than gene level. 

Information about the two salmon modes are found here: https://salmon.readthedocs.io/en/latest/salmon.html#

## Salmon-mapping based mode (Quasi-mapping) 

In this mode, salmon needs a transcriptome and a decoy file. 

A transcriptome was created from the reference genome + the unfiltered annotation using gffread (**create_transcript_unfiltered.sh**) 

I used the whole genome decoys approach (where genome sequences themselves serve as decoys for the transcripts), and generated a gentrome (all transcripts first, then the genome/decoy sequences), and a decoy .txt file (which includes the names/headers of the genome sequences) (**generate_decoys_whole_genome_unfiltered.sh**).

The salmon index was created from the gentrome and the decoy files (**salmon_index_whole_genome_unfiltered.sh**).

Then i ran salmon for mapping on dataset 1, using the flags:  
--qcBias  
which corrects for GC-content during quantification,  
--seqBias  
which corrects for sequence specific bias where fragments starting with certain motifs might get preferential sequencing,  
--ValidateMappings  
which is the selective alignment mode (which is now default),  
(**salmon_mapping_whole_genome_unfiltered.sh**). 

The alignments were transferred to R, where I;  
-ran txiimport on the transcript level,  
-filtered on ≥3 mean counts per sample in each sex,  
-DESeq2 for DE analysis based on male vs female,  
-used vst for count normalization with variance stabilization.  
(**salmon_map_dominance_unfiltered_WG_decoys.R**)  

I combined the results with the structural and functional annotations and imported them to Visual Studio Code for plotting (PCA and Volcano Plot).  
(**salmon_map_unfiltered_plotting_transcript.ipynb**)

## STAR (with featureCounts)

Created a STAR index with the flags:

--sjdbGTFfile,  
--sjdbOverhang 149 (the max read length -1)  
for making it splice junction aware. 

Then aligned reads and counts using the flags:

--outSAMtype BAM SortedByCoordinate  
sorted based on genomic location,  
--quantMode GeneCounts  
performs gene-level quantification,  
--twopassMode Basic 
can discover novel junctions not in the annotation,  
--outFilterMultimapNmax 20  
max # of multiple alignment locations per read. Multimapped reads are randomly assigned to one location. 

(**star_alignment_dominance.sh** (includes indexing), **star_alignment_dominance_continuation.sh**, **star_alignment_dominance_continuation_2.sh**).

Picard was used to mark read ruplicates by setting a flag in the bam files based on identical start positions (can be removed later), and samtools were used to index the bam files and ease downstream analyses (**run_picard_samtools.sh**). 

Subread featureCounts were used to summarize counts per gene/transcript. Four different modes were used:  
Mode 1: Standard counting
No multimapping, summarize counts by gene, used for differential expression analysis with the flags:  
-p -B -C   
(paired-end fragments/read pairs, both must map, does not count chimeric fragments)  
-g "gene_id"  
-t "exon"  
-s 2   
(reverse stranded)  

2.54 billion gene counts
Found 35489 genes 

Mode 2: gene level with multimapping. 
Still the flags -g "gene_id" and -t "exon", but also adding the -M and --fraction flags for multimapping and fractional counting (splitting reads between the multipe targets). 
2.67 billion gene counts (~5% are multimapped) 

Mode 3: transcript level counting. 
To compare with the dominance paper. They said "summarizing exons per transcript". Using -f flag to count on the exon level and then later sum to get transcript level, and -g for grouping by transcript_id instead. 
1.442 billion transcript counts. 

Mode 4: transcript level with multimapping 

-M \
--fraction \
-f \
-g "transcript_id" \
-t "exon" \

(**run_subread_featurecounts.sh**, **run_subread_featureCounts_transcript_multi.sh**)

This created the files: 	
gene_counts_standard.txt  
gene_counts_multimappers.txt  
transcript_counts.txt  
Transcript_counts_multimappers.txt

These files were imported to RStudio, where I: 
-Aggregate to transcript level by summing exon counts,  
-loaded the multimapped transcript counts,  
-filtered on ≥3 mean counts per sample in each sex,  
-DESeq2 for DE analysis based on male vs female,  
-used vst for count normalization with variance stabilization.  
(**star_DE_analysis.R**)  

I combined the results with the structural and functional annotations and imported them to Visual Studio Code for plotting (PCA and Volcano Plot).  
(**STAR_plotting.ipynb**)

## Salmon-alingment based mode 

First had to rerun STAR to be compatible with salmon and transcript alignment files.  
Still splice junction aware with twopassMode Basic.  
--quantMode TranscriptomeSAM (Outputs a BAM file aligned to transcript sequences. Salmon requires reads to be mapped to transcriptome coordinates). 
--outSAMtype BAM SortedByCoordinate 
(Sorted by reference coordinates. I saw conflicting ideas, but this should be compatible with Salmon)
--outSAMattributes NH HI AS nM XS GX GN
(Tags needed for Salmon.  
Nr alignmeds for a read, alignment index for multimappers, alignment score, nr of mismatches, strand information, Gene ID, Gene name.)  
--outFilterMultimapNmax 20
(keeps up to 20 alignments per read)  
--winAnchorMultimapNmax 100  
(How many "anchor points" per window. Max 100 regions)  

(**star_transcriptome_for_salmon_1.sh**, **star_transcriptome_for_salmon_2.sh**, **star_transcriptome_for_salmon_3.sh**)

Im waiting for the UPPMAX node to be up and running before I can run Salmon and downstream R and python analyses. 

# Mapping software comparison

## Correlation tests  
Average baseMean per sex per gene between the three softwares. 

# Paralog analyses  
I used the read data from ...  

## HOG size and sex bias 
First I looked at the mean logFoldChange (male vs female) within each Hierarchical Orthogroup (HOG) against the size of each HOG, with the hypothesis that the more paralogs each HOG have, the higher the average logFC will be. 

## Paralog ancestry and 
I will also look at the relative branch lengths within each HOG as a indicator of the age of each paralog using a mixed model.  



