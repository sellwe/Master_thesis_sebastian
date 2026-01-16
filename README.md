# Master_thesis_sebastian

Scripts that will be run:  
Convert gtf to gff3  
Run salmon-align  


Below follows the pipeline used. 

# Datasets

## Reference Genome and Annotaions 

The genome comes from the paper "Y-Linked Copy Number Polymorphism of Target of Rapamycin Is Associated with Sexual Size Dimorphism in Seed Beetles" by Kaufmann et. al 2023. I am using the the small male Y haplotype assembly from the paper, as it is more continuous, and the small Y haplotype is the most abundant haplotype in the population.

Most scripts are titeled "_unfiltered" as they are based the non-isoform filtered annotation files, as i want to conserv that information. 

I created a symbolic link to braker_proteins.aa in order to run eggNog to get functional annotation, which was combined with the structural annotation to create the "full annotation" in R (**run_eggnog.sh**).  

*Data from Orthofinder will also be combined ...*

## RNA Dataset 1

Dataset 1 is based on "Sex-Specific Dominance of Gene Expression in Seed Beetles" by Kaufmann et.al 2024. The data was downloaded from https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJEB70958&o=acc_s%3Aa with the script (**download_PRJEB70968.sh**). 

Run FastQC and MultiQC (**run-fastqc_multiqc.sh**).
Run fastp for trimming (**run_fastp_multiqc.sh**).
Ran fastqc and multiqc again to confirm improvements (**run_fastqc_multiqc_post_trim.sh**).

# Mapping methods

## Salmon-mapping based mode (Quasi-mapping) 

In this mode, salmon needs a transcriptome and a decoy file. 

A transcriptome was created from the reference genome using the unfiltered annotation using gffread (**create_transcript_unfiltered.sh**) 

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
I combined the results with the structural and functional annotations and imported them to Visual Studio Code for plotting (PCA and Volcano Plot). 

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

(**star_alignment_dominance.sh (includes indexing)**, **star_alignment_dominance_continuation.sh**, **star_alignment_dominance_continuation_2.sh**)

Picard was used to mark read ruplicates by setting a flag in the bam files based on identical start positions (can be removed later), and samtools were used to index the bam files and ease downstream analyses (**run_picard_samtools.sh**). 

Subread featureCounts were used to summarize counts per gene/transcript). Four different modes were used:  
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
I combined the results with the structural and functional annotations and imported them to Visual Studio Code for plotting (PCA and Volcano Plot). 

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

# Orthofinder work 



