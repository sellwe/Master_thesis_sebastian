# Master_thesis_sebastian

Below follows the pipeline used. 

# Datasets

## Reference Genome

The genome comes from the paper "Y-Linked Copy Number Polymorphism of Target of Rapamycin Is Associated with Sexual Size Dimorphism in Seed Beetles" by Kaufmann et. al 2023. I am using the the small male Y haplotype assembly from the paper, as it is more continuous, and the small Y haplotype is the most abundant haplotype in the population.

Most scripts are titeled "_unfiltered" as they are based the non-isoform filtered annotation files, as i want to conserv that information. 

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

