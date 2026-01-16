#!/bin/bash
#SBATCH -A uppmax2025-2-148
#SBATCH -p core
#SBATCH -n 16
#SBATCH -t 24:00:00
#SBATCH -J star_index_gene_counts_dominance
#SBATCH --output=%x.%j.out

module load bioinfo-tools
module load star/2.7.11a

GENOME_FA="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/C_maculatus_assembly.fna"
GTF="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/C_maculatus_annotation_unfiltered_fixed.gtf"
GENOME_DIR="/proj/naiss2023-6-65/Sebastian/data/reference_genome_c_mac/star_index"
RNA_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/fastp_multiqc"
OUT_DIR="/proj/naiss2023-6-65/Sebastian/data/rna_data_kaufmann_2024_dominance/results/star_gene_counts"

mkdir -p "$GENOME_DIR" "$OUT_DIR"

#Create the STAR genome index
STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir "$GENOME_DIR" \
     --genomeFastaFiles "$GENOME_FA" \
     --sjdbGTFfile "$GTF" \
     --sjdbOverhang 149 #max read length 150-1

#Align reads and produce gene counts
for R1 in "$RNA_DIR"/*_1.fastq.gz; do
  R2="${R1/_1.fastq.gz/_2.fastq.gz}"
  SAMPLE="$(basename "$R1" "_1.fastq.gz")"
  STAR --runThreadN 20 \
       --genomeDir "$GENOME_DIR" \
       --readFilesIn "$R1" "$R2" \
       --readFilesCommand zcat \
       --outSAMtype BAM SortedByCoordinate \
       --quantMode GeneCounts \
       --twopassMode Basic \
       --outFilterMultimapNmax 20 \
       --limitBAMsortRAM 20000000000 \
       --outFileNamePrefix "$OUT_DIR/${SAMPLE}_"
done
