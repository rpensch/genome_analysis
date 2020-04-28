#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 06:00:00 
#SBATCH -J rna_mapping 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  

# Load modules 
module load bioinfo-tools
module load bwa
module load samtools

# Your commands
bwa index analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_2.fastq.gz | samtools view -b | samtools sort -o rna_mapping_BH_paired_ERR1797972.bam

