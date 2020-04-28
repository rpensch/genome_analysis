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

# Your commands
bwa index analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_2.fastq.gz > analyses/6_rna_mapping/aln_BH_paired_ERR1797972.pe.sam

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_2.fastq.gz > analyses/6_rna_mapping/aln_BH_paired_ERR1797974.pe.sam

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_2.fastq.gz > analyses/6_rna_mapping/aln_BH_paired_ERR1797973.pe.sam

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_2.fastq.gz > analyses/6_rna_mapping/aln_Serum_paired_ERR1797969.pe.sam

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_2.fastq.gz > analyses/6_rna_mapping/aln_Serum_paired_ERR1797970.pe.sam

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_1.fastq.gz \
data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_2.fastq.gz > analyses/6_rna_mapping/aln_Serum_paired_ERR1797971.pe.sam

