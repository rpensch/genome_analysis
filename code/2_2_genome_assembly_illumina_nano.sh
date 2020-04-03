#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 03:00:00 
#SBATCH -J spades_assembly_1 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 
module load bioinfo-tools 
module load spades

# Your commands

spades.py --isolate \
-1 data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
-2 data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
--nanopore data/genomic_data/raw_data/Nanopore/E745_all.fasta.gz \
-k -o analyses/2_genome_assembly

#question: should I use --plasmid when there is 1 chromosome as well? 
#question: what kmer size should I choose?
