#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 03:00:00 
#SBATCH -J plasmid_identification 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 

module load bioinfo-tools 
module load spades

# Your commands

spades.py \
-1 data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz \
-2 data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz \
--plasmid --nanopore data/genomic_data/raw_data/Nanopore/E745_all.fasta.gz \
-k 33,55 -o analyses/8_plasmid_identification