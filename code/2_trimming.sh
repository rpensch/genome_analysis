#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 02:00:00 
#SBATCH -J trimming
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com 

# Load modules 

module load bioinfo-tools
module load trimmomatic

# Your commands

java -jar $TRIMMOMATIC_HOME/trimmomatic.jar PE data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz analyses/1_preprocessing/trimming/illumina_trimmed_paired_1 analyses/1_preprocessing/trimming/illumina_trimmed_unpaired_1 analyses/1_preprocessing/trimming/illumina_trimmed_paired_2 analyses/1_preprocessing/trimming/illumina_trimmed_unpaired_2 TRAILING:28 MINLEN:75

fastqc analyses/1_preprocessing/trimming/illumina_trimmed_paired_1 analyses/1_preprocessing/trimming/illumina_trimmed_paired_2 --outdir=analyses/1_preprocessing/fastqc_trimmed