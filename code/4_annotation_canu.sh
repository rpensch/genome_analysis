#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 00:30:00 
#SBATCH -J annotation 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 
module load bioinfo-tools 
module load prokka 

# Your commands
prokka --outdir analyses/4_annotation --prefix efaecium_annotation \
--genus Enterococcus --species faecium --strain E745 --gram pos \
analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta