#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 01:00:00 
#SBATCH -J mummerplot 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 
module load bioinfo-tools 
module load MUMmer

# Your commands
nucmer --mum -p 1_assembly_nucmer \  
analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
analyses/2_genome_assembly/2_efaecium_spades/efaecium_spades.contigs.fasta \


mummerplot -c -p 1_assembly_mummerplot -S 1_assembly_nucmer.delta
mummerplot -c -p 2_assembly_mummerplot 2_assembly_nucmer.delta
mummerplot -c -p 3_assembly_mummerplot -S -R -Q 1_assembly_nucmer.delta
mummerplot -c -p 4_assembly_mummerplot -R -Q 2_assembly_nucmer.delta


