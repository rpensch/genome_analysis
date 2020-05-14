module load bioinfo-tools 
module load MUMmer

mv analyses/2_genome_assembly/2_efaecium_spades/contigs.fasta \
analyses/2_genome_assembly/2_efaecium_spades/efaecium_spades.contigs.fasta #rename spades assembly file

nucmer --mum -p 1_assembly_nucmer analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta analyses/2_genome_assembly/2_efaecium_spades/efaecium_spades.contigs.fasta

mummerplot -p 1_assembly_mummerplot --png 1_assembly_nucmer.delta

#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 01:00:00 
#SBATCH -J quast 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 
module load bioinfo-tools 
module load quast

# Your commands

quast.py analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
analyses/2_genome_assembly/2_efaecium_spades/efaecium_spades.contigs.fasta \
-o analyses/3_assembly_evaluation/quast \
--gene-finding