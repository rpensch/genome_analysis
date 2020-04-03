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
mv analyses/2_genome_assembly/1_efaecium_pacbio/efaecium.contigs.fasta \
analyses/2_genome_assembly/1_efaecium_pacbio/efaecium_pacbio.contigs.fasta #rename pacbio assembly file

mv 



quast.py analyses/2_genome_assembly/1_efaecium_pacbio/efaecium_pacbio.contigs.fasta \
test_data/contigs_2.fasta \
-o analyses/3_assembly_evaluation/quast \
--gene-finding

