module load bioinfo-tools 
module load MUMmer

mv analyses/2_genome_assembly/1_efaecium_canu/efaecium.contigs.fasta \
analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta #rename canu assembly file

nucmer --mum -p 1_assembly_reference data/reference/GCF_001720945.1_ASM172094v1_genomic.fna analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

mummerplot -p 1_assembly_mummerplot_reference -R data/reference/GCF_001720945.1_ASM172094v1_genomic.fna -l -Q analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta --png 1_assembly_reference.delta 

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
data/reference/GCF_001720945.1_ASM172094v1_genomic.fna \
-o analyses/3_assembly_evaluation/quast \
--gene-finding
