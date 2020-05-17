#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 05:00:00 
#SBATCH -J synteny_comparison
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  

# Load modules 

module load bioinfo-tools
module load blast 

# Your commands

blastn –db nt –query analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta –out analyses/8_synteny/homology_search_results.out -outfmt "6 qseqid sseqid sacc saccvers sscinames evalue pident" -num_threads 2
blastn –db nt –query analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta -evalue 1 –out analyses/8_synteny/homology_search_results1.out  -outfmt "6 qseqid sseqid sacc saccvers sscinames evalue pident" -num_threads 2

blastn -query analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta -subject analyses/8_synteny/GCF_000391485.2_ASM39148v2_genomic.fna -out analyses/8_synteny/efaecium_efaecalis_blast.out -outfmt 6 -num_threads 2

wc -l homology_search_results.out
wc -l homology_search_results1.out