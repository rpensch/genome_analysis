#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 08:00:00 
#SBATCH -J rna_readcount
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  

# Load modules 
module load bioinfo-tools
module load htseq

# Your commands

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797972.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_BH_paired_ERR1797972.txt

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797974.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_BH_paired_ERR1797974.txt

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797973.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_BH_paired_ERR1797973.txt

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797969.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_Serum_paired_ERR1797969.txt

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797970.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_Serum_paired_ERR1797970.txt

htseq-count -f bam -r pos -t CDS analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797971.bam analyses/4_annotation/efaecium_annotation_noseq.gff > analyses/7_rna_readcount/htseq_readcount_Serum_paired_ERR1797971.txt