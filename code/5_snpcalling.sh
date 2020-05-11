module load bioinfo-tools
module load sratools

fastq-dump --split-files --skip-technical -O data/snp_extra --gzip SRR3306347

bwa index analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/snp_extra/...1 \
data/snp_extra/...2 | samtools view -b | samtools sort -o analyses/5_snp_calling/rna_mapping_snp_extra_paired.bam

#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 04:00:00 
#SBATCH -J snp_calling
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  

# Load modules 

module load bioinfo-tools
module load samtools
module load python/2.7.15
module load bcftools

# Your commands

samtools faidx analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

samtools mpileup -g -f analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797972.bam \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797973.bam \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797974.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797969.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797970.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797971.bam | bcftools view -bvcg | vcfutils.pl varFilter - > snp_calling.vcf