module load bioinfo-tools
module load sratools

fastq-dump --split-files --skip-technical -O data/snp_extra --gzip SRR3306347

#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 02:00:00 
#SBATCH -J snp_calling
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  

# Load modules 

module load bioinfo-tools
module load bwa
module load samtools

# Your commands

bwa index analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta
bwa mem analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
data/snp_extra/SRR3306347_1.fastq.gz \
data/snp_extra/SRR3306347_2.fastq.gz | samtools view -b | samtools sort -o analyses/5_snp_calling/rna_mapping_snp_extra_paired.bam -O BAM -T rna_mapping_snp_exra_tmp

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
module load samtools/0.1.19 
module load python/2.7.15


# Your commands

samtools faidx analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta

samtools mpileup -g -f analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797972.bam \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797973.bam \
analyses/6_rna_mapping/rna_mapping_BH_paired_ERR1797974.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797969.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797970.bam \
analyses/6_rna_mapping/rna_mapping_Serum_paired_ERR1797971.bam \
analyses/5_snpcalling/rna_mapping_snp_extra_paired.bam > my-raw.bcf

bcftools view -bvcg my-raw.bcf > my-var.bcf
bcftools view my-var.bcf | vcfutils.pl varFilter - > snp_calling.vcf



module load bioinfo-tools bcftools

grep -v '#' snp_calling.vcf | wc -l

bcftools view -v snps snp_calling.vcf > snps.vcf
grep -v '#' snps.vcf | wc -l

bcftools view -v indels snp_calling.vcf > indels.vcf
grep -v '#' snps.vcf | wc -l

bcftools view -i 'QUAL>10' snp_calling.vcf | grep -v '#' | wc -l
bcftools view -i 'QUAL>50' snp_calling.vcf | grep -v '#' | wc -l
bcftools view -i 'QUAL>100' snp_calling.vcf | grep -v '#' | wc -l
