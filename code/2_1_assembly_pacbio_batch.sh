#!/bin/bash -l   
#SBATCH -A g2020008 
#SBATCH -p core 
#SBATCH -n 2 
#SBATCH -t 06:00:00 
#SBATCH -J pensch_pacbio_assembly_1 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user rar.pensch@gmail.com  
 
# Load modules 
module load bioinfo-tools 
module load canu  

# Your commands
canu \
 -p efaecium -d efaecium_canu \
 genomeSize=3.2m \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.1.subreads.fastq.gz \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.1.subreads.fastq.gz \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.3.subreads.fastq.gz \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.1.subreads.fastq.gz \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.2.subreads.fastq.gz \
 -pacbio-raw data/genomic_data/raw_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.3.subreads.fastq.gz
 


