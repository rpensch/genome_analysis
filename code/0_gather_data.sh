#Gathering the data in the working directory by creating soft links

#Genomics data
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.1.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.2.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131023_233432_42174_c100519312550000001823081209281335_s1_X0.3.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.1.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.2.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/PacBio/m131024_200535_42174_c100563672550000001823084212221342_s1_p0.3.subreads.fastq.gz ../data/raw_data/genomic_data/PacBio

ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz ../data/raw_data/genomic_data/Illumina
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz ../data/raw_data/genomic_data/Illumina

ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/genomics_data/Nanopore/E745_all.fasta.gz ../data/genomic_data/raw_data/Nanopore

#Transcriptomics data
ln -s /proj/g2020008/nobackup/private/1_Zhang_2017/transcriptomics_data/ ../data/