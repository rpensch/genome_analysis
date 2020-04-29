#Quality control of raw sequences

module load bioinfo-tools
module load FastQC

#Illumina sequences with FastQC - genomic data

fastqc data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz data/genomic_data/raw_data/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz --outdir=analyses/1_preprocessing/fastqc_raw

#Illumina sequences with FastQC - transcriptomics data

fastqc data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797972_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

fastqc data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797973_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

fastqc data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_BH/trim_paired_ERR1797974_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

fastqc data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797969_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

fastqc data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797970_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

fastqc data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_1.fastq.gz data/transcriptomics_data/RNA-Seq_Serum/trim_paired_ERR1797971_pass_2.fastq.gz --outdir=analyses/1_preprocessing/fastqc_raw

