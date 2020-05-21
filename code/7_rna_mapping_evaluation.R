setwd('C:/Users/Raphaela/Documents/UU/GA/genome_analysis/results/8_rna_readcount/')

bh1 <- read.table("htseq_readcount_BH_paired_ERR1797972.txt", row.names = 1, col.names = c('gene', 'ERR1797972'))
bh2 <- read.table("htseq_readcount_BH_paired_ERR1797973.txt", row.names = 1, col.names = c('gene', 'ERR1797973'))
bh3 <- read.table("htseq_readcount_BH_paired_ERR1797974.txt", row.names = 1, col.names = c('gene', 'ERR1797974'))
serum1 <- read.table("htseq_readcount_Serum_paired_ERR1797969.txt", row.names = 1, col.names = c('gene', 'ERR1797969'))
serum2 <- read.table("htseq_readcount_Serum_paired_ERR1797970.txt", row.names = 1, col.names = c('gene', 'ERR1797970'))
serum3 <- read.table("htseq_readcount_Serum_paired_ERR1797971.txt", row.names = 1, col.names = c('gene', 'ERR1797971'))

countData <- data.frame(bh1, bh2, bh3, serum1, serum2, serum3)

total_reads_per_file <- colSums(countData)
not_aligned_per_file <- countData["__not_aligned",]
mapped_reads_per_file <- total_reads_per_file - not_aligned_per_file

total_reads <- sum(total_reads_per_file)
total_not_aligned <- sum(not_aligned_per_file)
total_mapped_reads <- sum(mapped_reads_per_file)

perc_mapped_reads <- total_mapped_reads * 100 / total_reads
perc_mapped_reads

n<-dim(countData)[1]
reads_on_features<-countData[1:(n-5),]

total_reads_on_features <- sum(colSums(reads_on_features))
perc_reads_on_features <- total_reads_on_features * 100 / total_reads
perc_reads_on_features

total_no_feature <- sum(countData["__no_feature",])
total_no_feature
