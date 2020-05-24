setwd("C:/Users/Raphaela/Documents/UU/GA/genome_analysis/results/8_rna_readcount")

bh1 <- read.table("htseq_readcount_BH_paired_ERR1797972.txt", row.names = 1, col.names = c('gene', 'ERR1797972'))
bh2 <- read.table("htseq_readcount_BH_paired_ERR1797973.txt", row.names = 1, col.names = c('gene', 'ERR1797973'))
bh3 <- read.table("htseq_readcount_BH_paired_ERR1797974.txt", row.names = 1, col.names = c('gene', 'ERR1797974'))
serum1 <- read.table("htseq_readcount_Serum_paired_ERR1797969.txt", row.names = 1, col.names = c('gene', 'ERR1797969'))
serum2 <- read.table("htseq_readcount_Serum_paired_ERR1797970.txt", row.names = 1, col.names = c('gene', 'ERR1797970'))
serum3 <- read.table("htseq_readcount_Serum_paired_ERR1797971.txt", row.names = 1, col.names = c('gene', 'ERR1797971'))

countData <- data.frame(bh1, bh2, bh3, serum1, serum2, serum3)

n<-dim(countData)[1]
countData<-countData[1:(n-5),]


hist(countData[,1], main="Histogram of ERR1797972 Count Data", xlab="Counts")

log2_counts <- log2(countData+1)

boxplot(log2_counts, main="Distribution of log2-normalized Read Counts")
hist(log2_counts[,1], main="Histogram of ERR1797972 log2-normalized Count Data", xlab="Counts")
