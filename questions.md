# Student Manual Questions and Answers

## Reads Quality Control

### What is the structure of a FASTQ file?
### How is the quality of the data stored in the FASTQ files? How are paired reads identified?
### How is the quality of your data?
### What can generate the issues you observe in your data? Can these cause any problems during subsequent analyses?

## Reads Preprocessing

### How many reads have been discarded after trimming?
### How can this affect your future analyses and results?
### How is the quality of your data after trimming?
### What do the LEADING, TRAILING and SLIDINGWINDOW options do?

## Genome Assembly 

### What information can you get from the plots and reports given by the assembler (if you get any)?
### What intermediate steps generate informative output about the assembly?
### How many contigs do you expect? How many do you obtain?
### What is the difference between a ‘contig’ and a ‘unitig’?
### What is the difference between a ‘contig’ and a ‘scaffold’?
### What are the k-mers? What k-mer(s) should you use? What are the problems and benefits of choosing a small kmer? And a big k-mer?
### Some assemblers can include a read-correction step before doing the assembly. What is this step doing?
### How different do different assemblers perform for the same data?
### Can you see any other letter appart from AGTC in your assembly? If so, what are those?

## Assembly Evaluation

### What do measures like N50, N90, etc. mean? How can they help you evaluate the quality of your assembly? Which measure is the best to summarize the quality of the assembly (N50, number of ORFs, completeness, total size, longest contig ...)
### How does your assembly compare with the reference assembly? What can have caused the differences?
### Why do you think your assembly is better/worse than the public one?
### When running metaQuast for a metagenome, it may happen that very few contigs map back to the reference genomes. Is this expected? Does that mean your assembly is bad? Why?

## Annotation

### What types of features are detected by the software? Which ones are more reliable a priori?
### How many features of each kind are detected in your contigs? Do you detect the same number of features as the authors? How do they differ?
### Why is it more difficult to do the functional annotation in eukaryotic genomes?
### How many genes are annotated as ‘hypothetical protein’? Why is that so? How would you tackle that problem?
### How can you evaluate the quality of the obtained functional annotation?
### How comparable are the results obtained from two different structural annotation softwares?

## Homology Search 

### How relevant is the output format that you choose?
### How do the resulting hits vary when you change the minimum e-value?
### How is the alignment score calculated?
### How important is the number of threads when you blast against a database, or against a particular sequence?

## Mapping 

### What percentage of your reads map back to your contigs? Why do you think that is?
### What potential issues can cause mRNA reads not to map properly to genes in the chromosome? Do you expect this to differ between prokaryotic and eukaryotic projects?
### What percentage of reads map to genes?
### How many reads do not map to genes? What does that mean? How does that relate to the type of sequencing data you are mapping?
### What do you interpret from your read coverage differences across the genome?
### Do you see big differences between replicates?

## Post-mapping analyses - SNP calling

### What is the structure of a SAM file, and how does it relate to a BAM file?
### What is the structure of vcf and bcf files?
### How many SNPs and INDELs do you get?
### How is the quality of those variants?
### What is the difference between the variant quality, the mapping quality and the fastq quality?
### How are these variantes distributed along the genome?

## Post-mapping analyses - Read counting

### What is the distribution of the counts per gene? Are most genes expressed? How many counts would indicate that a gene is expressed?
### In the metagenomics project, the data doesn’t offer enough statistical power for a differential expression analysis. Why not? What can you still tell from the data only from the read counts?

## Expression Analyses

### If your expression results differ from those in the published article, why could it be?
### How do the different samples and replicates cluster together?
### What effect and implications has the p###value selection in the expression results?
### What is the q-value and how does it differ from the p-value? Which one should you use to determine if the result is statistically significant?
### Do you need a normalization step? What would you normalize against? Does DESeq do it?
### What would you do to increase the statistical power of your expression analysis?
### In the metagenomics project, the data doesn’t offer enough statistical power for a differential expression analysis. Why not? What can you still tell from the data?