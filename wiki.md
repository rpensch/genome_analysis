# Genome Analysis
## Project Plan
### Introduction
In the course of this project I am aiming to reproduce the results of the following research paper: 

"RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum” (2017) by Xinglin Zhang et al. 

*Enterococcus faecium* is a bacterium that normally occupies the human digestive system but can often cause infections of the bloodstream in hospitalized patients. This is dangerous since some strains have acquired resistance against the most common antibiotics. Surviving and thriving in such suboptimal and nutrient-poor conditions as they are in the human blood is a challenge and the mechanisms that make it possible for *E. faecium* are to be analyzed in this study. In a first step, the question whether and how gene expression differs in human serum and conditions of a conventional nutrient-rich growth medium is of interest. This is analyzed by examining the transcriptomes of *E. faecium* grown in both environments. Later, genes that promote growth in human serum are identified by comparing transposon sequencing (Tn-seq) data of transposon mutant libraries of the bacterium from both serum and plain growth medium. 

The data that is available to try and answer these questions is for one sequencing data from Illumina HiSeq 100 bp paired-end sequencing, Pacific Biosciences RS II SMRT technology long-read sequencing and Oxford Nanopore Technologies MinION with R7 flowcell chemistry to assemble the genome of *Enterococcus faecium*. Transcriptome data is available from 100 bp paired-end sequencing with the Illumina HiSeq 2500 platform and Tn-seq data from 50 nt single-end sequencing with Illumina HiSeq 2500. 

### Methods
The steps and analyses required for carrying out this project successfully will be described below in roughly the order they will be performed.

#### Read Quality Control:
Before starting the project it is important to make sure the quality of the sequencing data is sufficient to go through with the analysis and figure out where problems might come up at later stages. With Illumina data, this will be done before and after preprocessing, with FastQC which is a quick and simple software for quality control and will provide an overview of the data. PacBio reads will later be assembled with Canu which also includes a quality control step.

#### Preprocessing:
The data I will be working with is already trimmed but might need some more preprocessing. In this case this will be done with Trimmomatic for Illumina sequences which can be used to trim sequencing data and remove adapters. PacBio sequences will be trimmed with Canu during assembly. 

#### Genome Assembly:
In a first step PacBio sequencing data will be assembled using Canu which was developed specifically for the assembly of long reads. The paper reports that this assembly lead to a gap in the assembled genome that could not be closed which is why Illumina and MinION reads will be assembled with Spades in a second step. Spades can deal with various different types of reads at the same time which is why it is a good fit here. Genome assembly is a very time-consuming step. Both steps can take several hours which is why careful planning and execution of these will be crucial in this project. 

#### Assembly Evaluation:
Next, the quality of the assembly needs to be evaluated. Quast will be used to determine the general quality. With MUMmerplot, the two assemblies can be compared by visualizing their whole-genome alignment. 

#### Annotation:
Prokka will be used for annotation. It is a software tool for the fast structural and functional annotation of prokaryotic genomes.

#### Synteny Comparison:
ACT will be used for visualiuation and synteny comparison with a closely related genome. 

#### Mapping:
BWA will be used to align paired-end RNA-seq as well as single-end Tn-seq reads to the assembled genome. 

#### SNP calling:
SNP calling identifies variable sites in the genome and can be carried out with BCFtools which finds SNPs and small insertions and deletions. 

#### Differential Expression Analysis:
Counting the mapped RNA reads for differential expression analysis will be done with Htseq, a Python package developed for exactly this purpose. This is a step that can take up to seven hours for the paired-end reads which will have to be considered. DESeq2 will be used for statistical analysis of the output.

### Timeframe
The official checkpoints for the project are as follows: genome assembly and annotation should be finished by April 17, comparative genomic analyses by April 29 and RNA mapping by May 8 with the final presentations on May 28 and 29. After finishing the plan for the project, work will begin in week 14 which leaves three weeks for assembly and annotation. Since assembly is a step that takes very long and I will be performing it with two different methods enough time should be assigned to complete it. I would like to be done with both methods including quality control and preprocessing at the end of week 15 to have a buffer and enough time to do annotation in week 16. In the following week I will be working on synteny comparison and SNP calling which leaves week 18 and 19 to complete mapping and start differential expression analysis. All analyses should be finished around the middle of week 20 to have enough time to organize and interpret all the results. 

### Organization
Organizing the working directory properly will be important to keep track of the data and all results. I have decided to follow the proposed structure and adjust it if necessary. This means there will be separate directories for analyses, code and data. In analyses, all the steps and methods will have their own folders where files will be numerically organized. Code will also be numerically sorted by method. The data directory will have folders for metadata, raw data and trimmed data. 

### References
Zhang, X., de Maat, V., Guzmán Prieto, A.M., Prajsnar, T.K., Bayjanov, J.R., de Been, M., Rogers, M.R.C., Bonten, M.J.M., Mesnage, S., Willems, R.J.L. & van Schaik, W. 2017, "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum", BMC Genomics, vol. 18, no. 1, pp. 893-12.

## Analyses

### Goals and Hypotheses

As mentioned previously in the project plan, the goal of this project is to reproduce the analyses and results of the paper "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant *Enterococcus faecium* during growth in human serum” (2017) by Xinglin Zhang et al.. This means preprocessing of the available reads which will then be used to generate a de novo assembly. This assembly will be evaluated and should be of high quality in order to be able proceed with the following analyses successfully. Since the authors managed to produce a good assembly and available softwares are very good and efficient, I expect there to be no problems with this step. After evaluation of the assembly, structural and functional annotation will be performed and synteny compared with a closely related genome. In order to do the later, BLAST will be used for a homology search to find a suitable species. Since the authors of the paper mention *Enterococcus faecalis* as the cause of infections similar to *E. faecium*, it might be an interesting option to look at if it comes up as a relevant hit in the results of the search. Then, RNA reads of both BH and human serum growth conditions will be preprocessed if necessary and aligned to the assembly. RNA reads aligned to the different features will be counted and finally differential expression analysis will be performed. If all the steps are successful this should result in similar findings as in the paper which means the discovery of a strongly differentially expressed gene cluster involved in the biosynthesis of purine. 

### 1 Reads Quality Control

Quality control of Illumina raw sequence data was carried out with FastQC. This program provides a quick and simple quality check and gives an overview about basic statistics, per base sequence quality, sequence content, GC-content and N-content, per sequence GC-content, sequence length distribution, sequence duplication levels, overrepresented sequences and Kmer-content. 

For both genomic Illumina read files, the quality checks came back very good, there do not seem to be any problems with low quality. According to the results of this quality control with FastQC, preprocessing of the data is not necessary. Trimming will be performed for learning purposes anyway but the untrimmed reads will be used in further analyses. There are no problems expected that correlate with the quality of the reads. 

FastQC did raise several warnings for all RNA sequencing data files on the other hand. The categories sequence duplication levels and overrepresented sequences produced failures in almost all the files as well as per base sequence content and per base GC content. Sequence length distributions and per tile sequence quality seemed to cause problems multiple times too. These warnings can mostly be explained simply by the nature of RNA sequencing data. It is, for example, expected that some transcripts occur more often than others which can lead to warnings for duplication levels and overrepresented sequences meaning that this does not actually pose a problem to the analysis. The length of RNA transcripts varies as well, so it is entirely normal that the sequence length distribution varies. Per base GC content and per base sequence content are prone to be biased by the overrepresented sequences in our data, so this is nothing to worry about either. Lastly, the FastQC manual suggests that warnings for per tile sequence quality can be ignored if they only seem to affect a very small number of tiles which is the case here. So, all in all the quality of the RNA sequencing data seems to be fine and should not cause problems further down the line which means that no further preprocessing is to be conducted. 

#### Other Questions:

 - What is the structure of a FASTQ file?

 A FASTQ file always starts with the sequence identifier which is preceded by the at sign (@). In the next line follows the actual sequence after which a plus sign (+) is inserted as a separator. The fourth line then includes the base call quality scores. 
 
 - How is the quality of the data stored in the FASTQ files? How are paired reads identified?

The quality scores are encoded with ASCII. This means that the numeric quality scores are stored in the FASTQ file as letters. For paired reads, there should normally be two files available (R1 and R2). 

### 2 Reads Preprocessing

Reads preprocessing of the genomic Illumina reads was performed with Trimmomatic. Since the reads have been preprocessed already, there was not much left to correct, but because the per base sequence quality still dropped a little bit towards the end of the reads, trailing was used to remove low quality bases at the ends. 

After trimming the quality was evaluated with FastQC one more time. 118589 sequences were discarded in each of the two files. This should not cause any problems since the Illumina reads would be used to assemble the genome together with Nanopore reads which means there would most likely be more than enough data available to generate the assembly. In general, the quality has improved in the way that the per base sequence quality score towards the end of the reads is now a lot better. FastQC does raise a warning in sequence length distribution which was expected since the reads will low quality at the ends were shortened and now the reads do not have the same length anymore. 

#### Other Questions:

 - What do the LEADING, TRAILING and SLIDINGWINDOW options do?

Leading removes low quality bases at the beginning and trailing at the end of the read. Slidingwindow trims the read once the average quality in the sliding window falls below a threshold. 


### 3 Genome Assembly

The first step in the assembly of the *E. faecium* genome was an assembly of PacBio sequence data with Canu. Canu is specifically designed for the assembly of long PacBio or Nanopore reads and includes a correction and trimming step. This means separate preprocessing of the PacBio data is not necessary. The genome size of *E. faecium* which needs to be included in the command for this step is derived from appendix 1 and amounts to 3.2 Mbp. Another assembly was later generated with Spades using Illumina and Nanopore data (see Extra Analyses).

The assembly that Canu produced from the PacBio reads has 12 contigs which is better than expected from the results of the paper. The authors have used Celera to assemble the genome which is related to Canu and generated an assembly with 15 contigs, a few of which were discarded due to low coverage.

#### Other Questions:

 - What information can you get from the plots and reports given by the assembler (if you get any)?

Probably the most informative output Canu provides is the report as it shows histogram of read lengths, the histogram of k-mers in the raw and corrected reads, the summary of corrected data, summary of overlaps, the summary of contig lengths and assembly statistics.

 - What intermediate steps generate informative output about the assembly?

Canu's correction and trimming steps generate output that might be interesting to look at and know exactly what happened to the reads. Looking at the graph and ambiguities in it could be informative, too. 

 - What is the difference between a ‘contig’ and a ‘unitig’?

A contig is a consensus sequence made up of a continuous stretch of reads that overlap without any gaps while 
a unitig is a short assembly of DNA sequence with the limitation that the reads cannot contain contradictory overlaps. So in a way, unitigs are considered to be higher quality contigs. 

 - What is the difference between a ‘contig’ and a ‘scaffold’?

A scaffold consists of several contigs separated by gaps that are put together using additional information about the relative position and orientation of the contigs in the genome. 

 - What are the k-mers? What k-mer(s) should you use? What are the problems and benefits of choosing a small kmer? And a big k-mer?

A k-mer is simply a sequence substring of defined length. For short-read assembly, k-mers can be used to correct errors in the assembly by counting the number of times each k-mer appears across the data set and replacing rare k-mers with common ones. More importantly, they are needed in order to build the De Bruijn graph. Choosing the best k-mer can be difficult. The size depends on the read length, sequencing depth and also on the genome size. It is generally recommended to go for k-mer sizes that are not larger than 80 % of the read length and the best results are generated when the algorithm tries multiple different sizes. The Velvetadvisor, for example, would estimate the k-mer size as follows:

k-mer size=e 1 + read length - (k-mer coverage * read length)/(genome coverage)
where, genome coverage = (total number of reads * read length)/ (Estimated genome size)

The k-mer size greatly affects the structure of the De Bruijn graph. A small k-mer size can lead to a very branched De Bruijn graph because repeats and duplications due to errors in reads that are longer than the k-mer size are not recognized and k-mers are wrongly connected which can make assembly difficult and lead to smaller contigs. A large k-mer size can recognize small repeats and therefore reduce the number of branches in the graph. On the other hand, longer k-mers can lead to a very sparse and fragmented graph which can hinder optimal assembly as well. 

 - Some assemblers can include a read-correction step before doing the assembly. What is this step doing?

Canu uses a hybrid error correction step prior to assembly. First, high-identity short-read sequences are mapped to all the long-read sequences. Then, repeats are resolved by placing each short-read sequence in its highest identity repeat copy. Lastly, chimera and trimming problems are resolved within the long-read sequences and a consensus sequence is generated for each long-read sequence based on a multiple alignment of the short-read sequences. 

 - How different do different assemblers perform for the same data?

Even though the assembly of long reads like these is often considered easier as short read assembly, there are problems with high error rates and differences in performance of assemblers that should be taken into account. When looking at performance, Canu is known to be rather slow compared to other tools such as Flye. Different assemblers have different strengths and weaknesses, some are more reliable than others, some have problems detecting plasmids or are known to produce incomplete assemblies. 

 - Can you see any other letter apart from AGTC in your assembly? If so, what are those?

There are no other letters present in the assembly.

### 4 Assembly Evaluation

Assembly evaluation was carried out with Quast and resulted in the following statistics: the number of contigs is 12 and the longest contig has a length of 2774867 bp. The total length of the assembly is 3143732 bp. N50 is 2774867 bp, N75 2774867 bp, L50 1 and
L75 1 as well. The number of N's per 100 kbp is 0 and GC content is 37.82 %. All in all, this means that the assembly is pretty good. The number of contigs is a little high, but since the largest contig almost makes up 90 % of the assembly, this should not be a problem. Shorter contigs could also be discarded here. The reference assembly that was used for comparison has 6 contigs with the largest having a length of 3130373 bp and the length of the assembly being 3259287 bp. N50 and N75 are 3130373 and L50 and L75 are both 1. The number of N's per 100 kbp is 0 here as well and the GC content is 37.69 %. There are definitely differences between the two assemblies, but these are not huge and therefore not concerning. They might be caused by different raw data, for example with a higher coverage or other sequencing technologies, more or less extensive preprocessing, different assembly parameters and polishing steps. Also, several sequences were unassembled by Canu. In this case, the assembly generated during this course project is not as good as the reference assembly. The number of contigs in the reference is only half as big and the longest contig of the reference assembly is significantly longer as well which means that N50 and N75 measures give better results too.

The alignment with MUMer (nucmer algorithm) and visualization using MUMmerplot shows two very long stretches of forward matches in the biggest contigs of both the Canu assembly and the reference assembly which is a good sign of quality. These are however interrupted by a big insertion in the Canu assembly and there are also disagreements and suggested inversion. All in all, this comparison does look good though since the biggest most important contigs match well. 

#### Other Questions

 - What do measures like N50, N90, etc. mean? How can they help you evaluate the quality of your assembly? Which measure is the best to summarize the quality of the assembly (N50, number of ORFs, completeness, total size, longest contig ...)

N50 - the contig length such that using equal or longer contigs sum up to 50% of the bases in the assembly

N75 - the contig length such that using equal or longer contigs sum up to 75% of the bases in the assembly

NG50 - similar to N50, but in relation to the genome length as opposed to the assembly length

L50 - the smallest number of contigs that cover 50 % of the assembly when their lengths are summed up

L75 - the smallest number of contigs that cover 50 % of the assembly when their lengths are summed up

These measures help evaluate the quality of the assembly by giving an overview and summary statistics and providing a basis for comparisons with other assemblies of related species or a reference. They are also useful to look at in order to evaluate if the assembly outcome is as expected. There are probably measures that make more sense to look at than others, for example NG50 is often preferred over N50 because it takes the actual genome size into account, but in general it is the best idea to look at all of them to get a broad picture of the assembly. For this specific assembly, I would suggest that the longest contig is the most informative and able to summarize the assembly well.  

### 5 Annotation

Annotation was performed on the assembly of PacBio reads with the Canu assembler. Prokka, a tool for both structural and functional annotation of prokaryotic genomes, was used to perform this step which resulted in the annotation of 3193 protein coding sequences as well as 140 signal peptides, 101 tRNA, 1 tmRNA and 23 rRNA sequences. The authors of the paper have detected 3217 transcription units. These numbers differ because the authors have counted transcriptional units which include protein-coding genes as well as regulatory non-coding RNAs like microRNAs by mapping RNA transcripts to the assembly. 

#### Other Questions:

 - What types of features are detected by the software? Which ones are more reliable a priori?

Prokka predicts protein-coding sequences, rRNA, tRNA, signal leader peptides and non-coding RNA. The prediction of rRNA genes is very reliable because they are highly conserved and relatively easy to detect and tRNA genes can be predicted with high confidence because of the structural information that is available for them. 

 - Why is it more difficult to do the functional annotation in eukaryotic genomes?

First of all, structural annotation of eukaryotic genomes is a lot more difficult because coding exons are short and separated by introns and the upstream region of the genes is more variable than in prokaryotes. Eukaryotic genomes are very repeat-rich which can lead to a huge amount of BLAST alignment when they are not masked properly and transposon ORFs can look like host genes for an annotation tool. Then, there is also alternative splicing in eukaryotes which means that even if tools are able to predict genes well, not all proteins are predicted which makes functional annotation very difficult. 

 - How many genes are annotated as ‘hypothetical protein’? Why is that so? How would you tackle that problem?

The number of genes annotated as hypothetical proteins is 1389. Annotation of hypothetical proteins happens when there are known homologs available but their function is still unknown (conserved hypothetical) or a gene was predicted but there are no known homologs (hypothetical). Some methods to deal with them and actually find out what their function might be could be comparative genomics, methods looking at protein-protein interactions, clustering approaches where similar genes that are grouped together are assumed to have the same function and genome context methods based on the analysis of fusion events, the conservation of gene neighborhood and the significant co-occurrence of genes across different species.

 - How can you evaluate the quality of the obtained functional annotation?

There are a few tools available that one can use to evaluate the quality of an annotation, but in general one can look at the annotations to see if they are as expected for the species, nothing important is missing and no strange genes are added that do not make sense. Furthermore, it is advisable to compare the annotation to a curated annotation of a reference genome and check if there are major differences. 

 - How comparable are the results obtained from two different structural annotation softwares?

In theory the results should be comparable but since gene prediction can be a difficult task and softwares are not perfect there can be differences in the results. Softwares use different approaches and are more or less conservative with their predictions, so it can make sense to use different tools for structural annotation and compare the results.


### 6 Homology Search and Synteny Comparison

In order to compare synteny with a closely related genome a homology search was carried out using BLASTN to search for homologous species. After that *Enterococcus faecalis* was chosen because it was the best hit with a high e-value and sequence identity. It was also mentioned by the authors of the paper as a strain that causes similar nosocomial infections to *E. faecium*, so a comparison seemed interesting. An alignment of an *E. faecalis* reference assembly with the E. faecium assembly generated with Canu was created by again using BLASTN. This alignment was used in ACT together with the .gff annotation file of the *E. faecium* and the Genbank file of *E. faecalis* to compare synteny. 

Synteny between two genomes exists when there are regions of genome sequence or a certain number of genes that are co-arranged together with the other genome. When comparing *E. faecalis* and *E. faecium* one can observe a great number of these blocks of synteny, some in the same orientation, some inverted. This is to be expected since they are so closely related. These blocks are rather short considering their relatedness, which could indicate that the quality of the Canu assembly is not sufficient. The analysis of differential expression later revealed a gene cluster that is differentially expressed for growth in BH and serum (see 9 Differential Expression). These genes play a role in purine biosynthesis and are suspected to be important for *E. faecium*'s ability to grow in human serum. Synteny of these upregulated genes was established between *E. faecium* and *E. faecalis* which indicates that they have a functional relationship and their order is conserved. 

#### Other Questions:

 - How relevant is the output format that you choose?

When using blastn it can be quite important to specify the output format since this will determine the results one gets. In this specific case, the output also needs to be readable by ACT.

 - How do the resulting hits vary when you change the minimum e-value?

Changing the minimum e-value is a way of filtering the quality of the hits that will be returned as results.

 - How is the alignment score calculated?

The BLASTN alignment score is calculated by assigning score for each aligned pair and then summing up the scores for the whole alignment. For identical letters the score is +2 and for nonidentical it is -3. Gap openings are penalized stronger than gap extensions.  

 - How important is the number of threads when you blast against a database, or against a particular sequence?

It is a lot more important to choose a higher number of threads when blasting against a database since this can take a lot more time and blasting against only one sequence is rather quick.

### 7 Mapping

Read mapping of the Illumina paired end RNA reads was performed with BWA and SAMtools. First, the assembly was indexed, then the reads were mapped to it with BWA and sorted with SAMtools. The results of the mapping step were used to count the reads (see 8 Read Counting) and do differential expression analysis (see 9 Differential Expression). Most reads, namely 98.6, map back to the contigs. Reasons why not all the reads map back could be that the assembly is incomplete or flawed. 70.6 % of the reads map to genes, so 16571329 reads do not map to features. When counting the reads with HTSeq (again, see 8 Read Counting), a filter was applied with the -t tag to only count reads that map to protein-coding sequences which would explain the large number of reads that did not map to features. These could for example consist of tRNA, mtRNA or rRNA. IGV was used to look at read coverage. To do this, index files were created with IGVtools that were then loaded together with the bam files. Read coverage across the genome is medium to low in most regions with a few regions that show very high coverage. Regions with very low to no coverage indicate non-coding regions or very lowly expressed genes, wheras highly expressed genes should have very high coverage. Coverage patterns are very similar across all replicates, also between the different growth conditions. As expected, the region where the differentially expressed pur-genes (see 9 Differential Expression) are located, showed strong differences in coverage between the different conditions. 

#### Other Questions: 

 - What potential issues can cause mRNA reads not to map properly to genes in the chromosome? Do you expect this to differ between prokaryotic and eukaryotic projects?

In Eukaryotes, RNA splicing is performed in order to remove intronic regions. This has to be considered when aligning mRNA to the genome and makes mapping more difficult. Another problem to deal with in eukaryotic projects is alternative splicing where some exons are not included in the final mRNA product.

### 8 Read Counting

Mapped reads were counted with the Python package HTSeq that uses the BAM file generated in the previous step as well as the gff file produced during annotation. The histogram visualizing the distribution of counts per gene shows that the read count for most genes is actually zero which means that most genes are not expressed. In theory, a gene is expressed if the read count is > 0, so if at least one transcript is present. In order to do proper statistical analysis, this does not provide a lot of information though and it is more informative to investigate if expression is significantely different between the two conditions. In order to account for the fact that the counts for most genes are zero, the counts are normalized. The histogram of the normalized counts shows something like a normal distribution. 


### 9 Differential Expression

Differential expression analysis of the RNA sequence data retrieved from cultures grown in BH and serum was performed with the R package DESeq2. After creating the DESeqDataSet object from the read counts the DESeq2 pipelinewas executed. Bevor inspecting the result table further, insignificant results were filtered out. Usually, results are understood as significant if the p-value is below a predefined threshold (often 0.05 or 0.01). In the case of differential expression analysis, there are thousands of genes involved and accepting a p-value of for example 0.05 as significant would mean a great risk of false positives. To avoid and minimize that risk, DESeq2 calculates an adjusted p-value that is corrected with the Benjamini–Hochberg procedure. This value is then used to filter significant results. In this case, removing all the results with an adjusted p-value > 0.1 resulted in a total of 664 significantely differentially expressed genes. Analyzing these results further and looking at the genes with the strongest up- or downregulation revealed a differentially expressed gene cluster involved in purine biosynthesis just as the results of the paper show as well.

In general, the different samples and different conditions do cluster together as shown in the heat map. The PCA performed on the rlog-transformed data indicates that 99 % of the variance can be explained by the the different growth mediums.

#### Other Questions:

 - How do the different samples and replicates cluster together?

The PCA performed on the rlog-transformed data shows that the samples cluster together according to the different growth conditions that were applied and 99 % of the variance is actually explained by these different conditions.

 - What is the q-value and how does it differ from the p-value? Which one should you use to determine if the result is statistically significant?

 The q-value is a measure of significance analogous to the p-value. They do, however, differ in the sense that the p-value is a measure of significance in terms of the false positive rate whereas the q-value is in terms of the false discovery rate. The difference here is that the false positive rate is rate that features that are really not significant are called significant while the false discovery rate is the rate that features that are called significant are actually not. So, the q-value is a measure that focuses on the features that are called significant, it takes the multiple testing problem into account and limits false positives but tries to minimize the number of false negatives at the same time which makes it a more suitable measure than the p-value for differential expression analysis. 

 - Do you need a normalization step? What would you normalize against? Does DESeq do it?

In general, it is necessary to perform a normalization step in order to make accurate comparisons of gene expression. However, DESeq2 needs raw counts as input data and it is therefore not recommended to do any sort of normalization beforehand. DESeq2 does its own internal normalization steps accounting for sequencing depth and RNA composition with the median of ratios method.

 - What would you do to increase the statistical power of your expression analysis?

Independent filtering is a way of increasing statistical power by removing weakly expressed genes from the analysis which then improves the multiple testing adjustments that are made. This means, more genes will be found to be significant. 

## Extra Analyses

### 10 Genome Assembly with Illumina and Nanopore Reads

Illumina and Nanopore reads were then assembled together using Spades which specializes in assembling data from different sequencing methods as well as long and short reads in one step. Illumina and Nanopore files were included in the command according to the manual and to reduce running time the Kmer size was set to 55 instead of trying various different ones. As stated above, it is recommended to use a k-mer size that is smaller than 80 % of the read lenght wich is 100 bp in the case of the short Illumina reads here. 

This second assembly using Spades turned out less good than the assembly generated with Canu. There are 52 contigs in this one which is quite a lot compared to the 12 that Canu produced. However, very short contigs should probably be discarded and statistics are calculated based on contigs that are larger than 500 bp only which makes the total number of contigs 18 which is still higher. The total length of the Spades assembly amounts to 3135284 bp which is shorter than the Canu assembly and therefore further away from the recerence assembly as well. The longest contig is 986763 bp, N50 is 840748, N75 266941, L50 2 and L75 is 4. The number of N's per 100 kbp is zero and the GC content is 37.67 %. The MUMmerplot comparison shows that there are long stretches of forward matches and especially the larger contigs match the Canu assembly very well, but there are still quite a lot of rather big insertions in the Canu assembly. The MUMmerplot comparison with the reference genome, however, is very similar to the comparison plot of the Canu assembly: There are two very long stretches of forward matches in the largest contigs with one big insertion in the spades assembly and another smaller one. The assembly could probably be improved by trying out various different k-mer sizes, but the overall quality is fine and since this assembly will not be used in any other analysis no further runs with different parameters will be carried out. 

### 11 Plasmid Identification

First, plasmid identification was attempted using the Spades plasmid algorithm on Illumina and Nanopore reads which did not produce the expected results. Instead of identifying six plamsmids as the authods of the paper did, it did not recognize any. Even after switching and trying several smaller k-mer sizes which is supposed to facilitate plasmid identification as well as removing the Nanopore reads that can cause problems with the plasmid algorithm, the results did not change. Then, PlasmidFinder was used on the assembly of PacBio reads which finally identified six plasmids. 

### 12 SNP calling

SNP calling was performed with SAMtools and BCFtools with the reads generated by the authors of this paper as well as another read file from the study "Real time application of whole genome sequencing for outbreak investigation – what is an achievable turnaround time?" by McGann et al. in 2015 and 2016. This was retrieved from SRA using fastq-dump. Analogous to the read mapping performed earlier, these reads were aligned to the assembly with BWA and sorted with SAMtools. Then the SNPs of all the reads were called with SAMtools and BCFtools and visualized with IGV. This resulted in a total of 9351 variants; 8645 of these are SNPs and 706 indels. 8079 of the variants have a quality score higher than 10, 5931 higher than 50 and 5549 higher than 100. Distribution of the variants can be observed in [this](https://github.com/rpensch/genome_analysis/blob/master/results/12_snp_calling/snps_overview.png) plot. There is one region towards the begnning of the assembly where a high number of variants accumulate as well as four to five other regions with notably elevated variant levels. 

#### Other Questions:

 - What is the structure of a SAM file, and how does it relate to a BAM file?

A SAM file generally stores read alignments against a reference sequence and consists of a header and the alignment. Lines in the header start with the at sign (@). The alignment section consists of eleven mandatory fields and other optional ones. The mandatory fields are query name of the read or read pair, bitwise flag, reference name, 1-based leftmost position of clipped alignment, mapping quality, extended CIGAR STRING (Concise Idiosyncratic Gapped Alignment Report), mate reference name, 1-based leftmost mate position, inferred insert size, query sequence on the same strand as the reference and query quality. A BAM file stores exactely the same information but in compressed binary representation. 

 - What is the structure of vcf and bcf files?

A vcf file contains lines of meta data starting with '##', a header starting with '#' and data lines with eight mandatory and several optional columns that store information about a position in the genome and genotype information on samples for each position. A bcf file is the binary, compressed equivalent of a vcf file.

 - What is the difference between the variant quality, the mapping quality and the fastq quality?

Fastq quality indicates how well supported a sequenced base is, mapping quality measures how confident we can be that a read has been aligned with the reference in the right place and variant quality estimates a likelihood for ovserving a call by chance.

### 13 Antibiotic Resistance Potential

Resistance to several different antibiotics is a known problem for *E. faecium* and a great concern. Analogous to the approach the authors of the paper took, ResFinder was used to discover genes involved in acquired antibiotic resistance. It found the vanHAX operon which is responsible for vancomycin resistance and macrolide resistance gene mrs(C). The dfrG gene which causes trimethoprim resistance that the authors found was not discovered here which could be due to quality differences of the assembly. ResFinder did, however, get a hit for an aminoglycoside restistance gene that was not mentioned in the paper. This could again be due to some error in the assembly or the gene might just have been added to the antibiotic restistance gene database that the tool relies on inbetween the analysis of the authors and this one. In conclusion, resistance to vancomycin and a macrolide resistance seems very probable wheras the aac(6')-Ii gene involved in aminoglycoside resistance need further investigation. 

### References

 - Deonier, R. C., Tavaré, S., & Waterman, M. S. (2005). Computational genome analysis: An introduction. New York: Springer. doi:10.1007/0-387-28807-4 
 - Cha, S., & Bird, D. M. (2016). Optimizing k-mer size using a variant grid search to enhance de novo genome assembly. Bioinformation, 12(2), 36–40. doi:10.6026/97320630012036
 - Mahadik, K., Wright, C., Kulkarni, M., Bagchi, S., & Chaterji, S. (2019). Scalable genome assembly through parallel de bruijn graph construction for multiple k-mers. London: Nature Publishing Group. doi:10.1038/s41598-019-51284-9
 - Koren, S., Walenz, B. P., Berlin, K., Miller, J. R., Bergman, N. H., & Phillippy, A. M. (2017). Canu: Scalable and accurate long-read assembly via adaptive k -mer weighting and repeat separation. United States: Cold Spring Harbor Laboratory Press. doi:10.1101/gr.215087.116
 - Koren, S., Schatz, M. C., Walenz, B. P., Martin, J., Howard, J. T., Ganapathy, G., . . . Joint Genome Institute (JGI). (2012). Hybrid error correction and de novo assembly of single-molecule sequencing reads. United States: Nature Publishing Group. doi:10.1038/nbt.2280
 - Wick, R. R., & Holt, K. E. (2019). Benchmarking of long-read assemblers for prokaryote whole genome sequencing [version 1; peer review: 4 approved]. F1000research, 8, 2138. doi:10.12688/f1000research.21782.1
 - Sahu, A., Li, N., Dunkel, I., & Chung, H. (2020). EPIGENE: Genome-wide transcription unit annotation using a multivariate probabilistic model of histone modifications. England: BMC. doi:10.1186/s13072-020-00341-z
 - Seemann, T. (2014). Prokka: Rapid prokaryotic genome annotation. England: doi:10.1093/bioinformatics/btu153
 - Sivashankari, S., & Shanmughavel, P. (2006). Functional annotation of hypothetical proteins - A review. Singapore: doi:10.6026/97320630001335
 - Yandell, M., & Ence, D. (2012). A beginner's guide to eukaryotic genome annotation. Nature Reviews. Genetics, 13(5), 329-342. doi:10.1038/nrg3174
 - Liu, D., Hunt, M., & Tsai, I. J. (2018). Inferring synteny between genome assemblies: A systematic evaluation. BMC Bioinformatics, 19(1), 26-13. doi:10.1186/s12859-018-2026-4
 - McGann, P., Bunin, J. L., Snesrud, E., Singh, S., Maybank, R., Ong, A. C., . . . Lesho, E. (2015;2016;). Real time application of whole genome sequencing for outbreak investigation – what is an achievable turnaround time? Diagnostic Microbiology & Infectious Disease, 85(3), 277-282. doi:10.1016/j.diagmicrobio.2016.04.020 
 - Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., . . . 1000 Genome Project Data Processing Subgroup. (2009). The sequence alignment map format and SAMtools. England: Oxford University Press. doi:10.1093/bioinformatics/btp352
 - Storey, J. D. (2003). The positive false discovery rate: A bayesian interpretation and the q-value. The Annals of Statistics, 31(6), 2013-2035. doi:10.1214/aos/1074290335
 - Storey, J. D., & Tibshirani, R. (2003). Statistical significance for genomewide studies. Proceedings of the National Academy of Sciences of the United States of America, 100(16), 9440-9445. doi:10.1073/pnas.1530509100