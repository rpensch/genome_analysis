# Genome Analysis
## Project Plan
### Introduction
In the course of this project I am aiming to reproduce the results of the following research paper: 

"RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum” (2017) by Xinglin Zhang et al. 

*Enterococcus faecium* is a bacterium that normally occupies the human digestive system but can often cause infections of the bloodstream in hospitalized patients. This is dangerous since some strains have acquired resistance against the most common antibiotics. Surviving and thriving in such suboptimal and nutrient-poor conditions as they are in the human blood is a challenge and the mechanisms that make it possible for *E. faecium* are to be analyzed in this study. In a first step, the question whether and how gene expression differs in human serum and conditions of a conventional nutrient-rich growth medium is of interest. This is analyzed by examining the transcriptomes of *E. faecium* grown in both environments. Later, genes that promote growth in human serum are identified by comparing transposon sequencing (Tn-seq) data of transposon mutant libraries of the bacterium from both serum and plain growth medium. 

The data that is available to try and answer these questions is for one sequencing data from Illumina HiSeq 100 bp paired-end sequencing, Pacific Biosciences RS II SMRT technology long-read sequencing and Oxford Nanopore Technologies MinION with R7 flowcell chemistry to assemble the genome of *Enterococcus faecium*. Transcriptome data is available from 100 bp paired-end sequencing with the Illumina HiSeq 2500 platform and Tn-seq data from 50 nt single-end sequencing with Illumina HiSeq 2500. 

### Methods
The steps and analyses required for carrying out this project successfully will be described below in roughly the order they will be performed.

#### Quality Control:
Before starting the project it is important to make sure the quality of the sequencing data is sufficient to go through with the analysis and figure out where problems might come up at later stages. With Illumina data, this will be done before and after preprocessing, with FastQC which is a quick and simple software for quality control and will provide an overview of the data. PacBio reads will later be assembled with Canu which also includes a quality control step.

#### Preprocessing:
The data I will be working with is already trimmed but might need some more preprocessing. In this case this will be done with Trimmomatic for Illumina sequences which can be used to trim sequencing data and remove adapters. PacBio sequences will be trimmed with Canu during assembly. 

#### Genome Assembly:
In a first step PacBio sequencing data will be assembled using Canu which was developed specifically for the assembly of long reads. The paper reports that this assembly lead to a gap in the assembled genome that could not be closed which is why Illumina and MinION reads will be assembled with Spades in a second step. Spades can deal with various different types of reads at the same time which is why it is a good fit here. Genome assembly is a very time-consuming step. Both steps can take several hours which is why careful planning and execution of these will be crucial in this project. 

#### Assembly Evaluation:
Next, the quality of the assembly needs to be evaluated. Quast will be used to determine the general quality. With MUMmerplot, the two assemblies can be compared by visualizing their whole-genome alignment. 

#### SNP calling:
SNP calling identifies variable sites in the genome and can be carried out with BCFtools which finds SNPs and small insertions and deletions. 

#### Annotation:
Prokka will be used for annotation. It is a software tool for the fast structural and functional annotation of prokaryotic genomes.

#### Synteny Comparison:
ACT will be used for visualiuation and synteny comparison with a closely related genome. 

#### Mapping:
BWA will be used to align paired-end RNA-seq as well as single-end Tn-seq reads to the assembled genome. 

#### Differential Expression Analysis:
Counting the mapped RNA reads for differential expression analysis will be done with Htseq, a Python package developed for exactly this purpose. This is a step that can take up to seven hours for the paired-end reads which will have to be considered. DESeq2 will be used for statistical analysis of the output.

### Timeframe
The official checkpoints for the project are as follows: genome assembly and annotation should be finished by April 17, comparative genomic analyses by April 29 and RNA mapping by May 8 with the final presentations on May 28 and 29. After finishing the plan for the project, work will begin in week 14 which leaves three weeks for assembly and annotation. Since assembly is a step that takes very long and I will be performing it with two different methods enough time should be assigned to complete it. I would like to be done with both methods including quality control and preprocessing at the end of week 15 to have a buffer and enough time to do annotation in week 16. In the following week I will be working on synteny comparison and SNP calling which leaves week 18 and 19 to complete mapping and start differential expression analysis. All analyses should be finished around the middle of week 20 to have enough time to organize and interpret all the results. 

### Organization
Organizing the working directory properly will be important to keep track of the data and all results. I have decided to follow the proposed structure and adjust it if necessary. This means there will be separate directories for analyses, code and data. In analyses, all the steps and methods will have their own folders where files will be numerically organized. Code will also be numerically sorted by method. The data directory will have folders for metadata, raw data and trimmed data. 

### References
Zhang, X., de Maat, V., Guzmán Prieto, A.M., Prajsnar, T.K., Bayjanov, J.R., de Been, M., Rogers, M.R.C., Bonten, M.J.M., Mesnage, S., Willems, R.J.L. & van Schaik, W. 2017, "RNA-seq and Tn-seq reveal fitness determinants of vancomycin-resistant Enterococcus faecium during growth in human serum", BMC Genomics, vol. 18, no. 1, pp. 893-12.

## Analyses

### Reads Quality Control

Quality control of Illumina raw sequence data was carried out with FastQC. This program provides a quick and simple qulaity check and gives an overview about basic statistics, per base sequence quality, sequence content, GC-content and N-content, per sequence GC-content, sequence length distribution, sequence duplication levels, overrepresented sequences and Kmer-content. 

#### Results

For both genomic Illumina read files, the quality checks came back very good, there do not seem to be any problems with low quality. According to the results of this quality control with FastQC, preprocessing of the data is not necessary and reads will therefore not be trimmed any further. There are no problems expected in further analses that correlate with the quality of the reads. 

FastQC did raise several warnings for all RNA sequencing data files on the other hand. The categories sequence duplication levels and overrepresented sequences produced failures in almost all the files as well as per base sequence content and per base GC content. Sequence length distributions and per tile sequence quality seemed to cause problems multiple times too. These warnings can mostly be explained simply by the nature of RNA sequencing data. It is, for example, expected that some transcripts occur more often than others which can lead to warnings for duplication levels and overrepresented sequences meaning that this does not actually pose a problem to the analysis. The length of RNA transcripts varies as well, so it is entirely normal that the sequence length distribution varies. Per base GC content and per base sequence content are prone to be biased by the overrepresented sequences in our data, so this is nothing to worry about either. Lastly, the FastQC manual suggests that warnings for per tile sequence quality can be ignored if they only seem to affect a very small number of tiles which is the case here. So, all in all the quality of the RNA sequencing data seems to be fine and should not cause problems further down the line which means that no further preprocessing is to be conducted. 

#### Other Questions:

 - What is the structure of a FASTQ file?

 A FASTQ file always starts with the sequence identifier which is preceded by the at sign (@). In the next line follows the actual sequence after which a plus sign (+) is inserted as a seperator. The fourth line then includes the base call quality scores. 
 
 - How is the quality of the data stored in the FASTQ files? How are paired reads identified?

The quality scores are encoded with ASCII. This means that the numeric quality scores are stored in the FASTQ file as letters. For paired reads, there should normally be two files available (R1 and R2). 


### Genome Assembly

The first step in the assembly of the *E. faecium* genome was an assembly of PacBio sequence data with Canu. Canu is specifically designed for the assembly of long PacBio or Nanopore reads and includes a correction and trimming step. This means seperate preprocessing of the PacBio data is not necessary. The genome size of *E. faecium* which needs to be included in the command for this step is derived from appendix 1 and amounts to 3.2 Mbp.

Illumina and Nanopore reads were then assembled together using Spades which specializes in assembling data from different sequencing methods as well as long and short reads in one step. Illumina and Nanopore files were included in the command according to the manual and to reduce running time the Kmer size was set to 55 instead of trying various different ones. 

#### Results

#### Other Questions:

 - What is the difference between a ‘contig’ and a ‘unitig’?

A contig is a consensus sequence made up of a continuous stretch of reads that overlap without any gaps while 
a unitig is a short assembly of DNA sequence with the limitation that the reads cannot contain contradictory overlaps. So in a way, unitigs are considered to be higher quality contigs. 

 - What is the difference between a ‘contig’ and a ‘scaffold’?

A scaffold consists of several contigs separated by gaps that are put together using additional information about the relative position and orientation of the contigs in the genome. 

 - What are the k-mers? What k-mer(s) should you use? What are the problems and benefits of choosing a small kmer? And a big k-mer?

A k-mer is simply a sequence substring of defined length. For short-read assembly, k-mers can be used to correct errors in the assembly by counting the number of times each k-mer appears across the data set and replacing rare k-mers with common ones. More importantly, they are needed in order to build the De Bruijn graph. Choosing the best k-mer can be difficult. The size depends on the read length, sequencing depth and also on the genome size. It is generally recommended to go for k-mer sizes that are not larger than 80 % of the read lenght and the best results are generated when the algorithm tries multiple different sizes. The Velvetadvisor, for example, would estimate the k-mer size as follows:

k-mer siz=e 1 + read length - (k-mer coverage * read length)/(genome coverage)
where, genome coverage = (total number of reads * read length)/ (Estimated genome size)

The k-mer size greatly affects the structure of the De Bruijn graph. A small k-mer size can lead to a very branched De Bruijn graph because repeats and duplications due to errors in reads that are longer than the k-mer size are not recognized and k-mers are wrongly connected which can make assembly difficult and lead to smaller contigs. A large k-mer size can recognize small repeats and therefore reduce the number of branches in the graph. On the other hand, longer k-mers can lead to a very sparse and fragmented graph which can hinder optimal assembly as well. 

 - Some assemblers can include a read-correction step before doing the assembly. What is this step doing?

 Canu uses a hybrid error correction step prior to assembly. First, high-identity short-read sequences are mapped to all the long-read sequences. Then, repeats are resolved by placing each short-read sequence in its highest identity repeat copy. Lastly, chimera and trimming problems are resolved within the long-read sequences and a consensus sequence is generated for each long-read sequence based on a multiple alignment of the short-read sequences. 
 - 



### References

 - Deonier, R. C., Tavaré, S., & Waterman, M. S. (2005). Computational genome analysis: An introduction. New York: Springer. doi:10.1007/0-387-28807-4 
 - Cha, S., & Bird, D. M. (2016). Optimizing k-mer size using a variant grid search to enhance de novo genome assembly. Bioinformation, 12(2), 36–40. https://doi.org/10.6026/97320630012036
 - Mahadik, K., Wright, C., Kulkarni, M., Bagchi, S., & Chaterji, S. (2019). Scalable genome assembly through parallel de bruijn graph construction for multiple k-mers. London: Nature Publishing Group. doi:10.1038/s41598-019-51284-9
 - Koren, S., Walenz, B. P., Berlin, K., Miller, J. R., Bergman, N. H., & Phillippy, A. M. (2017). Canu: Scalable and accurate long-read assembly via adaptive k -mer weighting and repeat separation. United States: Cold Spring Harbor Laboratory Press. doi:10.1101/gr.215087.116
 - Koren, S., Schatz, M. C., Walenz, B. P., Martin, J., Howard, J. T., Ganapathy, G., . . . Joint Genome Institute (JGI). (2012). Hybrid error correction and de novo assembly of single-molecule sequencing reads. United States: Nature Publishing Group. doi:10.1038/nbt.2280