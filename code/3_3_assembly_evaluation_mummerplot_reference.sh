
module load bioinfo-tools 
module load MUMmer

nucmer --mum -p 1_assembly_reference \  
data/reference/GCF_001720945.1_ASM172094v1_genomic.fna \
analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta 


mummerplot -p 1_assembly_mummerplot_reference --png 1_assembly_reference.delta 
