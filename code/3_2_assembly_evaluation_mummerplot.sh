
module load bioinfo-tools 
module load MUMmer

nucmer --mum -p 1_assembly_nucmer \  
analyses/2_genome_assembly/1_efaecium_canu/efaecium_canu.contigs.fasta \
analyses/2_genome_assembly/2_efaecium_spades/efaecium_spades.contigs.fasta

mummerplot -c -p 5_assembly_mummerplot 1_assembly_mummerplot_SNP_cov -S --png 1_assembly_nucmer.delta
mummerplot -c -p 6_assembly_mummerplot 1_assembly_mummerplot_cov --png 1_assembly_nucmer.delta 
mummerplot -p 7_assembly_mummerplot 1_assembly_mummerplot_SNP -S --png 1_assembly_nucmer.delta
mummerplot -p 8_assembly_mummerplot 1_assembly_mummerplot --png 1_assembly_nucmer.delta 

#Q: what about -R and -Q (student manual)?
