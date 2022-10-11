# kmer_creation
This is an R program to create kmers

# How to run?
Download this R program to any location where you have your fasta file and just execute using:

Rscript kmer_create.R size_of_kmer_length_you_want_to_create

#Dependencies:
R version >= 4.0

# Libraries
1.Biostrings
2.dplyr
3.doParallel
4.data.table
5.seqRFLP
6.stringr

#Input: Fasta Files. 
For example you have fasta files of 4 Fungal proteins:

1. Yeast.fasta
2. Armillaria_ostoyae.fasta
3. Trichoderma_atroviridae.fasta
4. Candida.fasta

#Output: final_"user_supplied_kmerlength"_mer_count.txt
