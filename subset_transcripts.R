#!/usr/bin/env Rscript
library(Biostrings)

# Script to subset GENCODE transcripts to 20,000 for the workshop dataset
# This reduces computational load while maintaining biological relevance

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input FASTA file).n", call.=FALSE)
}

fasta_file = args[1]
output_file = "gencode_20k_transcripts.fa"
target_size = 20000

print(paste("Reading transcripts from:", fasta_file))
fasta = readDNAStringSet(fasta_file)
total_transcripts = length(fasta)

print(paste("Total transcripts in file:", total_transcripts))

if (total_transcripts < target_size) {
  stop(paste("Error: FASTA file has only", total_transcripts, 
             "transcripts but you requested", target_size, "transcripts.n"))
}

# Randomly sample 20,000 transcripts
set.seed(42)  # For reproducibility
selected_indices = sample(1:total_transcripts, target_size, replace=FALSE)
fasta_subset = fasta[selected_indices]

print(paste("Subset contains:", length(fasta_subset), "transcripts"))

# Write subset to file
writeXStringSet(fasta_subset, filepath=output_file, format="fasta")

print(paste("Subset transcripts written to:", output_file))
print("Ready to use with create_simulated_dataset.R")
