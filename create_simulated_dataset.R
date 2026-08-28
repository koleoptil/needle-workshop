#!/usr/bin/env Rscript
library(Biostrings)

# Script to create a simulated RNA-seq dataset for needle workshop
# Minimalist version - small dataset for quick generation

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input FASTA file).n", call.=FALSE)
}

fastapath = args[1]

print("Reading transcripts...")
fasta = readDNAStringSet(fastapath)
num_transcripts = length(fasta)

# DRASTICALLY reduce dataset size
subset_size = 500  # Use only 500 transcripts instead of 20000
subset_indices = sample(1:num_transcripts, subset_size)
fasta_subset = fasta[subset_indices]

# Define cancer-associated genes
num_cancer_genes = 25
cancer_gene_indices = sample(1:subset_size, num_cancer_genes)

# Simulation parameters - keep simple
num_patients = 8  # Reduced from 24
num_replicates = 2
readlen = 75

risk_assignments = rep(c("low", "medium", "high"), length.out = num_patients)

print(paste("Using subset of transcripts:", subset_size))
print(paste("Total patients:", num_patients))

# Create output directories
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("data/cancer_gene_expression")) dir.create("data/cancer_gene_expression")

# Main loop
for (patient_id in 1:num_patients) {
  risk_level = risk_assignments[patient_id]
  
  fold_changes = rep(1, subset_size)
  
  if (risk_level == "low") {
    num_elevated = 1
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_changes[elevated_genes] = 2
  } else if (risk_level == "medium") {
    num_elevated = 2
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_changes[elevated_genes] = 4
  } else if (risk_level == "high") {
    num_elevated = 3
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_changes[elevated_genes] = 6
  }
  
  # Write cancer gene info
  elevated_indices = which(fold_changes > 1)
  if (length(elevated_indices) > 0) {
    elevated_genes_info = data.frame(
      gene_id = names(fasta_subset[elevated_indices]),
      fold_change = fold_changes[elevated_indices],
      risk_level = risk_level,
      patient_id = patient_id,
      stringsAsFactors = FALSE
    )
    write.table(elevated_genes_info, 
                file = paste0("data/cancer_gene_expression/patient_", 
                             sprintf("%02d", patient_id), "_genes.tsv"),
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # Create patient directory
  patient_dir = paste0("data/patient_", sprintf("%02d", patient_id), "_", risk_level)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Generate reads for each replicate
  for (rep in 1:num_replicates) {
    set.seed(patient_id * 100 + rep)
    
    fastq_lines = character()
    read_id = 0
    
    # Simple read generation
    for (tx_idx in 1:subset_size) {
      sequence = as.character(fasta_subset[[tx_idx]])
      seq_length = nchar(sequence)
      
      # Generate 10-50 reads per transcript (very small)
      num_reads = round(runif(1, 10, 50) * fold_changes[tx_idx])
      
      for (r in 1:num_reads) {
        read_id = read_id + 1
        pos = sample(1:(seq_length - readlen + 1), 1)
        read_seq = substr(sequence, pos, pos + readlen - 1)
        quality = strrep("I", readlen)
        
        fastq_lines = c(fastq_lines,
                       paste0("@read_", read_id),
                       read_seq,
                       "+",
                       quality)
      }
    }
    
    # Write FASTQ
    output_file = paste0(patient_dir, "/sample_", sprintf("%02d", rep), ".fq")
    writeLines(fastq_lines, output_file)
    system(paste("gzip -f", output_file))
    
    print(paste("Patient", sprintf("%02d", patient_id), "-", risk_level, 
                "- Rep", rep, "- Reads:", read_id))
  }
}

print("Complete!")
print("Output: data/patient_XX_[risk_level]/sample_0X.fq.gz")
