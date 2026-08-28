#!/usr/bin/env Rscript
library(Biostrings)

# Script to create a simulated RNA-seq dataset for needle workshop
# Ultra-fast version using vectorized operations

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input FASTA file).n", call.=FALSE)
}

fastapath = args[1]

# Read the reference transcriptome
print("Reading transcripts...")
fasta = readDNAStringSet(fastapath)
num_transcripts = length(fasta)
tx_lengths = width(fasta)

# Define cancer-associated genes (5% of transcripts)
num_cancer_genes = max(50, round(0.05 * num_transcripts))
cancer_gene_indices = sample(1:num_transcripts, num_cancer_genes)

# Simulation parameters
num_patients = 24
num_replicates = 2
coverage = 30
readlen = 75

# Risk assignments
risk_assignments = rep(c("low", "medium", "high"), length.out = num_patients)

print(paste("Total transcripts:", num_transcripts))
print(paste("Cancer-associated genes:", num_cancer_genes))
print(paste("Patient samples:", num_patients))

# Create output directories
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("data/cancer_gene_expression")) dir.create("data/cancer_gene_expression")

# Pre-calculate base reads per transcript for all coverages
base_readspertx = round(30 * tx_lengths / readlen)

# Main simulation loop
for (patient_id in 1:num_patients) {
  risk_level = risk_assignments[patient_id]
  
  # Initialize fold changes
  fold_changes = rep(1, num_transcripts)
  
  # Assign expression patterns based on risk level
  if (risk_level == "low") {
    num_elevated = sample(0:1, 1)
    if (num_elevated > 0) {
      elevated_genes = sample(cancer_gene_indices, num_elevated)
      fold_changes[elevated_genes] = sample(c(2, 3), num_elevated, replace = TRUE)
    }
  } else if (risk_level == "medium") {
    num_elevated = sample(2:3, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_changes[elevated_genes] = sample(c(3, 4, 5, 6), num_elevated, replace = TRUE)
  } else if (risk_level == "high") {
    num_elevated = sample(4:6, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_changes[elevated_genes] = sample(c(5, 6, 7, 8, 9, 10), num_elevated, replace = TRUE)
  }
  
  # Write cancer gene expression file
  elevated_indices = which(fold_changes > 1)
  if (length(elevated_indices) > 0) {
    elevated_genes_info = data.frame(
      gene_id = names(fasta[elevated_indices]),
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
  
  # Calculate adjusted reads per transcript for this coverage
  current_coverage = 30 + ((patient_id - 1) %/% 8) * 15
  readspertx = round(base_readspertx * (current_coverage / 30))
  
  # Create patient directory
  patient_dir = paste0("data/patient_", sprintf("%02d", patient_id), "_", risk_level)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Generate reads for each replicate
  for (rep in 1:num_replicates) {
    set.seed(patient_id * 1000 + rep)
    
    # Build all reads for this sample
    read_count = 0
    fastq_lines = character()
    
    for (tx_idx in 1:num_transcripts) {
      num_reads = round(readspertx[tx_idx] * fold_changes[tx_idx])
      
      if (num_reads > 0) {
        sequence = as.character(fasta[[tx_idx]])
        seq_length = nchar(sequence)
        
        # Sample read positions
        start_positions = sample(1:(seq_length - readlen + 1), num_reads, replace = TRUE)
        end_positions = start_positions + readlen - 1
        
        # Extract reads (vectorized)
        for (i in 1:num_reads) {
          read_count = read_count + 1
          read_seq = substr(sequence, start_positions[i], end_positions[i])
          quality = strrep("I", readlen)
          
          fastq_lines = c(fastq_lines,
                         paste0("@read_", read_count, "/1"),
                         read_seq,
                         "+",
                         quality)
        }
      }
    }
    
    # Write FASTQ file
    output_file = paste0(patient_dir, "/sample_", sprintf("%02d", rep), ".fq")
    writeLines(fastq_lines, output_file)
    system(paste("gzip -f", output_file))
    
    print(paste("Patient", sprintf("%02d", patient_id), "-", risk_level, 
                "- Rep", rep, "- Reads:", read_count, "- Coverage:", current_coverage))
  }
}

print("Dataset generation complete!")
print(paste("Output: data/patient_XX_[risk_level]/sample_0X.fq.gz"))
