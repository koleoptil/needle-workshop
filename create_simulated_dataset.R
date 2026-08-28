#!/usr/bin/env Rscript
library(Biostrings)
library(edgeR)

# Script to create a simulated RNA-seq dataset for needle workshop
# Optimized for speed - generates pooled reads per sample

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input FASTA file).n", call.=FALSE)
}

fastapath = args[1]

# Read the reference transcriptome
fasta = readDNAStringSet(fastapath)
num_transcripts = length(fasta)

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

print("Starting simulated dataset generation...")
print(paste("Total transcripts:", num_transcripts))
print(paste("Cancer-associated genes:", num_cancer_genes))
print(paste("Patient samples:", num_patients))
print(paste("Replicates per patient:", num_replicates))

# Create output directories
if (!dir.exists("data")) {
  dir.create("data")
}
if (!dir.exists("data/cancer_gene_expression")) {
  dir.create("data/cancer_gene_expression")
}

# Helper function to write FASTQ records
write_fastq_records <- function(sequences, qualities, output_file) {
  handle <- file(output_file, "w")
  
  for (i in 1:length(sequences)) {
    writeLines(c(
      paste0("@read_", i, "/1"),
      sequences[i],
      "+",
      qualities[i]
    ), handle)
  }
  
  close(handle)
  system(paste("gzip -f", output_file))
}

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
      fold_change_values = sample(c(2, 3), num_elevated, replace = TRUE)
      fold_changes[elevated_genes] = fold_change_values
    }
  } else if (risk_level == "medium") {
    num_elevated = sample(2:3, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_change_values = sample(c(3, 4, 5, 6), num_elevated, replace = TRUE)
    fold_changes[elevated_genes] = fold_change_values
  } else if (risk_level == "high") {
    num_elevated = sample(4:6, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_change_values = sample(c(5, 6, 7, 8, 9, 10), num_elevated, replace = TRUE)
    fold_changes[elevated_genes] = fold_change_values
  }
  
  # Write cancer gene expression file
  elevated_genes_info = data.frame(
    gene_id = names(fasta[which(fold_changes > 1)]),
    fold_change = fold_changes[which(fold_changes > 1)],
    risk_level = risk_level,
    patient_id = patient_id
  )
  
  write.table(elevated_genes_info, 
              file = paste0("data/cancer_gene_expression/patient_", 
                           sprintf("%02d", patient_id), "_genes.tsv"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Calculate reads per transcript
  tx_lengths = width(fasta)
  readspertx = round(coverage * tx_lengths / readlen)
  
  # Create patient directory
  patient_dir = paste0("data/patient_", sprintf("%02d", patient_id), "_", risk_level)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Generate reads for each replicate
  for (rep in 1:num_replicates) {
    set.seed(patient_id * 1000 + rep)
    
    all_sequences = character()
    all_qualities = character()
    
    # Pool reads from all transcripts for this replicate
    for (tx_idx in 1:num_transcripts) {
      num_reads = round(readspertx[tx_idx] * fold_changes[tx_idx])
      
      if (num_reads > 0) {
        sequence = as.character(fasta[[tx_idx]])
        seq_length = nchar(sequence)
        
        # Generate reads from random positions
        start_positions = sample(1:(seq_length - readlen + 1), num_reads, replace = TRUE)
        
        for (pos in start_positions) {
          read_seq = substr(sequence, pos, pos + readlen - 1)
          all_sequences = c(all_sequences, read_seq)
          all_qualities = c(all_qualities, paste(rep("I", readlen), collapse = ""))
        }
      }
    }
    
    # Write pooled FASTQ file
    output_file = paste0(patient_dir, "/sample_", sprintf("%02d", rep), ".fq")
    write_fastq_records(all_sequences, all_qualities, output_file)
    
    num_total_reads = length(all_sequences)
    print(paste("Generated patient", sprintf("%02d", patient_id), 
                "- Risk level:", risk_level,
                "- Replicate:", rep,
                "- Reads:", num_total_reads,
                "- Coverage:", coverage))
  }
  
  # Increase coverage every 8 patients
  if ((patient_id %% 8) == 0) {
    coverage = coverage + 15
    print(paste("Increased coverage to:", coverage))
  }
}

print("Dataset generation complete!")
print("Summary:")
print(paste("- Total patients:", num_patients))
print(paste("- Output structure: data/patient_XX_[risk_level]/"))
print(paste("- Cancer gene information: data/cancer_gene_expression/"))
print(paste("- Expected runtime for salmon/kallisto: ~15-25 minutes"))
print(paste("- Expected runtime for needle: ~2-5 minutes"))
