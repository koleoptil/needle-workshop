#!/usr/bin/env Rscript
library(Biostrings)
library(edgeR)
library(polyester)

# Script to create a simulated RNA-seq dataset for needle workshop
# Compares needle (fast) vs salmon/kallisto (slower) for expression estimation
# Includes cancer-associated genes with varying risk profiles across patients

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input FASTA file).n", call.=FALSE)
}

fastapath = args[1]

# Read the reference transcriptome
fasta = readDNAStringSet(fastapath)
num_transcripts = length(fasta)

# Define cancer-associated genes (we'll simulate these by selecting random transcripts)
num_cancer_genes = max(50, round(0.05 * num_transcripts))  # 5% of transcripts as "cancer genes"
cancer_gene_indices = sample(1:num_transcripts, num_cancer_genes)

# Simulation parameters
num_patients = 24  # Creates 24 patient samples (reasonable for a workshop)
num_replicates = 2  # Technical replicates per patient
coverage = 30      # Coverage depth (higher = more data, longer runtime for salmon/kallisto)
bp_length = 75     # Read length
readlen = 75

# Risk categories:
# - 8 low-risk patients: 0-1 cancer genes with high expression
# - 8 medium-risk patients: 2-3 cancer genes with high expression
# - 8 high-risk patients: 4-6 cancer genes with high expression

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

for (patient_id in 1:num_patients) {
  risk_level = risk_assignments[patient_id]
  
  # Initialize fold changes (baseline = 1 in all samples)
  fold_changes = matrix(rep(1, num_transcripts * num_replicates), 
                        nrow = num_transcripts, 
                        ncol = num_replicates)
  
  # Assign expression patterns based on risk level
  if (risk_level == "low") {
    # 0-1 cancer genes with moderate elevation (2-3x)
    num_elevated = sample(0:1, 1)
    if (num_elevated > 0) {
      elevated_genes = sample(cancer_gene_indices, num_elevated)
      fold_change_values = sample(c(2, 3), num_elevated, replace = TRUE)
      for (i in 1:num_elevated) {
        fold_changes[elevated_genes[i], ] = fold_change_values[i]
      }
    }
  } else if (risk_level == "medium") {
    # 2-3 cancer genes with higher elevation (3-6x)
    num_elevated = sample(2:3, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_change_values = sample(c(3, 4, 5, 6), num_elevated, replace = TRUE)
    for (i in 1:num_elevated) {
      fold_changes[elevated_genes[i], ] = fold_change_values[i]
    }
  } else if (risk_level == "high") {
    # 4-6 cancer genes with very high elevation (5-10x)
    num_elevated = sample(4:6, 1)
    elevated_genes = sample(cancer_gene_indices, num_elevated)
    fold_change_values = sample(c(5, 6, 7, 8, 9, 10), num_elevated, replace = TRUE)
    for (i in 1:num_elevated) {
      fold_changes[elevated_genes[i], ] = fold_change_values[i]
    }
  }
  
  # Write cancer gene expression file (for downstream analysis)
  elevated_genes_info = data.frame(
    gene_id = names(fasta[which(fold_changes[, 1] > 1)]),
    fold_change = fold_changes[which(fold_changes[, 1] > 1), 1],
    risk_level = risk_level,
    patient_id = patient_id
  )
  write.table(elevated_genes_info, 
              file = paste0("data/cancer_gene_expression/patient_", 
                           sprintf("%02d", patient_id), "_genes.tsv"),
              col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Calculate reads per transcript based on coverage
  readspertx = round(coverage * width(fasta) / bp_length)
  
  # Create patient-specific output directory
  patient_dir = paste0("data/patient_", sprintf("%02d", patient_id), "_", risk_level)
  dir.create(patient_dir, showWarnings = FALSE)
  
  # Simulate reads for this patient (with replicates)
  simulate_experiment(
    fastapath,
    reads_per_transcript = readspertx,
    num_reps = c(num_replicates, num_replicates),  # Two conditions with replicates each
    outdir = patient_dir,
    fold_changes = fold_changes,
    seed = (patient_id * 142),  # Different seed per patient
    readlen = readlen,
    paired = TRUE,
    gzip = TRUE
  )
  
  print(paste("Generated patient", sprintf("%02d", patient_id), 
              "- Risk level:", risk_level,
              "- Coverage:", coverage))
  
  # Increase coverage every 8 patients to make the dataset progressively larger
  # This ensures later patients have more data, increasing computational load
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
print(paste("- Expected runtime for salmon/kallisto: ~30 minutes"))
print(paste("- Expected runtime for needle: ~2-5 minutes"))
