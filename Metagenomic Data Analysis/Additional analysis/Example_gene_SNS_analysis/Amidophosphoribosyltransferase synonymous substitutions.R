# load packages
library(ggplot2)
library(readxl)
library(reshape2)
library(dplyr)
library(readr)
library(data.table)
library(stringr)
library(openxlsx)
library(tidyr)
library(plyr)
library(purrr)
library(Biostrings)

# Genome annotation
genome <- read_excel("data/gene_annotation_features.xlsx")

# SNP annotation
snp <- read.table("data/final_filtered_data_with_hypothetical_and_additional_annotation.csv", header = TRUE, sep = ",")

heatmap<- read.table("data/genome_final_list_wide_with_product.csv", header = TRUE, sep = ",")

hk <- subset(heatmap, product == c("Amidophosphoribosyltransferase"))

# Step 1: Select only numeric columns
numeric_hk <- hk[sapply(hk, is.numeric)]

# Step 2: Count the values greater than 1 in each row
count_over_1 <- rowSums(numeric_hk > 1, na.rm = TRUE)

# Step 3: Order the rows by the count of values greater than 1, in descending order
ranked_indices <- order(count_over_1, decreasing = TRUE)

# Step 4: Extract the genes in that order
ranked_genes <- hk$gene[ranked_indices]

# Step 5 (optional): Create a data frame showing the genes and their corresponding counts
ranked_genes_df <- data.frame(
  gene = hk$gene[ranked_indices],
  count_over_1 = count_over_1[ranked_indices]
)

gene_of_interest <- ranked_genes_df$gene[1]


# Subset SNP annotation for gene of interest

snp_hk <- subset(snp, gene== c(gene_of_interest))

snp_hk <- subset(snp_hk, class== c("SNS"))

snp_hk  <-subset(snp_hk, mutation_type == c("S"))

# Read the .faa file
faa_file <- "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/Bacteroides sp002491635/data/GUT_GENOME000231_6.faa"  # Replace with the actual path to your .faa file
sequences <- readAAStringSet(faa_file)

sequences_headers <- names(sequences)
sequences_headers<- as.data.frame(sequences_headers)


# Get protein identifier 

gene_id <- subset(genome, gene_id==(gene_of_interest))
protein_id<- gene_id$locus_tag

index <- grep(protein_id, names(sequences))
selected_sequence <- sequences[index]

sequence_names <- names(selected_sequence)
protein_sequences <- as.character(selected_sequence)

# Add protein sequence to snp file

snp_hk$protein_seq <- protein_sequences  
snp_hk$protein_id <- protein_id
snp_hk$start_position  <- gene_id$start_position
snp_hk$end_position  <- gene_id$end_position


hk_data <- snp_hk %>%
  select(protein_id,start_position,end_position,position,Location_ID,mutation,protein_seq,new_sample_label)

hk_data$pos_diff <- hk_data$position - hk_data$start_position
hk_data$codon <- hk_data$pos_diff / 3

all_aa<- nchar(hk_data$protein_seq[1])

# Get the total number of amino acids
total_amino_acids <- nchar(protein_sequences)
num_aa <- as.numeric(total_amino_acids)

#for + strand
hk_data$ref_aa <- hk_data$codon
hk_data$ref_aa_rounded <- ceiling(hk_data$ref_aa+ 0.5)

hk_data$extracted_aa <- mapply(function(seq, pos) substring(seq, pos, pos), 
                               hk_data$protein_seq, hk_data$ref_aa_rounded)

# Extract the first letter after 'N:' from the mutation column
hk_data$mutated_aa <- hk_data$extracted_aa 

# Compare the extracted amino acid with the mutated amino acid
hk_data$match <- hk_data$extracted_aa == hk_data$mutated_aa

# Optionally, convert the match column to a more user-friendly format (TRUE/FALSE)
hk_data$match <- ifelse(hk_data$match, "Match", "No Match")

hk_data$new_aa <-hk_data$extracted_aa 

# Function to mutate the sequence
mutate_sequence <- function(seq, pos, new_aa) {
  # Replace the amino acid at the given position
  if (pos <= nchar(seq)) {
    return(paste0(substring(seq, 1, pos - 1), new_aa, substring(seq, pos + 1)))
  } else {
    return(seq)  # If position is out of range, return the original sequence
  }
}

# Apply the function to create the mutated_seq column
hk_data$mutated_seq <- mapply(mutate_sequence, hk_data$protein_seq, hk_data$ref_aa_rounded, hk_data$new_aa)

hk_data$mutation_def <- paste("Position:", hk_data$ref_aa_rounded,
                              hk_data$mutated_aa, 
                              "->", 
                              hk_data$new_aa)


# Step 1: Get unique mutation definitions
unique_sequences <- unique(hk_data$mutation_def)

# Step 2: Create a new dataframe with unique mutation_defs as rows
presence_df <- data.frame(mutation_def = unique_sequences)

# Step 3: Get unique Location_IDs
location_ids <- unique(hk_data$new_sample_label)

# Step 4: Initialize columns for each Location_ID with 0
for (loc in location_ids) {
  presence_df[[loc]] <- 0  # Set initial values to 0
}

# Step 5: Populate the new dataframe with 1/0 based on mutation_def
for (i in 1:nrow(presence_df)) {
  seq <- presence_df$mutation_def[i]
  
  # Check for presence of the mutation_def in the original data
  present_locs <- unique(hk_data$new_sample_label[hk_data$mutation_def == seq])
  
  # Mark 1 for those Location_IDs present for this mutation_def
  presence_df[i, present_locs] <- 1
}

write.xlsx(hk_data, "data/Amidophosphoribosyltransferase synonymous substitutions.xlsx",rowNames = FALSE)



qsave(presence_df, "data/detected_mutations_synonymous.qs",
      preset = "archive")

qsave(hk_data, "data/hk_data_synonymous.qs",
      preset = "archive")




