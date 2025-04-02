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
library(pheatmap)

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

snp_hk  <-subset(snp_hk, mutation_type == c("N"))

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

#for + strand
hk_data$ref_aa <- hk_data$codon
hk_data$ref_aa_rounded <- ceiling(hk_data$ref_aa+ 0.5)

hk_data$extracted_aa <- mapply(function(seq, pos) substring(seq, pos, pos), 
                               hk_data$protein_seq, hk_data$ref_aa_rounded)

# Extract the first letter after 'N:' from the mutation column
hk_data$mutated_aa <- substring(hk_data$mutation, 3, 3)

# Compare the extracted amino acid with the mutated amino acid
hk_data$match <- hk_data$extracted_aa == hk_data$mutated_aa

# Optionally, convert the match column to a more user-friendly format (TRUE/FALSE)
hk_data$match <- ifelse(hk_data$match, "Match", "No Match")

hk_data$new_aa <- substring(hk_data$mutation, nchar(hk_data$mutation), nchar(hk_data$mutation))

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


### Write amino acid sequence for each location with all unique occuring mutations

# Define the function to apply mutations
apply_mutations_to_sequence <- function(protein_sequence, mutation_df) {
  # Convert the protein sequence to a vector of individual amino acids
  protein_sequence_vec <- unlist(strsplit(protein_sequence, ""))
  
  # Print original protein sequence
  print(paste("Original Protein Sequence:", protein_sequence))
  
  # Apply mutations
  for (i in 1:nrow(mutation_df)) {
    position <- mutation_df$ref_aa_rounded[i]  # Get the position for the mutation
    new_aa <- mutation_df$new_aa[i]             # Get the new amino acid
    
    # Check if position is valid
    if (position > 0 && position <= length(protein_sequence_vec)) {
      protein_sequence_vec[position] <- new_aa  # Apply the mutation
    } else {
      warning(paste("Skipping mutation at position:", position, "- Position exceeds sequence length"))
    }
  }
  
  # Recombine the vector back into a string after all mutations
  mutated_sequence <- paste(protein_sequence_vec, collapse = "")
  
  # Return the mutated protein sequence
  return(mutated_sequence)
}

# Rockefeller

rock_subset <- hk_data[grepl("^Rock", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
rock_mutations <- unique(rock_subset$mutation_def)

rock_subset  <- rock_subset  %>% 
  filter(mutation_def %in% rock_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- rock_subset

# Call the function
rock_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-Rockefeller.txt"

# Save the sequence to a text file
writeLines(rock_protein_sequence , file_path)


# Zurich

zur_subset <- hk_data[grepl("^Zur", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
zur_mutations <- unique(zur_subset$mutation_def)

zur_subset  <- zur_subset  %>% 
  filter(mutation_def %in% zur_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- zur_subset

# Call the function
zur_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-Zurich.txt"

# Save the sequence to a text file
writeLines(zur_protein_sequence , file_path)


# LO

lo_subset <- hk_data[grepl("^LO", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
lo_mutations <- unique(lo_subset$mutation_def)

lo_subset  <- lo_subset  %>% 
  filter(mutation_def %in% lo_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- lo_subset

# Call the function
lo_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-LO.txt"

# Save the sequence to a text file
writeLines(lo_protein_sequence , file_path)


# AN

an_subset <- hk_data[grepl("^AN", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
an_mutations <- unique(an_subset$mutation_def)

an_subset  <- an_subset  %>% 
  filter(mutation_def %in% an_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- an_subset

# Call the function
an_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-AN.txt"

# Save the sequence to a text file
writeLines(an_protein_sequence , file_path)


# SL

sl_subset <- hk_data[grepl("^SL", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
sl_mutations <- unique(sl_subset$mutation_def)

sl_subset  <- sl_subset  %>% 
  filter(mutation_def %in% sl_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- sl_subset

# Call the function
sl_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-SL.txt"

# Save the sequence to a text file
writeLines(sl_protein_sequence , file_path)

# CB

cb_subset <- hk_data[grepl("^CB", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
cb_mutations <- unique(cb_subset$mutation_def)

cb_subset  <- cb_subset  %>% 
  filter(mutation_def %in% cb_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- cb_subset

# Call the function
cb_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-CB.txt"

# Save the sequence to a text file
writeLines(cb_protein_sequence , file_path)

# DB

db_subset <- hk_data[grepl("^DB", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
db_mutations <- unique(db_subset$mutation_def)

db_subset  <- db_subset  %>% 
  filter(mutation_def %in% db_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- db_subset

# Call the function
db_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-DB.txt"

# Save the sequence to a text file
writeLines(db_protein_sequence , file_path)

# ES

es_subset <- hk_data[grepl("^ES", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
es_mutations <- unique(es_subset$mutation_def)

es_subset  <- es_subset  %>% 
  filter(mutation_def %in% es_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- es_subset

# Call the function
es_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-ES.txt"

# Save the sequence to a text file
writeLines(es_protein_sequence , file_path)

# MC

mc_subset <- hk_data[grepl("^MC", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
mc_mutations <- unique(mc_subset$mutation_def)

mc_subset  <- mc_subset  %>% 
  filter(mutation_def %in% mc_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- mc_subset

# Call the function
mc_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-MC.txt"

# Save the sequence to a text file
writeLines(mc_protein_sequence , file_path)


# Greece

gr_subset <- hk_data[grepl("^Greece", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
gr_mutations <- unique(gr_subset$mutation_def)

gr_subset  <- gr_subset  %>% 
  filter(mutation_def %in% gr_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- gr_subset

# Call the function
gr_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-Greece.txt"

# Save the sequence to a text file
writeLines(gr_protein_sequence , file_path)

# Istanbul

is_subset <- hk_data[grepl("^Istan", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
is_mutations <- unique(is_subset$mutation_def)

is_subset  <- is_subset  %>% 
  filter(mutation_def %in% is_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- is_subset

# Call the function
is_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-Istanbul.txt"

# Save the sequence to a text file
writeLines(is_protein_sequence , file_path)


# Leiden

lei_subset <- hk_data[grepl("^Leiden", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
lei_mutations <- unique(lei_subset$mutation_def)

lei_subset  <- lei_subset  %>% 
  filter(mutation_def %in% lei_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- lei_subset

# Call the function
lei_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-Leiden.txt"

# Save the sequence to a text file
writeLines(lei_protein_sequence , file_path)


# NAA

naa_subset <- hk_data[grepl("^NAA", hk_data$Location_ID), ]

# Sample the unique mutations at specific locaiton
naa_mutations <- unique(naa_subset$mutation_def)

naa_subset  <- naa_subset  %>% 
  filter(mutation_def %in% naa_mutations) %>% 
  distinct(mutation_def, .keep_all=TRUE)

mutation_df <- naa_subset

# Call the function
naa_protein_sequence <- apply_mutations_to_sequence(protein_sequences, mutation_df)

# Specify the file path (you can change the path)
file_path <- "Amidophosphoribosyltransferase/Amidophosphoribosyltransferase-NAA.txt"

# Save the sequence to a text file
writeLines(naa_protein_sequence , file_path)
### Prepare data for heatmap

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

write.xlsx(hk_data, "data/amidophosphoribosyltransferase.xlsx",rowNames = FALSE)

qsave(presence_df, "data/detected_mutations.qs",
      preset = "archive")

qsave(hk_data, "data/hk_data.qs",
      preset = "archive")
