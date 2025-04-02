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
library(broom)
library(dunn.test)
library(qs)

### Adapt these before running the script
genome_of_interest <- c("GUT_GENOME000221")
species_of_interest <- c("s__Bacteroides acidifaciens")
###

# Store current working directory
current_dir <- getwd()

# Select all contig identifiers
contigs <- read_excel("data/gene_annotation_features.xlsx")
# Extract unique text from both columns
unique_contigs <- unique(c(contigs$scaffold))

# Specify the directory where the file is located
base_dir <- "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table"

# Change working directory to the input directory
setwd(base_dir)

# Select all samples of interest
sample_data <- read_excel("IS.COMPARE_All_genomeWide_compare_Subset.xlsx")

# Subset according to species of interest
sample_data <- subset(sample_data, sample_data$`Strain Name`== species_of_interest)

# Extract unique text from both columns
unique_samples <- unique(c(
  gsub("\\.sorted\\.bam$", "", sample_data$name1),
  gsub("\\.sorted\\.bam$", "", sample_data$name2)
))

#Create metadata file for future analysis with the correct samples
metadata_all <- read.table("metadata_all.txt", header = TRUE, sep = "\t")

# Subset the data frame based on the Sample column
metadata <- subset(metadata_all, sample %in% unique_samples)                                                    

metadata <- metadata%>%
  mutate(Location = paste(country_abbreviation, institute_abbreviation, sep = "_"))

# Subset the dataframe to keep only specific columns using dplyr
metadata <- metadata %>%
  select(sample,Groupv2, Location, Gender, SamplingLocation,Young_Adult,sample_id,new_sample_label)

names(metadata)[names(metadata) == "Groupv2"] <- "Group"
names(metadata)[names(metadata) == "Young_Adult"] <- "Age"

# Save the updated DataFrame to a new CSV file with the name of your desired genome
write.csv(metadata, "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/Bacteroides acidifaciens/data/metadata.csv", row.names = FALSE)

qsave(metadata, "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/Bacteroides acidifaciens/data/metadata.qs",
      preset = "archive")

# List of specific folders to process #1
folders <- unique_samples

# Specify the directory where the file is located
base_dir <- "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis/combined_output/combined_output"

# Change working directory to the input directory
setwd(base_dir)

# Initialize an empty list to store data frames
data_list <- list()

# Loop through each folder and read data into individual data frames
for (folder in folders) {
  is_folder <- file.path(base_dir, paste0(folder, ".IS"))
  tsv_file <- file.path(is_folder, "output", paste0(folder, ".IS_gene_info.tsv"))
  
  if (dir.exists(is_folder) && file.exists(tsv_file)) {
    df <- fread(tsv_file, sep = "\t")
    df$Sample <- folder
    data_list[[folder]] <- df
  } else {
    message(paste("Folder or file does not exist:", is_folder, tsv_file))
  }
}

# Combine all data frames into one
combined_data <- bind_rows(data_list)

####  Filter the combined data
filtered_data <- combined_data %>%
  filter(scaffold %in% unique_contigs) %>%
  filter(coverage > 5, breadth > 0.8)

# Get unique genes
unique_gene <- filtered_data %>% distinct(gene, .keep_all = TRUE)

# return to the previous working directory
setwd(current_dir)

# Match and extract the unique genes for each folder, then combine the results
combined_unique_genes <- lapply(folders, function(folder) {
  df <- data_list[[folder]]
  df_gene <- df[match(unique_gene$gene, df$gene), ]
  return(df_gene)
})

# Combine the results into a single data frame
unique_genes_genome <- do.call(rbind, combined_unique_genes)

write.xlsx(unique_genes_genome, "data/unique_genes_genome.xlsx",rowNames = FALSE)

count_gene<-unique_genes_genome %>% 
  dplyr::count(gene)
count_gene <- count_gene  %>%  
  filter(n>20) # See if you have enough data after counting

unique_genes_genome_final <- unique_genes_genome  %>%  
  filter(gene %in% count_gene$gene)

write.xlsx(unique_genes_genome_final, "data/unique_genes_genome_final.xlsx",rowNames = FALSE)

setwd(current_dir)
### Create first data for dN/dS Heatmap

unique_genes_genome_final_list <- unique_genes_genome_final %>% dplyr::select(gene, dNdS_substitutions, Sample)

histogram <-ggplot(unique_genes_genome_final_list, aes(x = dNdS_substitutions)) +
  geom_histogram(binwidth = 0.1, fill = "grey",color="black", alpha = 0.7) +
  labs(title = "Distribution of dN/dS", x = "dN/dS Value", y = "Frequency") +
  theme_minimal()

histogram

ggsave("Plots/Histogram dN over dS.pdf", plot = histogram, width = 5, height = 3, useDingbats = FALSE, limitsize = FALSE)

location<- metadata  %>%
  select(sample,Location,Group)

names(location)[names(location) == "sample"] <- "Sample"

data_overview <-left_join(unique_genes_genome_final_list, location, by="Sample")

# Remove rows where location is NA
data_overview<- data_overview %>%
  filter(!is.na(Location))

histogram <- ggplot(data_overview, aes(x = dNdS_substitutions, fill = Group,color=Group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7) +
  labs(title = "Distribution of dN/dS (SPF vs Wild)", x = "dN/dS Value", y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ Location) +
  scale_fill_manual(values = c("SPF" = "#003399", "WildMice" = "#99CCFF"))+ # Custom color for outline
  scale_color_manual(values = c("SPF" = "black", "WildMice" = "black"))

histogram

ggsave("Plots/Histogram different vivaria dN over dS.pdf", plot = histogram, width = 10, height = 5, useDingbats = FALSE, limitsize = FALSE)


# Make summary statistics for different vivaria

data_overview_summary <- data_overview %>%
  group_by(gene, Location) %>%
  dplyr::summarise(
    mean_dN_dS = ifelse(all(is.na(dNdS_substitutions)), NA_real_, mean(dNdS_substitutions, na.rm = TRUE)),
    median_dN_dS = ifelse(all(is.na(dNdS_substitutions)), NA_real_, median(dNdS_substitutions, na.rm = TRUE)),
    sd_dN_dS = ifelse(all(is.na(dNdS_substitutions)), NA_real_, sd(dNdS_substitutions, na.rm = TRUE)),
    min_dN_dS = ifelse(all(is.na(dNdS_substitutions)), NA_real_, min(dNdS_substitutions, na.rm = TRUE)),
    max_dN_dS = ifelse(all(is.na(dNdS_substitutions)), NA_real_, max(dNdS_substitutions, na.rm = TRUE)),
    count = dplyr::n(),  # Number of samples in each group
    .groups = 'drop'  # Drop grouping structure after summarisation
  ) %>%
  ungroup()


# Remove rows where all summary variables are NA
data_overview_summary_clean <- data_overview_summary %>%
  filter(!is.na(mean_dN_dS) | !is.na(median_dN_dS) | !is.na(sd_dN_dS) | !is.na(min_dN_dS) | !is.na(max_dN_dS))

write.xlsx(data_overview_summary_clean, "data/dN_dS_summary_statistics_clean.xlsx",rowNames = FALSE)

# Function to perform Kruskal-Wallis test for each gene
kruskal_analysis <- function(df) {
  # Filter out locations with only missing values
  df_filtered <- df %>%
    group_by(Location) %>%
    filter(!all(is.na(dNdS_substitutions))) %>%
    ungroup()
  
  # Check if we have enough locations
  if (n_distinct(df_filtered$Location) < 2) {
    return(data.frame(
      gene = unique(df$gene),
      kruskal_p_value = NA
    ))
  }
  
  # Perform Kruskal-Wallis test
  kruskal_result <- tryCatch({
    kruskal.test(dNdS_substitutions ~ Location, data = df_filtered)
  }, error = function(e) {
    list(p.value = NA)
  })
  kruskal_p_value <- kruskal_result$p.value
  
  # Return results
  return(data.frame(
    gene = unique(df$gene),
    kruskal_p_value = kruskal_p_value
  ))
}

kruskal_results <- data_overview %>%
  group_by(gene) %>%
  group_split() %>%
  lapply(kruskal_analysis) %>%
  bind_rows()

kruskal_p_values <- kruskal_results %>%
  filter(!is.na(kruskal_p_value)) %>%
  select(gene, kruskal_p_value)

# Adjust Kruskal-Wallis p-values for multiple comparisons
adjusted_kruskal_p_values <- p.adjust(kruskal_p_values$kruskal_p_value, method = "holm")

kruskal_results <- kruskal_results %>%
  left_join(data.frame(gene = kruskal_p_values$gene, adjusted_kruskal_p_value = adjusted_kruskal_p_values), by = "gene")

kruskal_results <- kruskal_results %>%
  filter(!is.na(kruskal_p_value))

final_data_overview_summary <- data_overview_summary %>%
  left_join(kruskal_results, by = "gene")


# Remove rows where all summary variables are NA
final_data_overview_summary <- final_data_overview_summary %>%
  filter(!is.na(mean_dN_dS) | !is.na(median_dN_dS) | !is.na(sd_dN_dS) | !is.na(min_dN_dS) | !is.na(max_dN_dS))

write.xlsx(final_data_overview_summary, "data/dN_dS_summary_statistics_with_kruskal.xlsx",rowNames = FALSE)

# Identify all duplicate rows
duplicates_all <- unique_genes_genome_final_list %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))

# Check if duplicates are all NA
na_only_duplicates <- duplicates_all %>%
  filter(across(everything(), ~ is.na(.))) %>%
  distinct()  # Ensure unique rows

# Remove these NA-only duplicate rows from the original data
cleaned_data <- unique_genes_genome_final_list %>%
  anti_join(na_only_duplicates, by = colnames(na_only_duplicates))

# Identify all duplicate rows if there are still some
duplicates_all <- cleaned_data %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))

unique_genes_genome_final_list_wide <- cleaned_data %>%
  pivot_wider(names_from = Sample, values_from = dNdS_substitutions)

write.xlsx(unique_genes_genome_final_list_wide, "data/unique_genes_genome_final_list_wide_dN_dS.xlsx",rowNames = FALSE)

# Associate the data with the gene annotation
genome_data <- contigs %>%
  select(gene_id,product)

# Rename ID_complete to gene in genome_data for joining
genome_data <- genome_data %>%
  dplyr::rename(gene = gene_id)

genome_data <- genome_data %>%
  mutate(gene = as.character(gene))

# Perform the join operation
gut_genome_data <- unique_genes_genome_final_list_wide %>%
  dplyr::left_join(genome_data, by = "gene")

# Replace NA values in the product column with 'No product information'
gut_genome_data$product[is.na(gut_genome_data$product)] <- 'No product information'

# Remove rows where the 'gene' column is NA
gut_genome_data <- gut_genome_data[!is.na(gut_genome_data$gene), ]

# Save the updated DataFrame to a new CSV file with the name of your desired genome
write.csv(gut_genome_data, "data/genome_final_list_wide_with_product.csv", row.names = FALSE)

# Save dataframe for heatmap
subset_data <- gut_genome_data %>%
  select(-gene)

write.csv(subset_data, "data/data_for_heatmap_dN_dS.csv", row.names = FALSE)

qsave(subset_data, "data/data_for_heatmap_dN_dS.qs",
      preset = "archive")



#Now the same for pN over pS

unique_genes_genome_final_list <- unique_genes_genome_final %>% dplyr::select(gene, pNpS_variants, Sample)

histogram <-ggplot(unique_genes_genome_final_list, aes(x = pNpS_variants)) +
  geom_histogram(binwidth = 0.1, fill = "grey", color = "black", alpha = 0.7) +
  labs(title = "Distribution of pN/pS", x = "pN/pS Value", y = "Frequency") +
  theme_minimal()

histogram

ggsave("Plots/Histogram pN over pS.pdf", plot = histogram, width = 5, height = 3, useDingbats = FALSE, limitsize = FALSE)

data_overview <-left_join(unique_genes_genome_final_list, location, by="Sample")

# Remove rows where location is NA
data_overview<- data_overview %>%
  filter(!is.na(Location))

histogram <- ggplot(data_overview, aes(x = pNpS_variants, fill = Group,color=Group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7) +
  labs(title = "Distribution of pN/pS (SPF vs Wild)", x = "pN/pS Value", y = "Frequency") +
  theme_minimal() +
  facet_wrap(~ Location) +
  scale_fill_manual(values = c("SPF" = "#003399", "WildMice" = "#99CCFF"))+ # Custom color for outline
  scale_color_manual(values = c("SPF" = "black", "WildMice" = "black"))

histogram

ggsave("Plots/Histogram different vivaria pN over pS.pdf", plot = histogram, width = 10, height = 5, useDingbats = FALSE, limitsize = FALSE)


# Make summary statistics for different vivaria

data_overview_summary <- data_overview %>%
  group_by(gene, Location) %>%
  dplyr::summarise(
    mean_pN_pS = ifelse(all(is.na(pNpS_variants)), NA_real_, mean(pNpS_variants, na.rm = TRUE)),
    median_pN_pS = ifelse(all(is.na(pNpS_variants)), NA_real_, median(pNpS_variants, na.rm = TRUE)),
    sd_pN_pS = ifelse(all(is.na(pNpS_variants)), NA_real_, sd(pNpS_variants, na.rm = TRUE)),
    min_pN_pS = ifelse(all(is.na(pNpS_variants)), NA_real_, min(pNpS_variants, na.rm = TRUE)),
    max_pN_pS = ifelse(all(is.na(pNpS_variants)), NA_real_, max(pNpS_variants, na.rm = TRUE)),
    count = dplyr::n(),  # Number of samples in each group
    .groups = 'drop'  # Drop grouping structure after summarisation
  ) %>%
  ungroup()


# Remove rows where all summary variables are NA
data_overview_summary_clean <- data_overview_summary %>%
  filter(!is.na(mean_pN_pS) | !is.na(median_pN_pS) | !is.na(sd_pN_pS) | !is.na(min_pN_pS) | !is.na(max_pN_pS))

write.xlsx(data_overview_summary_clean, "data/pN_pS_summary_statistics_clean.xlsx",rowNames = FALSE)

# Function to perform Kruskal-Wallis test for each gene
kruskal_analysis <- function(df) {
  # Filter out locations with only missing values
  df_filtered <- df %>%
    group_by(Location) %>%
    filter(!all(is.na(pNpS_variants))) %>%
    ungroup()
  
  # Check if we have enough locations
  if (n_distinct(df_filtered$Location) < 2) {
    return(data.frame(
      gene = unique(df$gene),
      kruskal_p_value = NA
    ))
  }
  
  # Perform Kruskal-Wallis test
  kruskal_result <- tryCatch({
    kruskal.test(pNpS_variants ~ Location, data = df_filtered)
  }, error = function(e) {
    list(p.value = NA)
  })
  kruskal_p_value <- kruskal_result$p.value
  
  # Return results
  return(data.frame(
    gene = unique(df$gene),
    kruskal_p_value = kruskal_p_value
  ))
}

kruskal_results <- data_overview %>%
  group_by(gene) %>%
  group_split() %>%
  lapply(kruskal_analysis) %>%
  bind_rows()

kruskal_p_values <- kruskal_results %>%
  filter(!is.na(kruskal_p_value)) %>%
  select(gene, kruskal_p_value)

# Adjust Kruskal-Wallis p-values for multiple comparisons
adjusted_kruskal_p_values <- p.adjust(kruskal_p_values$kruskal_p_value, method = "holm")

kruskal_results <- kruskal_results %>%
  left_join(data.frame(gene = kruskal_p_values$gene, adjusted_kruskal_p_value = adjusted_kruskal_p_values), by = "gene")

kruskal_results <- kruskal_results %>%
  filter(!is.na(kruskal_p_value))

final_data_overview_summary <- data_overview_summary %>%
  left_join(kruskal_results, by = "gene")


# Remove rows where all summary variables are NA
final_data_overview_summary <- final_data_overview_summary %>%
  filter(!is.na(mean_pN_pS) | !is.na(median_pN_pS) | !is.na(sd_pN_pS) | !is.na(min_pN_pS) | !is.na(max_pN_pS))

write.xlsx(final_data_overview_summary, "data/pN_pS_summary_statistics_with_kruskal.xlsx",rowNames = FALSE)


# Identify all duplicate rows
duplicates_all <- unique_genes_genome_final_list %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))

# Check if duplicates are all NA
na_only_duplicates <- duplicates_all %>%
  filter(across(everything(), ~ is.na(.))) %>%
  distinct()  # Ensure unique rows

# Remove these NA-only duplicate rows from the original data
cleaned_data <- unique_genes_genome_final_list %>%
  anti_join(na_only_duplicates, by = colnames(na_only_duplicates))

# Identify all duplicate rows
duplicates_all <- cleaned_data %>%
  filter(duplicated(.) | duplicated(., fromLast = TRUE))

unique_genes_genome_final_list_wide <- cleaned_data %>%
  pivot_wider(names_from = Sample, values_from = pNpS_variants)

write.xlsx(unique_genes_genome_final_list_wide, "data/unique_genes_genome_final_list_wide_pNpS.xlsx",rowNames = FALSE)

# Associate the data with the gene annotation
genome_data <- contigs %>%
  select(gene_id,product)

# Rename ID_complete to gene in genome_data for joining
genome_data <- genome_data %>%
  dplyr::rename(gene = gene_id)

genome_data <- genome_data %>%
  mutate(gene = as.character(gene))

# Perform the join operation
gut_genome_data <- unique_genes_genome_final_list_wide %>%
  dplyr::left_join(genome_data, by = "gene")

# Replace NA values in the product column with 'No product information'
gut_genome_data$product[is.na(gut_genome_data$product)] <- 'No product information'

# Remove rows where the 'gene' column is NA
gut_genome_data <- gut_genome_data[!is.na(gut_genome_data$gene), ]

# Remove columns that contain only zeros
gut_genome_data <- gut_genome_data %>% select_if(~!all(. == 0))

# Save the updated DataFrame to a new CSV file with the name of your desired genome
write.csv(gut_genome_data, "data/genome_final_list_wide_with_product_pNpS.csv", row.names = FALSE)

# Save dataframe for heatmap
subset_data <- gut_genome_data %>%
  select(-gene)

write.csv(subset_data, "data/data_for_heatmap_pN_pS.csv", row.names = FALSE)

qsave(subset_data, "data/data_for_heatmap_pN_pS.qs",
      preset = "archive")


### Check for gene of interest

unique_genes_genome <- subset(unique_genes_genome, gene == c("GUT_GENOME000221_54_12"))
unique_genes_genome <- subset(unique_genes_genome, dNdS_substitutions > 1)



                      