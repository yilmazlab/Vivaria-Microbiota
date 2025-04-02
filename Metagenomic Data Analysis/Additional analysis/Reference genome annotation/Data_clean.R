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

### Adapt these before running the script
genome_of_interest <- c("GUT_GENOME000221")
species_of_interest <- c("s__Bacteroides_acidifaciens")
###

# Store current working directory
current_dir <- getwd()

### Cleaning of annotation table
# Import prokka annotation and example gene table

# Specify the directory where the file is located
base_dir <- "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table"

# Change working directory to the input directory
setwd(base_dir)

gene_table_example <- read.table("genes_table_Sweden4.csv", header = TRUE, sep = ",")
gene_table_subset <- subset(gene_table_example, grepl(paste0("^", genome_of_interest), gene_table_example$scaffold))
rm(gene_table_example)

# Rename columns
gene_table_subset <- gene_table_subset %>%
  dplyr::rename(strand = direction,start_position=start,end_position=end)

gene_table_subset <- gene_table_subset %>%
  select(-X)

gene_table_subset <- gene_table_subset %>%
  mutate(strand = recode(strand, `1` = "+", `-1` = "-"))

# Add 1 bp to all start and end positions to map gff file
gene_table_subset <- gene_table_subset %>%
  mutate(start_position = start_position + 1)

# Add 1 bp to all start and end positions to map gff file
gene_table_subset <- gene_table_subset %>%
  mutate(end_position = end_position + 1)

setwd(current_dir)

file_path <- "/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/Bacteroides acidifaciens/data/combined_annotations.gff"

# Check if the file exists
if (file.exists(file_path)) {
  cat("File exists.\n")
} else {
  stop("File does not exist at the specified path.")
}

# Read and display the first few lines of the file
file_lines <- readLines(file_path, n = 10)
cat(file_lines, sep = "\n")

# Read the GFF file using base R
genome <- read.table(file_path, sep = "\t", header = FALSE, comment.char = "#", fill = TRUE, quote = "")

# Step 1: Separate the key-value pairs into multiple rows
df_separated <- genome %>%
  separate_rows(V9, sep = ";")


# Step 2: Separate the key and value into different columns
df_key_value <- df_separated %>%
  separate(V9, into = c("key", "value"), sep = "=")

# Step 3: Filter out rows where the key is empty
df_filtered <- df_key_value %>%
  filter(key != "")

# Step 4: Pivot the data frame to wide format
genome_clean <- df_filtered %>%
  pivot_wider(names_from = key, values_from = value, values_fn = list)

# Split and keep first value

get_first_element <- function(x) {
  if (length(x) > 0) x[1] else NA
}

genome_clean$ID <- sapply(genome_clean$ID, get_first_element)
genome_clean$inference <- sapply(genome_clean$inference, get_first_element)
genome_clean$locus_tag <- sapply(genome_clean$locus_tag, get_first_element)
genome_clean$product <- sapply(genome_clean$product, get_first_element)
genome_clean$eC_number <- sapply(genome_clean$eC_number, get_first_element)
genome_clean$gene <- sapply(genome_clean$gene, get_first_element)
genome_clean$db_xref <- sapply(genome_clean$db_xref, get_first_element)

# Rename columns
genome_clean <- genome_clean %>%
  dplyr::rename(scaffold = V1,start_position=V4,end_position=V5, ftype=V3, strand=V7, source=V2)

prokka_annotation <- genome_clean 

prokka_annotation <- prokka_annotation %>%
  select(-V6,-V8, -rpt_family, -rpt_type, -rpt_unit_seq, -ID, -Name)

# Clean environment
rm(list = setdiff(ls(), c("prokka_annotation", "gene_table_subset","genome_of_interest","species_of_interest")))


matched_df <- prokka_annotation %>%
  inner_join(gene_table_subset, by = c("scaffold", "start_position", "end_position", "strand"))

matched_df <- matched_df %>%
  dplyr::rename(gene_name= gene.x,gene_id=gene.y)

matched_df  <- matched_df [, c("gene_id", setdiff(names(matched_df ), "gene_id"))]

gene_annotation <- matched_df

write.xlsx(gene_annotation, "data/gene_annotation_prokka.xlsx",rowNames = FALSE)

# Find rows in another_df that do not have a match in prokka_annotation
unmatched_df <- gene_table_subset %>%
  anti_join(prokka_annotation, by = c("scaffold", "start_position", "end_position", "strand"))

write.xlsx(unmatched_df, "data/genes_without_prokka_annotation.xlsx",rowNames = FALSE)

# Clean environment
rm(list = setdiff(ls(), c("gene_annotation","genome_of_interest","species_of_interest")))

### Import PFAM annotation
pfam_annotation <- read_excel("data/Pfam_results_annotated.xlsx")
pfam_annotation <- pfam_annotation %>%
  dplyr::rename(gene_id = gene,pfam_score =score,pfam_boundaries= boundaries, pfam_cond_evalue="cond-evalue", pfam_indp_evalue="indp-evalue", pfam_description=DESC,pfam_accession=ACC)

pfam_annotation <- pfam_annotation %>%
  select(-resolved)

# Filter for the best match of each gene_id based on the score
pfam_annotation  <- pfam_annotation  %>%
  group_by(gene_id) %>%
  filter(pfam_score == max(pfam_score)) %>%
  ungroup()

gene_annotation_features <- left_join(gene_annotation, pfam_annotation, by= "gene_id")

### Import CAZymes annotation
cazymes <- read.table("data/GUT_GENOME000221_CAZymes_results.csv", header = TRUE, sep = ",")

cazymes <- cazymes %>%
  dplyr::rename(gene_id=gene,CAZyme_class = class, CAZyme_family = family, CAZyme_subfamily = subfamily, CAZyme_coverage=Coverage, CAZyme_evalue=E.value)

cazymes <- cazymes %>%
  select(-HMM_length,-Query_length,-Query_start,-Query_end,-HMM_start,-HMM_end)

cazymes  <- cazymes  %>%
  group_by(gene_id) %>%
  filter(CAZyme_evalue == min(CAZyme_evalue)) %>%
  filter(CAZyme_coverage == max(CAZyme_coverage)) %>%
  ungroup()

gene_annotation_features <- left_join(gene_annotation_features, cazymes, by= "gene_id")

### Import KEGG annotation
kegg <- read_excel("data/kofamscan_results.xlsx")

kegg <- kegg %>%
  dplyr::rename(gene_id = gene,ko_score =score, ko_thrshld=thrshld, ko_evalue=e_value)

kegg <- kegg %>%
  group_by(gene_id) %>%
  filter(ko_score == max(ko_score)) %>%
  filter(ko_evalue == min(ko_evalue)) %>%
  ungroup()

gene_annotation_features <- left_join(gene_annotation_features, kegg, by= "gene_id")

### Import Antibiotic resistance annotation

card <- read_excel("data/card_results.xlsx")

card <- card %>%
  dplyr::rename(gene_id = gene,card_score =bit_score, card_evalue="e-value")

card <- card %>%
  select(-percentID,-alignment_length,-mm,-gaps,-querry_start,-querry_end,-target_start,-target_end,-protein_seq_accession,-target,-ARO_category_accessions)

gene_annotation_features <- left_join(gene_annotation_features, card, by= "gene_id")

write.xlsx(gene_annotation_features, "data/gene_annotation_features.xlsx",rowNames = FALSE)

## Check if products are matching with the additional feature annotation
compare_products <- gene_annotation_features %>%
  select(product,pfam_description,KO_definition,ARO_description)










