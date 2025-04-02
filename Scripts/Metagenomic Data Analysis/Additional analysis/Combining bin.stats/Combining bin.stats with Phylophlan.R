# Load the data.table library for efficient data manipulation
library(data.table)

# Define the path to the BIN_REASSEMBLY folder
folder_path <- "~/Desktop/Phylophlan/BIN_REASSEMBLY/"

# Get a list of all BIN_REASSEMBLY folders for different categories using specific patterns
folders_bahti <- list.files(folder_path, pattern = "BIN_REASSEMBLY_Bahti_H.*", full.names = TRUE)
folders_human <- list.files(folder_path, pattern = "BIN_REASSEMBLY_Human.*", full.names = TRUE)
folders_Crick <- list.files(folder_path, pattern = "BIN_REASSEMBLY_Crick.*", full.names = TRUE)
folders_Greece <- list.files(folder_path, pattern = "BIN_REASSEMBLY_Greece.*", full.names = TRUE)
folders_Poland <- list.files(folder_path, pattern = "BIN_REASSEMBLY_Poland.*", full.names = TRUE)

# Combine all folder paths into a single vector
folders <- c(folders_bahti, folders_human, folders_Crick, folders_Greece, folders_Poland)

# Initialize an empty data.table to store the combined data from all folders
combined_data <- data.table()

# Iterate over each folder to process the data
for (folder in folders) {
  # Construct the path to the reassembled_bins.stats file and read it
  stats_file <- file.path(folder, "reassembled_bins.stats")
  stats_data <- fread(stats_file)  # Read the stats file into a data.table
  
  # Construct the path to the phylophlan.tsv file and read it
  phylophlan_file <- file.path(folder, "reassembled_bins", "phylophlan.tsv")
  phylophlan_data <- fread(phylophlan_file)  # Read the PhyloPhlAn file into a data.table
  
  # Rename columns in stats_data for clarity and consistency
  setnames(stats_data, c("bin", "completeness", "contamination", "GC", "lineage", "N50", "size"))
  
  # Rename columns in phylophlan_data for clarity and consistency
  setnames(phylophlan_data, c("bin", "phylophlan_column"))
  
  # Merge the stats_data and phylophlan_data tables based on the "bin" column
  merged_data <- merge(stats_data, phylophlan_data, by = "bin")
  
  # Add a new column to indicate the origin folder of the data
  merged_data[, BIN_REASSEMBLY_origin := folder]
  
  # Append the merged data to the combined_data data.table
  combined_data <- rbind(combined_data, merged_data)
}

# Define the output file path for the combined data
output_file <- "~/Desktop/Phylophlan/combined_data.csv"

# Write the combined data to a CSV file without row names
write.csv(combined_data, output_file, row.names = FALSE)