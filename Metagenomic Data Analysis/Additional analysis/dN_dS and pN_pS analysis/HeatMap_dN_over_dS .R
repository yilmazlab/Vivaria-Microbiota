# Load necessary libraries
library(pheatmap)
library(viridis)  # For the viridis color palette
library(tibble)
library(dplyr)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(cluster)

Taxa <- read.table("data/data_for_heatmap_dN_dS.csv", header = TRUE, sep = ",")

# Read the metadata and data files
annotation_table <- read.table("data/metadata.csv", header = TRUE, sep = ",")

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$sample_id 

# Only show rows where at least half the values are not NA
num_columns <- ncol(Taxa)
Taxa <- Taxa[rowSums(!is.na(Taxa)) >= (num_columns / 2), ]

# Keep the first column as protein names
protein_names <- Taxa$product
Taxa <- Taxa %>%
  select(-"product")

# Convert the dataframe to a matrix
Taxa <- as.matrix(Taxa)

# Add protein names as row names
rownames(Taxa) <- protein_names

# Ensure annotation_table is ordered by the samples in Taxa
# This will make sure that sample_id and sample order matches
annotation_table_ordered <- annotation_table[match(colnames(Taxa), annotation_table$sample), ]

# Set the sample_id as column names of Taxa
colnames(Taxa) <- annotation_table_ordered$sample_id


# Define annotation colors
metadata_all <- read.table("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/Run1_5_metadata.txt", header = TRUE, sep = "\t")
location <- unique(metadata_all$Location)
location_colors <- viridis(length(location))
location_color_assignment <- setNames(location_colors, location)

unique_locations <- unique(annotation_table$Location)
filtered_location_color_assignment <- location_color_assignment[unique_locations]

ann_colors <- list(
  Gender = c(Male = "#498dcb", Female = "#f37c79", Unknown ="grey"),
  SamplingLocation = c(Feces = "#3d280a", Cecum = "#805212"),
  Group = c(SPF = "#f37c79", WildMice = "#498dcb"),
  Location = filtered_location_color_assignment,
  Age= c(Young = "#CCFFCC", Adult = "#99CCFF")
)

annotation_table <- annotation_table %>%
  select(-sample, -sample_id)

# Flatten the dataframe to a vector, ignoring NA values
Taxa_no_na <- unlist(Taxa)

# Calculate min and max for the whole dataframe
min_value <- min(Taxa_no_na, na.rm = TRUE)
max_value <- max(Taxa_no_na, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, ensuring we cover the range of data
myBreaks <- c(
  seq(min_value, 0.9999, length.out = ceiling(paletteLength / 2)),  # Values below 1
  1,  # Exact break for 1
  seq(1.0001, max_value, length.out = floor(paletteLength / 2))  # Values above 1
)

# Define the color palette for values below and above 1, with grey for exactly 1
my_palette <- c(
  colorRampPalette(c("#003399", "#99CCFF"))(paletteLength / 2 - 1),  # Values below 1 (minus 1 for the grey)
  "black",  # Exact color for the value 1
  colorRampPalette(c("#FFCC66", "#CC3300"))(paletteLength / 2)  # Values above 1
)

# Custom distance function to handle NA values
custom_distance <- function(x) {
  dist_matrix <- as.matrix(dist(x, method = "euclidean"))  # Compute Euclidean distance
  dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)  # Replace NA with the max distance
  dist_matrix
}

# Compute custom distances
dist_rows <- custom_distance(Taxa)
dist_cols <- custom_distance(t(Taxa))

# Perform hierarchical clustering
hclust_rows <- hclust(as.dist(dist_rows))
hclust_cols <- hclust(as.dist(dist_cols))

# Generate the heatmap with improved visualization settings
p <- pheatmap(Taxa, 
              fontsize_row = 8,
              fontsize_col = 11,           # Increase column font size
              na_col = "white",            # Color for NA values
              cluster_cols = hclust_cols,  # Enable column clustering
              cluster_rows = hclust_rows,  # Enable row clustering
              color = my_palette,          # Color palette
              breaks = myBreaks,           # Breaks for color scale
              cellwidth = 15,              # Adjust cell width
              cellheight = 8,              # Adjust cell height
              cutree_cols = 1,             # Number of clusters for columns
              annotation_col = annotation_table,  # Column annotations (if needed)
              annotation_colors = ann_colors,     # Colors for annotations
              border_color = "grey60")     # Border color for heatmap cells

p

# Save the plot
ggsave("Heatmaps/Bacteroides acidifaciens heatmap for dN over dS.pdf", plot = p, width = 80, height = 250, useDingbats = FALSE, limitsize = FALSE)


# Remove Hypothetical Proteins --------------------------------------------
# Filter out hypothetical proteins
# Filter out rows where the 'row_name_column' contains 'hypothetical protein' or no product information 

Taxa <- read.table("data/data_for_heatmap_dN_dS.csv", header = TRUE, sep = ",")

# Filter out rows where 'row_name_column' contains 'hypothetical protein'
Taxa <- Taxa %>%
  filter(!grepl('hypothetical protein', product, ignore.case = TRUE))

Taxa <- Taxa %>%
  filter(!grepl('No product information', product, ignore.case = TRUE))


# Read the metadata and data files
annotation_table <- read.table("data/metadata.csv", header = TRUE, sep = ",")

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$sample_id 

# Only show rows where at least half the values are not NA
num_columns <- ncol(Taxa)
Taxa <- Taxa[rowSums(!is.na(Taxa)) >= (num_columns / 2), ]

# Keep the first column as protein names
protein_names <- Taxa$product
Taxa <- Taxa %>%
  select(-"product")

# Convert the dataframe to a matrix
Taxa <- as.matrix(Taxa)

# Add protein names as row names
rownames(Taxa) <- protein_names

# Ensure annotation_table is ordered by the samples in Taxa
# This will make sure that sample_id and sample order matches
annotation_table_ordered <- annotation_table[match(colnames(Taxa), annotation_table$sample), ]

# Set the sample_id as column names of Taxa
colnames(Taxa) <- annotation_table_ordered$sample_id

# Define annotation colors
metadata_all <- read.table("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/Run1_5_metadata.txt", header = TRUE, sep = "\t")
location <- unique(metadata_all$Location)
location_colors <- viridis(length(location))
location_color_assignment <- setNames(location_colors, location)

unique_locations <- unique(annotation_table$Location)
filtered_location_color_assignment <- location_color_assignment[unique_locations]

ann_colors <- list(
  Gender = c(Male = "#498dcb", Female = "#f37c79", Unknown ="grey"),
  SamplingLocation = c(Feces = "#3d280a", Cecum = "#805212"),
  Group = c(SPF = "#f37c79", WildMice = "#498dcb"),
  Location = filtered_location_color_assignment,
  Age= c(Young = "#CCFFCC", Adult = "#99CCFF")
)

annotation_table <- annotation_table %>%
  select(-sample, -sample_id)

# Flatten the dataframe to a vector, ignoring NA values
Taxa_no_na <- unlist(Taxa)

# Calculate min and max for the whole dataframe
min_value <- min(Taxa_no_na, na.rm = TRUE)
max_value <- max(Taxa_no_na, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, ensuring we cover the range of data
myBreaks <- c(
  seq(min_value, 0.9999, length.out = ceiling(paletteLength / 2)),  # Values below 1
  1,  # Exact break for 1
  seq(1.0001, max_value, length.out = floor(paletteLength / 2))  # Values above 1
)

# Define the color palette for values below and above 1, with grey for exactly 1
my_palette <- c(
  colorRampPalette(c("#003399", "#99CCFF"))(paletteLength / 2 - 1),  # Values below 1 (minus 1 for the grey)
  "black",  # Exact color for the value 1
  colorRampPalette(c("#FFCC66", "#CC3300"))(paletteLength / 2)  # Values above 1
)

# Custom distance function to handle NA values
custom_distance <- function(x) {
  dist_matrix <- as.matrix(dist(x, method = "euclidean"))  # Compute Euclidean distance
  dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)  # Replace NA with the max distance
  dist_matrix
}

# Compute custom distances
dist_rows <- custom_distance(Taxa)
dist_cols <- custom_distance(t(Taxa))

# Perform hierarchical clustering
hclust_rows <- hclust(as.dist(dist_rows))
hclust_cols <- hclust(as.dist(dist_cols))

# Generate the heatmap with improved visualization settings
p <- pheatmap(Taxa, 
              fontsize_row = 8,
              fontsize_col = 11,           # Increase column font size
              na_col = "white",            # Color for NA values
              cluster_cols = hclust_cols,  # Enable column clustering
              cluster_rows = hclust_rows,  # Enable row clustering
              color = my_palette,          # Color palette
              breaks = myBreaks,           # Breaks for color scale
              cellwidth = 20,              # Adjust cell width
              cellheight = 8,              # Adjust cell height
              cutree_cols = 1,             # Number of clusters for columns
              annotation_col = annotation_table,  # Column annotations (if needed)
              annotation_colors = ann_colors,     # Colors for annotations
              border_color = "grey60")     # Border color for heatmap cells

p

# Save the plot
ggsave("Heatmaps/Bacteroides acidifaciens heatmap for dN over dS without hypothetical proteins.pdf", plot = p, width = 35, height = 150, useDingbats = FALSE, limitsize = FALSE)


# Remove Hypothetical Proteins + only show genes with positive selection--------------------------------------------
# Filter out hypothetical proteins
# Filter out rows where the 'row_name_column' contains 'hypothetical protein' or no product information 

Taxa <- read.table("data/data_for_heatmap_dN_dS.csv", header = TRUE, sep = ",")

# Filter out rows where 'row_name_column' contains 'hypothetical protein'
Taxa <- Taxa %>%
  filter(!grepl('hypothetical protein', product, ignore.case = TRUE))

Taxa <- Taxa %>%
  filter(!grepl('No product information', product, ignore.case = TRUE))


# Exclude the 'product' column
Taxa_no_product <- Taxa[, !(names(Taxa) %in% "product")]

# Filter rows where at least one value in the remaining columns is greater than 1, ignoring NA
filtered_rows <- apply(Taxa_no_product, 1, function(row) any(row > 1, na.rm = TRUE))

# Subset the original dataframe using the filtered rows
Taxa <- Taxa[filtered_rows, ]


# Read the metadata and data files
annotation_table <- read.table("data/metadata.csv", header = TRUE, sep = ",")

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$sample_id 

# Only show rows where at least half the values are not NA
num_columns <- ncol(Taxa)
Taxa <- Taxa[rowSums(!is.na(Taxa)) >= (num_columns / 2), ]

# Keep the first column as protein names
protein_names <- Taxa$product
Taxa <- Taxa %>%
  select(-"product")

# Convert the dataframe to a matrix
Taxa <- as.matrix(Taxa)

# Add protein names as row names
rownames(Taxa) <- protein_names

# Ensure annotation_table is ordered by the samples in Taxa
# This will make sure that sample_id and sample order matches
annotation_table_ordered <- annotation_table[match(colnames(Taxa), annotation_table$sample), ]

# Set the sample_id as column names of Taxa
colnames(Taxa) <- annotation_table_ordered$sample_id

# Define annotation colors
metadata_all <- read.table("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/Run1_5_metadata.txt", header = TRUE, sep = "\t")
location <- unique(metadata_all$Location)
location_colors <- viridis(length(location))
location_color_assignment <- setNames(location_colors, location)

unique_locations <- unique(annotation_table$Location)
filtered_location_color_assignment <- location_color_assignment[unique_locations]

ann_colors <- list(
  Gender = c(Male = "#498dcb", Female = "#f37c79", Unknown ="grey"),
  SamplingLocation = c(Feces = "#3d280a", Cecum = "#805212"),
  Group = c(SPF = "#f37c79", WildMice = "#498dcb"),
  Location = filtered_location_color_assignment,
  Age= c(Young = "#CCFFCC", Adult = "#99CCFF")
)

annotation_table <- annotation_table %>%
  select(-sample, -sample_id)

# Flatten the dataframe to a vector, ignoring NA values
Taxa_no_na <- unlist(Taxa)

# Calculate min and max for the whole dataframe
min_value <- min(Taxa_no_na, na.rm = TRUE)
max_value <- max(Taxa_no_na, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, ensuring we cover the range of data
myBreaks <- c(
  seq(min_value, 0.9999, length.out = ceiling(paletteLength / 2)),  # Values below 1
  1,  # Exact break for 1
  seq(1.0001, max_value, length.out = floor(paletteLength / 2))  # Values above 1
)

# Define the color palette for values below and above 1, with grey for exactly 1
my_palette <- c(
  colorRampPalette(c("#003399", "#99CCFF"))(paletteLength / 2 - 1),  # Values below 1 (minus 1 for the grey)
  "black",  # Exact color for the value 1
  colorRampPalette(c("#FFCC66", "#CC3300"))(paletteLength / 2)  # Values above 1
)

# Custom distance function to handle NA values
custom_distance <- function(x) {
  dist_matrix <- as.matrix(dist(x, method = "euclidean"))  # Compute Euclidean distance
  dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)  # Replace NA with the max distance
  dist_matrix
}

# Compute custom distances
dist_rows <- custom_distance(Taxa)
dist_cols <- custom_distance(t(Taxa))

# Perform hierarchical clustering
hclust_rows <- hclust(as.dist(dist_rows))
hclust_cols <- hclust(as.dist(dist_cols))

# Generate the heatmap with improved visualization settings
p <- pheatmap(Taxa, 
              fontsize_row = 8,
              fontsize_col = 11,           # Increase column font size
              na_col = "white",            # Color for NA values
              cluster_cols = hclust_cols,  # Enable column clustering
              cluster_rows = hclust_rows,  # Enable row clustering
              color = my_palette,          # Color palette
              breaks = myBreaks,           # Breaks for color scale
              cellwidth = 20,              # Adjust cell width
              cellheight = 8,              # Adjust cell height
              cutree_cols = 1,             # Number of clusters for columns
              annotation_col = annotation_table,  # Column annotations (if needed)
              annotation_colors = ann_colors,     # Colors for annotations
              border_color = "grey60")     # Border color for heatmap cells

p

# Save the plot
ggsave("Heatmaps/Bacteroides acidifaciens heatmap for dN over dS without hypothetical proteins only positively selected.pdf", plot = p, width = 40, height = 25, useDingbats = FALSE, limitsize = FALSE)



#Only show genes with positive selection--------------------------------------------

Taxa <- read.table("data/data_for_heatmap_dN_dS.csv", header = TRUE, sep = ",")

Taxa <- Taxa[rowSums(!is.na(Taxa)) >= 26, ]

# Exclude the 'product' column
Taxa_no_product <- Taxa[, !(names(Taxa) %in% "product")]

# Filter rows where at least one value in the remaining columns is greater than 1, ignoring NA
filtered_rows <- apply(Taxa_no_product, 1, function(row) any(row > 1, na.rm = TRUE))

# Subset the original dataframe using the filtered rows
Taxa <- Taxa[filtered_rows, ]

# Read the metadata and data files
annotation_table <- read.table("data/metadata.csv", header = TRUE, sep = ",")

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$sample_id 

# Only show rows where at least half the values are not NA
num_columns <- ncol(Taxa)
Taxa <- Taxa[rowSums(!is.na(Taxa)) >= (num_columns / 2), ]

# Keep the first column as protein names
protein_names <- Taxa$product
Taxa <- Taxa %>%
  select(-"product")

# Convert the dataframe to a matrix
Taxa <- as.matrix(Taxa)

# Add protein names as row names
rownames(Taxa) <- protein_names

# Ensure annotation_table is ordered by the samples in Taxa
# This will make sure that sample_id and sample order matches
annotation_table_ordered <- annotation_table[match(colnames(Taxa), annotation_table$sample), ]

# Set the sample_id as column names of Taxa
colnames(Taxa) <- annotation_table_ordered$sample_id

# Define annotation colors
metadata_all <- read.table("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/Run1_5_metadata.txt", header = TRUE, sep = "\t")
location <- unique(metadata_all$Location)
location_colors <- viridis(length(location))
location_color_assignment <- setNames(location_colors, location)

unique_locations <- unique(annotation_table$Location)
filtered_location_color_assignment <- location_color_assignment[unique_locations]

ann_colors <- list(
  Gender = c(Male = "#498dcb", Female = "#f37c79", Unknown ="grey"),
  SamplingLocation = c(Feces = "#3d280a", Cecum = "#805212"),
  Group = c(SPF = "#f37c79", WildMice = "#498dcb"),
  Location = filtered_location_color_assignment,
  Age= c(Young = "#CCFFCC", Adult = "#99CCFF")
)

annotation_table <- annotation_table %>%
  select(-sample, -sample_id)

# Flatten the dataframe to a vector, ignoring NA values
Taxa_no_na <- unlist(Taxa)

# Calculate min and max for the whole dataframe
min_value <- min(Taxa_no_na, na.rm = TRUE)
max_value <- max(Taxa_no_na, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, ensuring we cover the range of data
myBreaks <- c(
  seq(min_value, 0.9999, length.out = ceiling(paletteLength / 2)),  # Values below 1
  1,  # Exact break for 1
  seq(1.0001, max_value, length.out = floor(paletteLength / 2))  # Values above 1
)

# Define the color palette for values below and above 1, with grey for exactly 1
my_palette <- c(
  colorRampPalette(c("#003399", "#99CCFF"))(paletteLength / 2 - 1),  # Values below 1 (minus 1 for the grey)
  "black",  # Exact color for the value 1
  colorRampPalette(c("#FFCC66", "#CC3300"))(paletteLength / 2)  # Values above 1
)

# Custom distance function to handle NA values
custom_distance <- function(x) {
  dist_matrix <- as.matrix(dist(x, method = "euclidean"))  # Compute Euclidean distance
  dist_matrix[is.na(dist_matrix)] <- max(dist_matrix, na.rm = TRUE)  # Replace NA with the max distance
  dist_matrix
}

# Compute custom distances
dist_rows <- custom_distance(Taxa)
dist_cols <- custom_distance(t(Taxa))

# Perform hierarchical clustering
hclust_rows <- hclust(as.dist(dist_rows))
hclust_cols <- hclust(as.dist(dist_cols))

# Generate the heatmap with improved visualization settings
p <- pheatmap(Taxa, 
              fontsize_row = 8,
              fontsize_col = 11,           # Increase column font size
              na_col = "white",            # Color for NA values
              cluster_cols = hclust_cols,  # Enable column clustering
              cluster_rows = hclust_rows,  # Enable row clustering
              color = my_palette,          # Color palette
              breaks = myBreaks,           # Breaks for color scale
              cellwidth = 20,              # Adjust cell width
              cellheight = 8,              # Adjust cell height
              cutree_cols = 1,             # Number of clusters for columns
              annotation_col = annotation_table,  # Column annotations (if needed)
              annotation_colors = ann_colors,     # Colors for annotations
              border_color = "grey60")     # Border color for heatmap cells

p

# Save the plot
ggsave("Heatmaps/Bacteroides acidifaciens heatmap for dN over dS only positively selected.pdf", plot = p, width = 40, height = 65, useDingbats = FALSE, limitsize = FALSE)







