# Load necessary libraries
library(pheatmap)
library(viridis)  # For the viridis color palette
library(tibble)
library(dplyr)
library(readr)
library(RColorBrewer)
library(ggplot2)
library(cluster)
library(qs)

# Heatmap for nonsynonymous substitutions

hk_data <- qread(("data/hk_data.qs"))

presence_df <- qread(("data/detected_mutations.qs"))

# Keep the first column as mutation names
mutation <-presence_df$mutation_def
presence_df <- presence_df %>%
  select(-"mutation_def")


samples_for_heatmap <- unique(hk_data$new_sample_label)

annotation_table <- qread(("data/metadata.qs"))

annotation_table <- annotation_table %>%
  filter(new_sample_label %in% samples_for_heatmap)

annotation_ids <-  unique(annotation_table$new_sample_label)

# Keep only the columns in filtered_annotation_table that are in annotation_ids
presence_df <- presence_df %>%
  select(one_of(annotation_ids))

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$new_sample_label

# Convert the dataframe to a matrix
presence_df <- as.matrix(presence_df)

# Add protein names as row names
rownames(presence_df) <- mutation

# Set the sample_id as column names of Taxa
colnames(presence_df) <- annotation_table$new_sample_label

# Define annotation colors
metadata_all <-  qread(("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/metadata_all.qs"))
metadata_all <- metadata_all %>%
  mutate(Location = paste(country_abbreviation, institute_abbreviation, sep = "_"))

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
  select(-sample, -sample_id,-new_sample_label)


# Calculate min and max for the whole dataframe
min_value <- min(presence_df, na.rm = TRUE)
max_value <- max(presence_df, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, covering the range of data
myBreaks <- seq(min_value, 1, length.out = paletteLength)

# Define the color palette for values from 0 to 1, with blue for 0 and red for 1
my_palette <- colorRampPalette(c("#000099", "#CC0000"))(paletteLength)

library(pheatmap)
# Generate the heatmap with improved visualization settings
p <- pheatmap(presence_df, 
              fontsize_row = 8,  # Increase row font size
              fontsize_col = 15, # Increase column font size
              cluster_cols = TRUE,  # Enable column clustering for better separation
              cluster_rows = TRUE,
              color = my_palette,
              breaks = myBreaks,
              cellwidth = 20,  # Adjust cell width
              cellheight = 8, # Adjust cell height
              annotation_col = annotation_table,
              annotation_colors = ann_colors, 
              border_color = "grey60")
p
str(p)
# Save the plot
ggsave("Heatmaps/Amidophosphoribosyltransferase - distribution of nonsynonymous substitutions.pdf", plot = p, width = 25, height = 25, useDingbats = FALSE, limitsize = FALSE)

# Heatmap for synonymous substitutions

hk_data <- qread(("data/hk_data_synonymous.qs"))

presence_df <- qread(("data/detected_mutations_synonymous.qs"))

# Keep the first column as mutation names
mutation <-presence_df$mutation_def
presence_df <- presence_df %>%
  select(-"mutation_def")


samples_for_heatmap <- unique(hk_data$new_sample_label)

annotation_table <- qread(("data/metadata.qs"))

annotation_table <- annotation_table %>%
  filter(new_sample_label %in% samples_for_heatmap)

annotation_ids <-  unique(annotation_table$new_sample_label)

# Keep only the columns in filtered_annotation_table that are in annotation_ids
presence_df <- presence_df %>%
  select(one_of(annotation_ids))

# Set the first column as rownames
rownames(annotation_table) <- annotation_table$new_sample_label

# Convert the dataframe to a matrix
presence_df <- as.matrix(presence_df)

# Add protein names as row names
rownames(presence_df) <- mutation

# Set the sample_id as column names of Taxa
colnames(presence_df) <- annotation_table$new_sample_label

# Define annotation colors
metadata_all <-  qread(("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/metadata_all.qs"))
metadata_all <- metadata_all %>%
  mutate(Location = paste(country_abbreviation, institute_abbreviation, sep = "_"))

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
  select(-sample, -sample_id,-new_sample_label)


# Calculate min and max for the whole dataframe
min_value <- min(presence_df, na.rm = TRUE)
max_value <- max(presence_df, na.rm = TRUE)

# Define palette length
paletteLength <- 100  # Example length, adjust as needed

# Create breaks, covering the range of data
myBreaks <- seq(min_value, 1, length.out = paletteLength)

# Define the color palette for values from 0 to 1, with blue for 0 and red for 1
my_palette <- colorRampPalette(c("#000099", "#CC0000"))(paletteLength)

library(pheatmap)
# Generate the heatmap with improved visualization settings
p <- pheatmap(presence_df, 
              fontsize_row = 8,  # Increase row font size
              fontsize_col = 15, # Increase column font size
              cluster_cols = TRUE,  # Enable column clustering for better separation
              cluster_rows = TRUE,
              color = my_palette,
              breaks = myBreaks,
              cellwidth = 20,  # Adjust cell width
              cellheight = 8, # Adjust cell height
              annotation_col = annotation_table,
              annotation_colors = ann_colors, 
              border_color = "grey60")
p

# Save the plot
ggsave("Heatmaps/Amidophosphoribosyltransferase - distribution of synonymous substitutions.pdf", plot = p, width = 25, height = 25, useDingbats = FALSE, limitsize = FALSE)
