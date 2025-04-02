# Load required libraries
library(pheatmap)  # For creating heatmaps
library(RColorBrewer)  # For color palettes

# Set working directory
setwd("/Users/bahti/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Metagenomic/Phylophlan/2000")

# Load annotation data for columns and rows
annotation_table = read.table("metadata_anno.txt")  # Metadata for columns
annotation_row = read.table("annotation_row.txt")  # Metadata for rows

# Load taxa data
Taxa <- read.table("phylophlan_sgb_pres_abs.txt")  # Full taxa presence/absence data
Taxa <- read.table("phylophlan_sgb_pres_abs_subset.txt")  # Subset of taxa data

# Define annotation colors for heatmap
ann_colors = list(
  Gender = c(Female = "#cc6666", Male = "#B5BBE3", Unknown = "#f4fcfc"),  # Gender categories
  Content = c(Cecum = "#e8a188", Feces = "brown", Feces_Cecum = "#361010"),  # Sample content types
  Group = c(WildMice = "#38af60", SPF = "#3399FF", Human = "#efb94c"),  # Experimental groups
  SGB = c(kSGB = "black", uSGB = "yellow"),  # SGB categories
  Kingdom = c(Bacteria = "#41b6c4", Archaea = "orange"),  # Kingdom-level categories
  Phylum = c(  # Phylum-level categories
    Actinobacteria = "#E07B91", Bacteroidetes = "#8595E1", Firmicutes = "#38af60",
    Cyanobacteria = "cyan", Deferribacterota = "purple", Euryarchaeota = "orange",
    Fusobacteria = "#D33F6A", Lentisphaerota = "#023b10", Mycoplasmatota = "#aaf2ff",
    Proteobacteria = "yellow", Spirochaetota = "brown", Tenericutes = "#E6AFB9",
    Thermodesulfobacteria = "#abf5be", Unknown = "black", Verrucomicrobia = "#E07B91"
  )
)

# Define color palette for heatmap
# rwbcols <- c("#D33F6A","#E07B91","#E6AFB9","#f4fcfc","#B5BBE3", "#8595E1", "#4A6FE3")  # Alternative palette
rwbcols <- c("#f4fcfc", "#8595E1")  # Simplified palette

# Set the number of colors in the palette
paletteLength <- 1000  # Adjust as needed

# Generate the heatmap
a <- pheatmap(
  Taxa,  # Input data
  fontsize_row = 2, fontsize_col = 20,  # Font sizes for rows and columns
  cluster_cols = TRUE, cluster_rows = TRUE,  # Enable clustering for rows and columns
  color = rwbcols,  # Color palette
  cellwidth = 50, cellheight = 2,  # Cell dimensions
  # cutree_cols = 6, cutree_rows = 11,  # Uncomment to cut dendrograms into clusters
  annotation_col = annotation_table,  # Column annotations
  annotation_row = annotation_row,  # Row annotations
  annotation_colors = ann_colors,  # Annotation colors
  border_color = "grey60"  # Border color for cells
)

# Display the heatmap
a