# Load necessary libraries
library(dplyr)
library(tidyr)
library(data.table)
library(readxl)
library(plotly)
library(ggplot2)
library(ggforce)
library(circlize)


# Select genome annotation
gene_annotation <- read_excel("data/gene_annotation_features.xlsx")

names(gene_annotation)[names(gene_annotation) == "gene_id"] <- "gene"

gene_annotation  <- gene_annotation %>%
  select(gene, scaffold, start_position, end_position)

# Select SNS
sns <- read.table("data/filtered_combined_data_coverage_10_final.csv", header = TRUE, sep = ",")

annotation_table <- read.table("/Users/isabel/Library/CloudStorage/OneDrive-UniversitaetBern/Wild mice/Gene variants analysis updated/example_gene_table/metadata_all.txt", header = TRUE, sep = "\t")

names(annotation_table)[names(annotation_table) == "sample"] <- "Sample"

annotation_table <- annotation_table %>%
  mutate(Location = paste(country_abbreviation, institute_abbreviation, sep = "_"))


annotation_table <- annotation_table %>%
  select(Sample,Location)

gene_annotation$scaffold <- as.numeric(sub(".*_", "", gene_annotation$scaffold))

gene_annotation$scaffold <- as.factor(gene_annotation$scaffold)

sns <- left_join(sns,annotation_table,by="Sample")

sns <- left_join(sns,gene_annotation,by="gene")

sns <- sns %>% filter(!is.na(start_position))

sns$scaffold.x <- as.factor(sns$scaffold.x)

# Extract the numeric part after the last underscore and convert to numeric
sns$scaffold_num <- as.numeric(sub(".*_", "", sns$scaffold.x))

sns$scaffold_num <- as.factor(sns$scaffold_num)

all_locations <- unique(sns$Location)


### Make a seperate plot for all locations
# Jena

sns_got <- subset(sns, Location == c("DE_FSU"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Jena_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Jena_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# SL

sns_got <- subset(sns, Location == c("DE_CAU_Schoemberg"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/SL_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/SL_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# CB

sns_got <- subset(sns, Location == c("DE_CAU_Cologne"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/CB_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/CB_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# RiodeJaneiro

sns_got <- subset(sns, Location == c("BR_UFRJ"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/RiodeJaneiro_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/RiodeJaneiro_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})
# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()




# FR_CAU_Angers

sns_got <- subset(sns, Location == c("FR_CAU_Angers"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/AN_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/AN_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()



# CH_UNIBE

sns_got <- subset(sns, Location == c("CH_UNIBE"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Bern_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Bern_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# FR_CAU_DivonneLesBains

sns_got <- subset(sns, Location == c("FR_CAU_DivonneLesBains"))

print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/DB_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/DB_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()

### Make a seperate plot for all locations
# Espelette

sns_got <- subset(sns, Location == c("FR_CAU_Espelette"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/ES_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Es_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()



### Make a seperate plot for all locations
# Freiburg

sns_got <- subset(sns, Location == c("DE_ALUF"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Freiburg_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Freiburg_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()

# Geneva

sns_got <- subset(sns, Location == c("CH_UNG"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Geneva_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Geneva_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()



# Greece 

sns_got <- subset(sns, Location == c("GR_BSRC"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Greece_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Greece_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()



# Humanitas 

sns_got <- subset(sns, Location == c("IT_HUNIMED"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Humanitas_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Humanitas_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# Istanbul

sns_got <- subset(sns, Location == c("TR_IMU"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Istanbul_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Istanbul_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# Louan

sns_got <- subset(sns, Location == c("FR_CAU_Louan"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Louan_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Louan_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()

# Massif central

sns_got <- subset(sns, Location == c("FR_CAU_MassifCentral"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/MC_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/MC_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()




# Rockefeller

sns_got <- subset(sns, Location == c("US_RU"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Rockefeller_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Rockefeller_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()


# Zürich

sns_got <- subset(sns, Location == c("CH_UZH"))


print(unique(sns_got$Sample))
unique_sample_count <- sns_got %>%
  summarise(unique_samples = n_distinct(Sample)) %>%
  pull(unique_samples)

location  <-  unique(as.character(sns_got$Location))[1]

sns_got <- sns_got %>%
  select(scaffold_num, position, mutation_type, class)

names(sns_got)[names(sns_got) == "scaffold_num"] <- "scaffold"

# prepare sns.csv

sns_got_SNS <- subset(sns_got, class == c("SNS"))

sns_got_SNS   <- sns_got_SNS   %>%
  select(scaffold,class,position,mutation_type)

names(sns_got_SNS)[names(sns_got_SNS) == "class"] <- "sns"

sns_got_SNS_S <- subset(sns_got_SNS, mutation_type == c("S"))
sns_got_SNS_N <- subset(sns_got_SNS, mutation_type == c("N"))

# prepare con_snv.csv

con_snv_got_SNS <- subset(sns_got, class == c("con_SNV"))

con_snv_got_SNS  <- con_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(con_snv_got_SNS)[names(con_snv_got_SNS) == "class"] <- "con_snv"
con_snv_got_SNS_S <- subset(con_snv_got_SNS, mutation_type == c("S"))
con_snv_got_SNS_N <- subset(con_snv_got_SNS, mutation_type == c("N"))

# prepare pop_snv.csv

pop_snv_got_SNS <- subset(sns_got, class == c("pop_SNV"))

pop_snv_got_SNS  <- pop_snv_got_SNS  %>%
  select(scaffold,class,position,mutation_type)

names(pop_snv_got_SNS)[names(pop_snv_got_SNS) == "class"] <- "pop_snv"
pop_snv_got_SNS_S <- subset(pop_snv_got_SNS, mutation_type == c("S"))
pop_snv_got_SNS_N <- subset(pop_snv_got_SNS, mutation_type == c("N"))

# create plot
library(circlize)

# Assume 'contig' column exists in both gene_annotation and sns_got data
scaffold <- unique(gene_annotation$scaffold)

# Define the total length of each contig
scaffold_lengths <- aggregate(gene_annotation$end_position, by = list(gene_annotation$scaffold), FUN = max)
names(scaffold_lengths) <- c("scaffold", "length")


##Plot the Substitutions
# Define a helper function to add jitter
jitter_y <- function(n, spread = 0.5) {
  runif(n, min = 0.5 - spread, max = 0.5 + spread)
}

# Open a PDF device with specified dimensions
pdf("Circos plots/Zürich_SNP_profile_multi_contig_SNS.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_S data for the current scaffold/sector
  sns_S <- subset(sns_got_SNS_S, scaffold == sector.index)
  
  if (nrow(sns_S) > 0) {
    circos.points(
      sns_S$position, jitter_y(nrow(sns_S), spread = 0.5), 
      col = "#0099CC",  # Color for SNS_S SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Track for SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the SNS_N data for the current scaffold/sector
  sns_N <- subset(sns_got_SNS_N, scaffold == sector.index)
  
  if (nrow(sns_N) > 0) {
    circos.points(
      sns_N$position, jitter_y(nrow(sns_N), spread = 0.5), 
      col = "#FF9900",  # Color for SNS_N SNVs
      pch = 19,         # Filled circle
      cex = 0.05        # Small point size
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 4.1,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "SNS-synonymous" = "#0099CC",
  "SNS-nonsynonymous" = "#FF9900"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)


dev.off()


##Plot the Variants 

# Open a PDF device with specified dimensions
pdf("Circos plots/Zürich_SNP_profile_multi_contig_SNV.pdf", width = 10, height = 10)

# Set up padding to reduce space between sectors
circos.par(cell.padding = c(0.001, 0.001, 0.001, 0.001))  # Adjust the padding in x-direction

# Initialize the circular plot with one sector per scaffold
circos.initialize(factors = gene_annotation$scaffold, xlim = cbind(rep(0, nrow(scaffold_lengths)), scaffold_lengths$length))

# Plot the genes in the first track
circos.trackPlotRegion(factors = gene_annotation$scaffold, ylim = c(0, 1),track.height = 0.15, panel.fun = function(x, y) {
  current_scaffold <- get.cell.meta.data("sector.index")
  
  # Get gene data for the current scaffold
  contig_genes <- subset(gene_annotation, scaffold == current_scaffold)
  
  # Loop over the genes and plot them as rectangles in the current sector
  for (i in 1:nrow(contig_genes)) {
    circos.rect(
      contig_genes$start_position[i], 0, contig_genes$end_position[i], 1, 
      col = "lightblue", border = "lightgrey", lwd = 0.0000005
    )
  }
  
})


# Track for con_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_S data for the current scaffold/sector
  con_sns_S <- subset(con_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(con_sns_S) > 0) {
    circos.points(
      con_sns_S$position, jitter_y(nrow(con_sns_S), spread = 0.5), 
      col = "#0066CC",   # Color for con SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for con_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.08, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the con_snv_got_SNS_N data for the current scaffold/sector
  con_sns_N <- subset(con_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(con_sns_N) > 0) {
    circos.points(
      con_sns_N$position, jitter_y(nrow(con_sns_N), spread = 0.5), 
      col = "#CC0000",   # Color for con SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_S SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_S data for the current scaffold/sector
  pop_sns_S <- subset(pop_snv_got_SNS_S, scaffold == sector.index)
  
  if (nrow(pop_sns_S) > 0) {
    circos.points(
      pop_sns_S$position, jitter_y(nrow(pop_sns_S), spread = 0.5), 
      col = "#000066",   # Color for population SNS_S SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Track for pop_snv_got_SNS_N SNVs with jitter
circos.track(ylim = c(0, 1), track.height = 0.1, bg.col = "white", panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  
  # Filter the pop_snv_got_SNS_N data for the current scaffold/sector
  pop_sns_N <- subset(pop_snv_got_SNS_N, scaffold == sector.index)
  
  if (nrow(pop_sns_N) > 0) {
    circos.points(
      pop_sns_N$position, jitter_y(nrow(pop_sns_N), spread = 0.5), 
      col = "#990000",   # Color for population SNS_N SNVs
      pch = 19,          # Filled circle
      cex = 0.1          # Larger point size for visibility
    )
  }
})

# Iterate through the scaffold lengths to add labels
for (i in seq_len(nrow(scaffold_lengths))) {
  circos.text(
    x = scaffold_lengths$length[i] / 2,  # Center the label in the sector
    y = 6.2,                             # Position it outside the plot
    labels = scaffold_lengths$scaffold[i],  # Use scaffold name from scaffold_lengths
    sector.index = scaffold_lengths$scaffold[i],  # Associate with the correct sector
    facing = "bending.outside",          # Orient text outward for readability
    niceFacing = TRUE,                   # Smooth orientation
    adj = c(0.5, 0.5),                   # Center alignment
    cex = 0.25                            # Text size
  )
}
# Ensure the ring is closed and finalize the plot
circos.clear()

# Ensure the ring is closed and finalize the plot
circos.clear()


library(ComplexHeatmap)

# Define colors for SNV types
color_list <- list(
  "con_SNV-synonymous" = "#0066CC",
  "con_SNV-nonsynonymous" = "#CC0000",
  "pop_SNV-synonymous" = "#000066",
  "pop_SNV-nonsynonymous" = "#990000"
)

# Create the legend with a multi-line title
lgd_list_vertical <- Legend(
  at = names(color_list),  # Names of the color list are used as labels
  title = "SNP Type",  # Title with newline character
  direction = "vertical",
  legend_gp = gpar(fill = unlist(color_list))  # Use 'fill' in gpar to specify colors
)

# Draw the legend in the plot
draw(lgd_list_vertical, x = unit(5, "mm"), y = unit(5, "mm"), just = c("left", "bottom"))

# Assuming location and unique_sample_count are already defined
location_text <- paste("Location:", location)
sample_count_text <- paste("n:", unique_sample_count)

library(grid)
# Add a second legend (center text for location and sample count)
grid.text(
  paste(location_text, sample_count_text, sep = "\n"),  # Combine location and sample count in one label
  x = unit(0.5, "npc"),  # Position in the center (x-axis)
  y = unit(0.5, "npc"),  # Position in the center (y-axis)
  just = c("center", "center"),  # Centered text
  gp = gpar(fontsize = 12, fontface = "bold")  # Set text style
)

dev.off()

