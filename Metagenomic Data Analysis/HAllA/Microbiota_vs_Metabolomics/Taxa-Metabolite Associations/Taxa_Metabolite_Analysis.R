
# Load required libraries
library(ggplot2)
library(readxl)
library(reshape2)
library(pheatmap)
library(openxlsx)

# Load input files
associations_file <- 'Associations.xlsx'
metabolites_file <- 'Metabolitesv2.xlsx'
tsv_file <- 'X_original.tsv'

# Load data
associations_df <- read_excel(associations_file, sheet = 'Sheet1')
metabolites_df <- read_excel(metabolites_file, sheet = 'Sheet1')
tsv_data <- read.csv(tsv_file, sep = '\t')

# Analysis 1: Generate Taxa vs Metabolite Associations Bar Plot
taxa_metabolite_counts <- aggregate(Metabolite ~ Taxa, data=associations_df, FUN=length)

# Plot bar graph
pdf("Taxa_Metabolite_Associations.pdf", width=15, height=8)
ggplot(taxa_metabolite_counts, aes(x=reorder(Taxa, Metabolite), y=Metabolite)) +
  geom_bar(stat="identity", fill="steelblue") +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  xlab("Individual Taxa (sorted by number of associations)") +
  ylab("Total Number of Metabolite Associations") +
  ggtitle("Taxa Stratified by Metabolite Associations") + theme_bw()
dev.off()

# Analysis 2: Dot-Line Plot
pdf("Taxa_Metabolite_DotLine_Plot.pdf", width=15, height=8)
ggplot(taxa_metabolite_counts, aes(x=reorder(Taxa, Metabolite), y=Metabolite)) +
  geom_line(aes(group=1), color="blue") +
  geom_point(color="blue") +
  theme(axis.text.x = element_text(angle=-45, vjust=0.5, hjust=1)) +
  xlab("Individual Taxa (sorted by number of associations)") +
  ylab("Total Number of Metabolite Associations") +
  ggtitle("Dot-Line Plot: Taxa Stratified by Metabolite Associations") + theme_bw()
dev.off()

library(ggplot2)

# Example plot with corrected x-axis label rotation
ggplot(taxa_metabolite_counts, aes(x=reorder(Taxa, Metabolite), y=Metabolite)) +
  geom_line(aes(group=1), color="blue") +
  geom_point(color="blue") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +  # Correct x-axis label rotation
  xlab("Individual Taxa (sorted by number of associations)") +
  ylab("Total Number of Metabolite Associations") +
  ggtitle("Dot-Line Plot: Taxa Stratified by Metabolite Associations") +
  theme_bw()


# Analysis 3: Metabolite-Taxa Associations
metabolite_taxa_counts <- aggregate(Taxa ~ Metabolite, data=associations_df, FUN=length)

# Dot-Line Plot for Metabolites
pdf("Metabolite_Taxa_DotLine_Plot.pdf", width=15, height=8)
ggplot(metabolite_taxa_counts, aes(x=reorder(Metabolite, Taxa), y=Taxa)) +
  geom_line(aes(group=1), color="red") +
  geom_point(color="red") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +  # Correct x-axis label rotation
  xlab("Individual Metabolites (sorted by number of taxa associations)") + theme_bw()
  ylab("Total Number of Taxa Associations") +
  ggtitle("Dot-Line Plot: Metabolites Stratified by Number of Taxa Associations")
dev.off()

# Save Metabolite-Taxa Associations Excel File
write.xlsx(associations_df[, c("Taxa", "Metabolite")], "Metabolite_Taxa_Distribution.xlsx")

# Analysis 4: Generate Heatmap of Taxa Presence Across Vivaria
limited_taxa_list <- unique(associations_df$Taxa)
filtered_taxa_abundance <- tsv_data[tsv_data$MicrobiotaName %in% limited_taxa_list,]

# Convert to presence/absence format
binary_taxa_abundance <- as.data.frame(lapply(filtered_taxa_abundance[,-1], function(x) as.integer(x > 0)))
row.names(binary_taxa_abundance) <- filtered_taxa_abundance$MicrobiotaName

# Generate heatmap
pdf("Limited_Taxa_Presence_Heatmap.pdf", width=20, height=12)
pheatmap(binary_taxa_abundance, color = colorRampPalette(c("white", "blue"))(50), cluster_rows=FALSE, cluster_cols=FALSE,
         main="Corrected Presence-Absence Heatmap of Limited Taxa Across Vivaria")
dev.off()

# Print completion message
cat("Analysis completed. Files generated:\n")
cat("1. Taxa_Metabolite_Associations.pdf\n")
cat("2. Taxa_Metabolite_DotLine_Plot.pdf\n")
cat("3. Metabolite_Taxa_DotLine_Plot.pdf\n")
cat("4. Metabolite_Taxa_Distribution.xlsx (Excel)\n")
cat("5. Limited_Taxa_Presence_Heatmap.pdf\n")
