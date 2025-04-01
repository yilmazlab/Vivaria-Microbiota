
# Load necessary libraries
library(readxl)
library(tidyverse)

# Load the data
metabolite_list <- read_excel("MetaboliteList.xlsx", col_names = FALSE)$...1
associations_data <- read_excel("Associations.xlsx")
x_original_data <- read_tsv("X_original.tsv")

# Initialize an empty list to store concatenated abundance data
concatenated_abundance_data <- list()

# Iterate over each metabolite in the list
for (metabolite in metabolite_list) {
    # Filter for associated bacteria
    associated_bacteria <- associations_data %>% filter(Metabolite == metabolite) %>% pull(Taxa)
    
    # Concatenate abundance values for each bacterium
    abundance_values <- associated_bacteria %>%
        map(function(bacterium) {
            if (bacterium %in% x_original_data$MicrobiotaName) {
                x_original_data %>% filter(MicrobiotaName == bacterium) %>% select(-MicrobiotaName)
            } else {
                NULL
            }
        }) %>%
        compact() %>%
        bind_rows()

    # Store the concatenated values
    concatenated_abundance_data[[metabolite]] <- abundance_values %>% unlist() %>% na.omit()
}

# Convert to data frame and write to file
concatenated_abundance_df <- as.data.frame(concatenated_abundance_data)
write.xlsx(concatenated_abundance_df, "ConcatenatedRelativeAbundanceData.xlsx")

# Generate dot plot
pdf("ConcatenatedRelativeAbundanceDotPlot.pdf", width = 20, height = 10)
plot(1, type = "n", xlim = c(1, length(concatenated_abundance_df)), ylim = c(-10, 10), xaxt = "n",
     ylab = "Log10(Relative Abundance)", main = "Dot Plot of Concatenated Relative Abundance")
axis(1, at = 1:length(concatenated_abundance_df), labels = names(concatenated_abundance_df), las = 2)
for (i in 1:ncol(concatenated_abundance_df)) {
    points(rep(i, nrow(concatenated_abundance_df)), log10(concatenated_abundance_df[, i] + 1e-10), col = "blue", pch = 16)
}
dev.off()
