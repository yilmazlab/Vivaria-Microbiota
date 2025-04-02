
# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

# Load the data
data <- read_excel("Associations.xlsx", sheet = "Sheet1")

# Determine association type
data <- data %>%
  mutate(association_type = ifelse(association > 0, "positive", "negative"))

# Count positive and negative associations for each Y_feature
association_counts <- data %>%
  group_by(Y_features, association_type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = association_type, values_from = count, values_fill = 0)

# Save the results to an Excel file
write_xlsx(association_counts, "Corrected_Association_Counts.xlsx")
