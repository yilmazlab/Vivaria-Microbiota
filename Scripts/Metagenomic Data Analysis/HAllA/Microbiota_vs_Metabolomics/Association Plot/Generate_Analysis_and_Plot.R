
# Load necessary libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

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

# Create a stacked bar plot
association_counts <- association_counts %>%
  mutate(Y_features = factor(Y_features, levels = unique(Y_features)))

ggplot(association_counts, aes(x = Y_features)) +
  geom_bar(aes(y = positive, fill = "Positive"), stat = "identity") +
  geom_bar(aes(y = -negative, fill = "Negative"), stat = "identity") +
  scale_fill_manual(values = c("Positive" = "red", "Negative" = "green")) +
  labs(x = "Y_features", y = "Count", title = "Stacked Bar Plot of Positive and Negative Associations") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Save the plot as a PDF
ggsave("Positive_Negative_Associations.pdf", width = 18, height = 10)
