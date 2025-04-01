
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readxl)
library(splines)

# Load the data
data <- read_excel("Associations.xlsx")

# Calculate the unique X_features per Y_feature
unique_x_features <- data %>%
  group_by(Y_features) %>%
  summarise(unique_x_features = n_distinct(X_features)) %>%
  arrange(unique_x_features)

# Add cumulative frequency
unique_x_features <- unique_x_features %>%
  mutate(cumulative_frequency = cumsum(unique_x_features) / sum(unique_x_features) * 100)

# Plot the data
ggplot(unique_x_features, aes(x = unique_x_features, y = cumulative_frequency)) +
  geom_line(stat = "smooth", method = "loess", se = TRUE, color = "blue") +
  geom_point(color = "red") +
  labs(
    title = "Smoothed Cumulative Frequency Distribution with 95% CI and Data Points",
    x = "Total Number of Associations (Unique X_features per Y_feature)",
    y = "Cumulative Frequency of X_features (%)"
  ) +
  theme_minimal()
