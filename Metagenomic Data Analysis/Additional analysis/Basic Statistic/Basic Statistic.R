# Load required libraries
library(readxl)      # For reading Excel files
library(ggplot2)     # For creating plots
library(dplyr)       # For data manipulation
library(ggstatsplot) # For statistical plots

# Set working directory
setwd("/Users/bahti/Dropbox/Ongoing Analysis/Wild-Type Mice vs Wild Mice Analysis/2022/Metagenomic/Phylophlan/")

# Step 1: Read the combined_data.xlsx file
data <- read_excel("combined_data.xlsx")

# Step 2: Filter data based on completeness and contamination thresholds
# Keep rows where completeness > 50 and contamination < 10
data_filtered <- data %>% filter(completeness > 50, contamination < 10)

# Step 3: Define color codes for the Groupv2 variable
# Assign specific colors to each group for consistent visualization
color_codes <- c("Human" = "darkgreen", "WildMice" = "darkblue", "SPF" = "darkred")

# Step 4: Create a box plot for completeness with jittered points
plot <- ggplot(data_filtered, aes(x = Groupv2, y = completeness, fill = Groupv2)) +
  geom_boxplot() +  # Add box plot
  geom_jitter(color = "black", width = 0.15, size = 1, alpha = 0.5) +  # Add jittered points
  labs(x = "Groupv2", y = "Completeness") +  # Label axes
  ggtitle("Completeness Box Plot") +  # Add title
  scale_fill_manual(values = color_codes) +  # Use defined color codes
  theme_bw()  # Apply a clean theme
plot

# Step 5: Create a violin plot for completeness with jittered points
plot <- ggplot(data_filtered, aes(x = Groupv2, y = completeness, fill = Groupv2)) +
  geom_boxplot() +  # Add box plot
  geom_jitter(aes(color = Location), size = 1, width = 0.3, alpha = 0.8) +  # Add jittered points
  labs(x = "Groupv2", y = "Completeness") +  # Label axes
  ggtitle("Completeness according to group") +  # Add title
  scale_fill_manual(values = color_codes) +  # Use defined color codes
  theme_bw() +  # Apply a clean theme
  facet_wrap(~ Groupv2, scales = "free")  # Create separate panels for each group
plot
ggsave("Completeness according to group.pdf", plot, width = 15, height = 9)  # Save the plot

# Step 6: Create a statistical plot for contamination
plt <- ggbetweenstats(
  data = data_filtered,
  x = Groupv2,
  y = contamination
)
plt
ggsave("Completeness according to group_violin plot.pdf", plt, width = 15, height = 9)  # Save the plot

# Step 7: Create a violin plot for contamination with jittered points
plot <- ggplot(data_filtered, aes(x = Groupv2, y = contamination, fill = Groupv2)) +
  geom_boxplot() +  # Add box plot
  geom_jitter(aes(color = Location), size = 1, width = 0.3, alpha = 0.8) +  # Add jittered points
  labs(x = "Groupv2", y = "Contamination") +  # Label axes
  ggtitle("Contamination according to group") +  # Add title
  scale_fill_manual(values = color_codes) +  # Use defined color codes
  theme_bw() +  # Apply a clean theme
  facet_wrap(~ Groupv2, scales = "free")  # Create separate panels for each group
plot
ggsave("Contamination according to group_violin plot.pdf", plt, width = 15, height = 9)  # Save the plot

# Step 8: Perform statistical analysis
# Calculate mean and standard deviation of completeness for each group
stats <- data_filtered %>%
  filter(Groupv2 %in% c("Human", "WildMice", "SPF")) %>%
  group_by(Groupv2) %>%
  summarise(mean_completeness = mean(completeness), sd_completeness = sd(completeness))

# Perform one-way ANOVA to test for differences in contamination between groups
anova_result <- aov(contamination ~ Groupv2, data = data_filtered)
anova_summary <- summary(anova_result)

# Extract p-value from ANOVA results
p_value <- anova_summary[[1]]["Pr(>F)"]

# Perform Tukey HSD test for pairwise comparisons
tukey_result <- TukeyHSD(anova_result)

# Extract Tukey comparison results
tukey_comparison <- data.frame(as.matrix(tukey_result$Groupv2))

# Add significance stars based on p-value thresholds
stars <- ifelse(p_value < 0.001, "***",
                ifelse(p_value < 0.01, "**",
                       ifelse(p_value < 0.05, "*", "")))

# Step 9: Annotate the box plot with Tukey comparison results
plot <- ggplot(data_filtered, aes(x = Groupv2, y = completeness, fill = Groupv2)) +
  geom_boxplot() +  # Add box plot
  geom_jitter(aes(color = Location), size = 1, width = 0.3, alpha = 0.8) +  # Add jittered points
  labs(x = "Groupv2", y = "Completeness") +  # Label axes
  ggtitle("Completeness according to group") +  # Add title
  scale_fill_manual(values = color_codes) +  # Use defined color codes
  theme_bw() +  # Apply a clean theme
  facet_wrap(~ Groupv2, scales = "free") +  # Create separate panels for each group
  geom_text(data = tukey_comparison, aes(x = comparison, y = y, label = stars), 
            size = 6, color = "black", nudge_y = 1)  # Add significance stars
plot
ggsave("boxplot_with_comparison.png", plot, width = 8, height = 6)  # Save the plot

# Step 10: Categorize data into quality levels based on completeness and contamination
data_filtered <- data_filtered %>% 
  mutate(Quality = ifelse(completeness >= 90 & completeness < 100 & contamination <= 5, "High-Quality with low contamination", 
                          ifelse(completeness >= 90 & completeness < 100 & contamination > 5 & contamination < 10, "High-Quality with high contamination",
                                 ifelse(completeness >= 50 & completeness < 90 & contamination <= 5, "Medium-Quality with low contamination",
                                        ifelse(completeness >= 50 & completeness < 90 & contamination > 5 & contamination < 10, "Medium-Quality with high contamination",
                                               ifelse(completeness == 100 & contamination == 0, "Completed with no contamination",
                                                      ifelse(completeness == 100 & contamination > 0.00001 & contamination <= 5, "Completed with low contamination",
                                                             ifelse(completeness == 100 & contamination > 5 & contamination < 10, "Completed with high contamination", "Other"))))))))

# Step 11: Define color and size codes for quality levels
color_codes <- c("Medium-Quality with low contamination" = "gray",
                 "Medium-Quality with high contamination" = "darkgray",
                 "High-Quality with low contamination" = "#0080FF",
                 "High-Quality with high contamination" = "#004C99",
                 "Completed with high contamination" = "yellow",
                 "Completed with low contamination" = "orange",
                 "Completed with no contamination" = "#472b02")

size_codes <- c("Medium-Quality with low contamination" = 2, 
                "Medium-Quality with high contamination" = 3, 
                "High-Quality with low contamination" = 4, 
                "High-Quality with high contamination" = 5, 
                "Completed with high contamination" = 6,
                "Completed with low contamination" = 7,
                "Completed with no contamination" = 8)

# Step 12: Create a scatter plot for contamination vs completeness
plot <- ggplot(data_filtered, aes(x = completeness, y = contamination, color = Quality)) +
  geom_point(alpha = 0.6, aes(size = Quality)) +  # Add points with varying sizes
  scale_color_manual(values = color_codes, drop = FALSE) +  # Use defined color codes
  scale_size_manual(values = size_codes) +  # Use defined size codes
  ggtitle("Contamination vs Completeness") +  # Add title
  theme_bw()  # Apply a clean theme
plot
ggsave("Contamination vs Completeness.pdf", plot, width = 20, height = 8)  # Save the plot

# Step 13: Count the number of each quality level
quality_counts <- data_filtered %>% 
  count(Quality)

# Print the counts
print(quality_counts)
