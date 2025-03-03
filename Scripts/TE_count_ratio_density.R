# Load required libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Read the data
data <- read.table("TE_ratios.txt", header = TRUE, sep = "\t")

# Reshape the data to long format
data_long <- melt(
  data,
  id.vars = "ID",                      # Specify the ID column
  variable.name = "Pairwise_Comparison", # Column for original column names
  value.name = "Value"                  # Column for values
)

# Set custom factor levels to control layering order
data_long$Pairwise_Comparison <- factor(
  data_long$Pairwise_Comparison,
  levels = c("PestPkai", "PantPkai", "PantPest") # Specify the desired order
)

# Bootstrapping function to calculate confidence intervals for the median
bootstrap_median <- function(values, n = 1000, conf = 0.95) {
  medians <- replicate(n, median(sample(values, replace = TRUE)))
  lower <- quantile(medians, (1 - conf) / 2)
  upper <- quantile(medians, 1 - (1 - conf) / 2)
  c(lower, upper)
}

# Calculate medians and confidence intervals for each pairwise comparison
results <- data_long %>%
  group_by(Pairwise_Comparison) %>%
  summarise(
    Median = median(Value),
    LowerCI = bootstrap_median(Value)[1],
    UpperCI = bootstrap_median(Value)[2],
    y_position = as.numeric(factor(Pairwise_Comparison, levels = c("PestPkai", "PantPkai", "PantPest"))) * -0.1 # Set y values for display
  )

# Print results
print("Medians with 95% Confidence Intervals:")
print(results)

# Create the density plot
p <- ggplot(data_long, aes(x = Value, fill = Pairwise_Comparison, color = Pairwise_Comparison)) +
  geom_density(alpha = 0.75) + # (alpha = 0.5)
  scale_x_log10(
    breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 100, 1000), # Specify custom tick marks
    labels = scales::comma                    # Format tick labels nicely
  ) +
#  geom_vline(data = medians, aes(xintercept = Value, color = Pairwise_Comparison), 
#             linetype = "dashed", size = 0.8) + # Add medians with matching colors
  # Add confidence intervals as brackets
  geom_segment(data = results, aes(
    x = LowerCI, xend = UpperCI, 
    y = y_position, yend = y_position, 
    color = Pairwise_Comparison
  ), size = 0.8) +
  # Add medians as points
  geom_point(data = results, aes(
    x = Median, y = y_position, 
    color = Pairwise_Comparison
  ), size = 2) +  
  
  theme_minimal() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1), # Add a border
    legend.position = c(0.95, 0.95),                                  # Position legend
    legend.justification = c(1, 1),                                    # Align legend top-right
    legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid')
  ) +
  labs(
    x = "Relative TE abundance",
    y = "Density"
  ) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9")) + # Set custom fill colors
#  scale_color_manual(values = c("#999999","#E69F00", "#56B4E9")) # Set custom line colors
  scale_color_manual(values = c("black", "black", "black"))
