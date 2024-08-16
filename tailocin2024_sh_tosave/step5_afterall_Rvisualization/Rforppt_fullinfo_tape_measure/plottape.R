# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Read the sample names from treetopologyorder.txt
tree_topology_order <- readLines("treetopologyorder.txt")
tree_order_data <- data.frame(sample = tree_topology_order, stringsAsFactors = FALSE)

# Read the new file with haplotype names and lengths
haplotype_lengths <- read.table("alltapemeasureforfishlength.txt", header = FALSE, stringsAsFactors = FALSE, col.names = c("haplotype", "length"))

# Extract sample names from haplotype names
haplotype_lengths$sample <- sub("_tape$", "", haplotype_lengths$haplotype)

# Merge tree order data with haplotype lengths
combined <- merge(tree_order_data, haplotype_lengths, by.x = "sample", by.y = "sample", all.x = TRUE)

# Ensure all samples from file 1 are included and in the correct order
combined$sample <- factor(combined$sample, levels = rev(unique(tree_order_data$sample)))

# Set length to 0 for samples without corresponding lengths in the haplotype file
combined$length[is.na(combined$length)] <- 0

# Expand the data to plot each position
expanded_data <- combined %>%
  group_by(sample) %>%
  mutate(position = ifelse(length > 0, list(0:(length - 1)), list(0))) %>%
  unnest(cols = c(position))

# Define colors for different lengths
unique_lengths <- unique(expanded_data$length)
num_colors <- length(unique_lengths)
palette <- brewer.pal(min(num_colors, 12), "Set3") # Choose a perceptually uniform palette
length_colors <- setNames(palette, unique_lengths)
length_colors <- c("0" = "white", length_colors)

# Create the plot
p <- ggplot(expanded_data, aes(x = position, y = sample, color = as.factor(length))) +
  geom_point(position = position_jitter(height = 0.3)) +
  scale_color_manual(values = length_colors, na.translate = FALSE) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "",
       x = "",
       y = "",
       color = "Length") 
# Save the plot with proportional shrink on x-axis
ggsave("tapehaplotype_lengths_plot.pdf", plot = p, width = 5, height = 10, units = "in")