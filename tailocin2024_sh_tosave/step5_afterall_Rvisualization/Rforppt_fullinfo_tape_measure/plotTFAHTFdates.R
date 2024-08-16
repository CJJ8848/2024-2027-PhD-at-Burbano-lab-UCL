# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# Read the sample names from treetopologyorder.txt
tree_topology_order <- readLines("treetopologyorder.txt")
tree_order_data <- data.frame(sample = tree_topology_order, stringsAsFactors = FALSE)
# Read the formatted HTF data
formatted_HTF <- readLines("formatted_HTF.txt")
# Read the HTF order data
HTF_order <- readLines("allHTFnames.txt")

# Process the HTF order data
order_data <- data.frame(sample = character(), haplotype = character(), stringsAsFactors = FALSE)
for (line in HTF_order) {
  elements <- unlist(strsplit(sub(">", "", line), "\\|"))
  if (length(elements) == 2) {
    order_data <- rbind(order_data, data.frame(sample = elements[1], haplotype = elements[2], stringsAsFactors = FALSE))
  } else if (length(elements) == 1) {
    order_data <- rbind(order_data, data.frame(sample = elements[1], haplotype = NA, stringsAsFactors = FALSE))
  }
}

# Merge tree order data with HTF order data
merged_data <- merge(tree_order_data, order_data, by = "sample", all.x = TRUE)



# Process the formatted HTF data
HTF_lengths <- data.frame(haplotype = character(), start = numeric(), end = numeric(), stringsAsFactors = FALSE)
for (line in formatted_HTF) {
  parts <- unlist(strsplit(line, ":"))
  haplotype <- parts[1]
  range <- unlist(strsplit(parts[2], "-"))
  start <- as.numeric(range[1])
  end <- as.numeric(range[2])
  HTF_lengths <- rbind(HTF_lengths, data.frame(haplotype = haplotype, start = start, end = end, stringsAsFactors = FALSE))
}

HTF_lengths$length <- HTF_lengths$end - HTF_lengths$start + 1

# Merge the order data with haplotype lengths
combined <- merge(merged_data, HTF_lengths, by.x = "haplotype", by.y = "haplotype", all.x = TRUE)
combined$length[is.na(combined$length)] <- 0

# Ensure all samples from file 1 are included and in the correct order
combined$sample <- factor(combined$sample, levels = rev(unique(tree_order_data$sample)))

# Expand the data to plot each position
expanded_data <- combined %>%
  group_by(sample) %>%
  mutate(position = ifelse(length > 0, list(0:(length - 1)), list(0))) %>%
  unnest(cols = c(position))

# Ensure correct handling of NA values and case consistency
expanded_data$haplotype <- as.character(expanded_data$haplotype)
expanded_data$haplotype[is.na(expanded_data$haplotype)] <- "NA"

# Define colors for haplotypes, white for NA
unique_haplotypes <- unique(expanded_data$haplotype)
num_colors <- length(unique_haplotypes) - 1
palette <- brewer.pal(min(num_colors, 12), "Set3") # Choose a perceptually uniform palette
haplotype_colors <- setNames(palette, unique_haplotypes[unique_haplotypes != "NA"])
haplotype_colors <- c("NA" = "white", haplotype_colors)

# Create the plot
p <- ggplot(expanded_data, aes(x = position, y = sample, color = haplotype)) +
  geom_point(position = position_jitter(height = 0.3)) +
  scale_color_manual(values = haplotype_colors, na.translate = FALSE) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "",
       x = "",
       y = "",
       color = "")

# Save the plot with proportional shrink on x-axis
ggsave("HTF_lengths_plot.pdf", plot = p, width = 5, height = 10, units = "in")














# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)

# Read the sample names from treetopologyorder.txt
tree_topology_order <- readLines("treetopologyorder.txt")
tree_order_data <- data.frame(sample = tree_topology_order, stringsAsFactors = FALSE)

# Read the sample dates from sampleanddate68.txt
sample_dates <- read.table("sampleanddate68.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(sample_dates) <- c("sample", "date")

# Merge the tree order data with sample dates
merged_data <- merge(tree_order_data, sample_dates, by = "sample", all.x = TRUE)

# Ensure all samples from file 1 are included and in the correct order
merged_data$sample <- factor(merged_data$sample, levels = rev(unique(tree_order_data$sample)))

# Process the dates to create lengths
merged_data <- merged_data %>%
  mutate(length = date ) %>%
  group_by(sample) %>%
  mutate(position = ifelse(!is.na(length), list(0:(length - 1)), list(0))) %>%
  unnest(cols = c(position))

# Plotting
p <- ggplot(merged_data, aes(x = position, y = sample)) +
  geom_point(color='grey') +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Sample Dates as Lengths",
       x = "Position",
       y = "Sample")

# Save the plot with proportional shrink on x-axis
ggsave("sample_dates_plot.png", plot = p, width = 5, height = 10, units = "in")