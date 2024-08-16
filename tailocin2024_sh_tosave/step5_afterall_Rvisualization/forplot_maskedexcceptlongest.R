# Load necessary libraries
library(ggplot2)
library(dplyr)
library(readr)

# Set working directory
setwd('/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/tailocin46')

# Define file paths
merged_coverage_file <- "all_87fa_samples_tailocin_nonN_markexceptlongest.txt"
cumulative_positions_file <- "cumulative_positions.tsv"
nonOTU5_samples <- c("HB0828", "HB0863", "PL0066", "PL0108", "PL0203", "PL0258", "DC3000", "p12.A9", "p12.H7", 
                     "p13.D5", "p13.F1", "p13.F5", "p20.B10", "p23.A5", "p23.B4", "p26.C10", "p26.F6", "p27.C5", 
                     "p3.G11", "p4.A8", "p4.C5", "p4.D11", "p4.H3", "p5.A5", "p5.F2", "p6.B5", "p6.D10", "p6.E9", 
                     "p6.F1", "p6.G2", "p6.G3", "p7.F2", "p8.D11", "p8.G10", "p8.G2", "p9.H10")

# Read merged coverage file
merged_coverage <- read_tsv(merged_coverage_file, col_names = c("Sample", "Position", "Coverage"))

# Read cumulative positions file
cumulative_positions <- read.table(cumulative_positions_file, header = TRUE)

# Adjust positions in merged_coverage by adding gaps between chunks
chunk_separation <- 1000
adjusted_positions <- merged_coverage

# Incrementally add gaps after each segment
for(i in 2:nrow(cumulative_positions)) {
  start <- cumulative_positions$Start[i]+ (i-2)*chunk_separation
  end <-cumulative_positions$End[i]+ (i-2)*chunk_separation
  adjusted_positions$Position[adjusted_positions$Position >= start] <- adjusted_positions$Position[adjusted_positions$Position >= start] + chunk_separation
}

# Calculate coverage proportion per sample
coverage_proportion <- adjusted_positions %>%
  group_by(Sample) %>%
  summarize(CoverageProportion = sum(!is.na(Coverage)) / 18057)

# Filter samples with coverage proportion >= 0.65
filtered65 <- coverage_proportion %>% filter(CoverageProportion >= 0.65)

# Plot coverage proportion
coverage_plot <- ggplot(coverage_proportion, aes(x = reorder(Sample, CoverageProportion), y = CoverageProportion)) +
  geom_col() +
  labs(title = "All 87 Samples Coverage Proportion", x = "", y = "") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  geom_hline(yintercept = 0.65, color = "red", linetype = "dashed")

ggsave('all87allcoveredproportion_markexceptlongest.png', plot = coverage_plot)

# Write filtered coverage proportion to file
write.table(filtered65, 'filter65average_covered_30h_57m_tailocin_markexceptlongest.txt', quote = FALSE, row.names = TRUE, col.names = FALSE)

# Order samples by coverage proportion (descending)
ordered_samples <- coverage_proportion %>%
  arrange((CoverageProportion)) %>%
  pull(Sample)

# Reorder merged_coverage dataframe based on ordered samples
adjusted_positions <- adjusted_positions %>%
  mutate(Sample = factor(Sample, levels = ordered_samples))

# Base plot for alignment coverage with custom colors for nonOTU5 samples
coverage_base_plot <- ggplot(adjusted_positions, aes(x = Position, y = Sample)) +
  geom_point(size = 1, aes(color = case_when(
    Sample %in% nonOTU5_samples ~ 'red',
    grepl('^p', Sample) ~ 'lightgreen',
    TRUE ~ 'lightblue'
  ))) +
  scale_color_identity() +
  labs(title = "Alignments of tailocin (p25.C2), the longest TFA and HFT", x = "Position", y = "Coverage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12), 
        axis.text.y = element_text(size = 12),
        plot.title = element_text(size = 16))

# Add vertical dashed lines and labels for chunks
for(i in seq_along(cumulative_positions$Segment)) {
  start <- cumulative_positions$Start[i] + (i-1) * chunk_separation
  end <- cumulative_positions$End[i] + (i-1) * chunk_separation
  chunk_name <- cumulative_positions$Segment[i]
  
  coverage_base_plot <- coverage_base_plot +
    geom_vline(xintercept = start, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = end, linetype = "dashed", color = "grey") +
    annotate("text", x = (start + end) / 2, y = max(as.numeric(as.factor(adjusted_positions$Sample))) - 25, 
             label = chunk_name, angle = 90, vjust = 1, color = "grey", size = 5)
}

# Add the specific grey lines and labels at 2216, 2728, 2739, and 4541
x_coords <- c(2216, 2728, 2739, 4541)

for (i in 1:length(x_coords)) {
  coverage_base_plot <- coverage_base_plot +
    geom_vline(xintercept = x_coords[i], linetype = "dashed", color = "grey") +
    annotate("text", x = x_coords[i], y = max(as.numeric(as.factor(adjusted_positions$Sample))) - 25, 
             label = labels[i], angle = 90, vjust = 1, color = "grey", size = 5)
}

# Save the final plot with chunks
ggsave("mandh_alignment_coverage_plot_with_chunks_markexceptlongest.png", plot = coverage_base_plot, width = 15, height = 12)