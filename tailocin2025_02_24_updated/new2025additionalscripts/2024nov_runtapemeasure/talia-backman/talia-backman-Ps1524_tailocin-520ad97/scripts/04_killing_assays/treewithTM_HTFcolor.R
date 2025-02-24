# Load required libraries
library(ape)  # For phylogenetic tree manipulation
library(ggtree)  # For tree visualization
library(ggplot2)
library(gridExtra)  # For combining plots
library(grid)
# Set working directory (adjust accordingly)
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/talia-backman-Ps1524_tailocin-520ad97/')

# Read the tree file
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")

# Read the OTU file
otu_data <- read.table("OTU52526and2403and2280.txt", header = TRUE, sep = "\t")

# Remove the "Seq####_" prefix from OTU sample names
otu_data$Clean_Sample <- sub("Seq[0-9]+_", "", otu_data$Sample)

# Replace missing lengths with 0 for unmatched samples
all_samples <- data.frame(Clean_Sample = tree$tip.label)  # Include all samples from the tree
matched_otu_data <- merge(all_samples, otu_data, by = "Clean_Sample", all.x = TRUE)  # Keep all tree labels
matched_otu_data$Length[is.na(matched_otu_data$Length)] <- 0  # Replace missing lengths with 0

# Adjust the Length column (subtract 2000, ensuring no negative values)
matched_otu_data$Length <- pmax(0, matched_otu_data$Length - 2000)

# Extract the order of tree tip labels from ggtree
tree_plot_obj <- ggtree(tree)  # Temporary object for label extraction
tree_label_order <- get_taxa_name(tree_plot_obj)  # Extract the actual plotting order

# Reorder matched_otu_data based on the plotting order of the tree
matched_otu_data <- matched_otu_data[match(tree_label_order, matched_otu_data$Clean_Sample), ]
row.names(matched_otu_data) <- seq(1, nrow(matched_otu_data))
# Create the tree plot without vertical lines or tip labels
tree_plot <- ggtree(tree) +
  theme_tree2() +
  theme(axis.line.x = element_blank(),  # Remove vertical axis lines
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  geom_tiplab(size = 0, align = FALSE)  # Hide tip labels

# Create the bar plot
# Create the bar plot
# Create the bar plot with reversed order
# Create the bar plot with no background and a color gradient
bar_plot <- ggplot(matched_otu_data, aes(x = factor(Clean_Sample, levels = rev(matched_otu_data$Clean_Sample)), y = Length, fill = Length)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +  # Flip coordinates to match tree orientation
  labs(y = "Adjusted Length", x = "") +
  theme_minimal() +  # Minimal theme for no background
  scale_fill_gradient(low = "lightblue", high = "darkblue") +  # Better color gradient
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank()  # Remove plot background
  )

# Combine the plots side by side
combined_plot <- grid.arrange(tree_plot, bar_plot, ncol = 2, widths = c(1.5, 1))

# Save the combined figure to a PDF
pdf("tree_with_length_barplot_updated.pdf", width = 12, height = 8)
grid.draw(combined_plot)
dev.off()
library(ggtree)

# Create the tree plot with labels
tree_plot_with_labels <- ggtree(tree) +
  theme_tree2() +
  geom_tiplab(size = 1, align = F, offset = 0.5) +  # Align labels and add offset to prevent overlap
  theme(axis.line.x = element_blank(),  # Remove vertical axis lines
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# Save the tree plot with labels to a PDF
pdf("tree_with_labels_only.pdf", width = 60, height = 30)  # Adjust size as needed
print(tree_plot_with_labels)
dev.off()

#with pega hyplotypers:

# Load required libraries
library(ape)       # For phylogenetic tree manipulation
library(ggtree)    # For tree visualization
library(ggplot2)
library(gridExtra) # For combining plots




#HTF

# Load the HTF data
HTF <- read.csv('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/talia-backman-Ps1524_tailocin-520ad97/data/supplemental_data/HTF_aa_haplotypes.csv')

# Trim the names in the 'Strain' column to match the desired format
HTF$Clean_Strain <- sub(".*_(p[0-9]+\\.[A-Z0-9]+)", "\\1", HTF$Strain)

# Merge the HTF data with the matched OTU data
matched_otu_data <- merge(matched_otu_data, HTF, by.x = "Clean_Sample", by.y = "Clean_Strain", all.x = TRUE)

# Assign colors for each unique haplotype in the HTF table
htf_haplotypes <- unique(HTF$HP12_haplotype)
htf_colors <- setNames(colorRampPalette(c("blue", "green", "red", "orange", "purple", "pink"))(length(htf_haplotypes)), htf_haplotypes)

# Assign color for each bar based on the HP12_haplotype column
matched_otu_data$Color <- ifelse(is.na(matched_otu_data$HP12_haplotype), "grey", matched_otu_data$HP12_haplotype)

# Reorder matched_otu_data based on the plotting order of the tree
matched_otu_data <- matched_otu_data[match(tree_label_order, matched_otu_data$Clean_Sample), ]
row.names(matched_otu_data) <- seq(1, nrow(matched_otu_data))
# Update the bar plot with HTF haplotype-based coloring
bar_plot_htf <- ggplot(matched_otu_data, aes(
  x = factor(Clean_Sample, levels = rev(matched_otu_data$Clean_Sample)),
  y = Length, fill = Color)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  labs(y = "Adjusted Length", x = "", fill = "HTF Haplotype") +
  theme_minimal() +
  scale_fill_manual(values = c(htf_colors, "grey" = "grey")) +  # Include grey for unmatched
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Combine the tree plot and HTF haplotype-colored bar plot side by side
combined_plot_htf <- grid.arrange(tree_plot, bar_plot_htf, ncol = 2, widths = c(1.5, 1))

# Save the updated combined figure to a PDF
pdf("tree_with_htf_haplotype_barplot.pdf", width = 12, height = 8)
grid.draw(combined_plot_htf)
dev.off()




#length

# Load necessary libraries
library(RColorBrewer)

# Define the length data with HTFlength column
length_data <- data.frame(
  HP12_haplotype = c("p12.E4", "p20.D8", "p21.F9", "p23.B3", "p25.A12", 
                     "p25.C2", "p25.G1", "p26.D6", "p5.D5", "p7.B1"),
  HTFlength = c(1383, 1383, 1245, 1383, 1383, 1803, 1803, 1803, 1803, 1830)
)

# Clean the `HP12_haplotype` column in matched_otu_data to extract haplotype names
matched_otu_data$Cleaned_Haplotype <- sub(",.*", "", matched_otu_data$HP12_haplotype)
matched_otu_data$Cleaned_Haplotype <- sub(" .*", "", matched_otu_data$Cleaned_Haplotype)

# Merge the length data with matched_otu_data using the cleaned haplotype column
matched_otu_data <- merge(matched_otu_data, length_data, 
                          by.x = "Cleaned_Haplotype", by.y = "HP12_haplotype", 
                          all.x = TRUE)

# Assign a color palette for each unique HTFlength using RColorBrewer
unique_lengths <- unique(length_data$HTFlength)
htf_colors <- setNames(brewer.pal(n = min(length(unique_lengths), 9), name = "Set1"), unique_lengths)

# Add a new column for bar colors based on unique HTFlength
matched_otu_data$Length_Color <- factor(matched_otu_data$HTFlength, levels = unique_lengths)
# Reorder matched_otu_data based on the plotting order of the tree
matched_otu_data <- matched_otu_data[match(tree_label_order, matched_otu_data$Clean_Sample), ]
row.names(matched_otu_data) <- seq(1, nrow(matched_otu_data))
# Create the bar plot with colors based on unique HTFlength
bar_plot_htf <- ggplot(matched_otu_data, aes(
  x = factor(Clean_Sample, levels = rev(matched_otu_data$Clean_Sample)),
  y = Length, fill = Length_Color)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  labs(y = "HTFlength", x = "", fill = "Unique HTFlength") +
  theme_minimal() +
  scale_fill_manual(values = htf_colors) +  # Map unique lengths to colors
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Combine the tree plot and HTF-length colored bar plot side by side
combined_plot_htf <- grid.arrange(tree_plot, bar_plot_htf, ncol = 2, widths = c(1.5, 1))

# Save the updated combined figure to a PDF
pdf("tree_with_htf_haplotype_barplot_colorwithlength.pdf", width = 12, height = 8)
grid.draw(combined_plot_htf)
dev.off()






#lengthbar
# Ensure HTFlength and Length are numeric for proper sorting
matched_otu_data$HTFlength <- as.numeric(as.character(matched_otu_data$HTFlength))
matched_otu_data$Length <- as.numeric(as.character(matched_otu_data$Length))

# Reorder data: First by HTFlength (descending), then by Length (descending)
matched_otu_data <- matched_otu_data[order(-matched_otu_data$HTFlength, -matched_otu_data$Length), ]

# Update the Clean_Sample factor levels based on the new order
matched_otu_data$Clean_Sample <- factor(matched_otu_data$Clean_Sample, levels = matched_otu_data$Clean_Sample)

# Create the bar plot with colors based on unique HTFlength
bar_plot_htf <- ggplot(matched_otu_data, aes(
  x = Clean_Sample,
  y = Length,
  fill = Length_Color
)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  labs(y = "Length", x = "", fill = "Unique HTFlength") +
  theme_minimal() +
  scale_fill_manual(values = htf_colors) +  # Map unique lengths to colors
  theme(
    axis.text.y = element_text(size = 10),  # Display y-axis labels for better readability
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Save the plot as a longer plot
pdf("bar_plot_htf_sorted_by_HTFlength_Length.pdf", width = 8, height = 12)  # Longer plot
print(bar_plot_htf)
dev.off()







#box plot
# Check if Length and HTFlength columns exist
if ("Length" %in% colnames(matched_otu_data) && "HTFlength" %in% colnames(matched_otu_data)) {
  # Filter the data for Length 280 and 526
  filtered_data <- subset(matched_otu_data, Length %in% c(280, 526))
  
  # Ensure Length is treated as a factor
  filtered_data$Length <- as.factor(filtered_data$Length)
  
  # Perform statistical test: Wilcoxon rank-sum test
  test_result <- wilcox.test(HTFlength ~ Length, data = filtered_data)  # Non-parametric test for two groups
  p_value <- test_result$p.value  # Extract the p-value
  
  # Create the box plot with points and p-value annotation
  box_plot <- ggplot(filtered_data, aes(x = Length, y = HTFlength)) +
    geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +  # Box plot
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "darkblue") +   # Overlay points
    labs(x = "Length", y = "HTFlength", title = "Box Plot of HTFlength by Length (280 vs. 526)") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    ) +
    annotate("text", 
             x = 1.5, y = max(filtered_data$HTFlength, na.rm = TRUE) * 1.05, 
             label = paste("p =", signif(p_value, 3)), 
             size = 5, color = "red", fontface = "italic")
  
  # Save the plot to a PDF
  pdf("box_plot_HTFlength_vs_Length_280_526.pdf", width = 8, height = 6)
  print(box_plot)
  dev.off()
  
  # Display the p-value in the console
  cat("P-value from the Wilcoxon rank-sum test: ", signif(p_value, 3), "\n")
} else {
  cat("Required columns ('Length' and 'HTFlength') are missing in matched_otu_data.\n")


}

# Check if Length and HTFlength columns exist
if ("Length" %in% colnames(matched_otu_data) && "HTFlength" %in% colnames(matched_otu_data)) {
  # Ensure HTFlength is treated as a factor
  matched_otu_data$HTFlength <- as.factor(matched_otu_data$HTFlength)
  matched_otu_data <- subset(matched_otu_data, Length %in% c(280, 526))
  
  # Perform statistical test: Kruskal-Wallis test
  test_result <- kruskal.test(Length ~ HTFlength, data = matched_otu_data)  # Non-parametric test for multiple groups
  p_value <- test_result$p.value  # Extract the p-value
  
  # Create the box plot with points and p-value annotation
  box_plot <- ggplot(matched_otu_data, aes(x = HTFlength, y = Length)) +
    geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +  # Box plot
    geom_jitter(width = 0.2, size = 2, alpha = 0.7, color = "darkblue") +  # Overlay points
    labs(x = "HTFlength Groups", y = "Length", title = "Box Plot of Length by HTFlength Groups") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 14, face = "bold"),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 16, face = "bold")
    ) 
  
  # Save the plot to a PDF
  pdf("box_plot_Length_vs_HTFlength_groups.pdf", width = 8, height = 6)
  print(box_plot)
  dev.off()
  
  # Display the p-value in the console
  cat("P-value from the Kruskal-Wallis test: ", signif(p_value, 3), "\n")
} else {
  cat("Required columns ('Length' and 'HTFlength') are missing in matched_otu_data.\n")
}



