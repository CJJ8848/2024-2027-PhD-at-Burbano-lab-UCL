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

# Load haplotype data
haplo_2280 <- read.table("2280sample_haplotype_mapping.txt", header = TRUE, sep = "\t")
haplo_2526 <- read.table("2526sample_haplotype_mapping.txt", header = TRUE, sep = "\t")

# Combine haplotype data
haplo_combined <- rbind(haplo_2280, haplo_2526)
haplo_combined$Clean_Sample <- haplo_combined$Sample  # Ensure consistent column naming

# Merge haplotype info with matched OTU data
matched_otu_data <- merge(matched_otu_data, haplo_combined, by = "Clean_Sample", all.x = TRUE)
matched_otu_data <- matched_otu_data[match(tree_label_order, matched_otu_data$Clean_Sample), ]
row.names(matched_otu_data) <- seq(1, nrow(matched_otu_data))

# Assign colors: Use distinct colors for haplotypes and grey for others
haplotype_colors <- c(
  "2280_IV" = "#1F78B4",  # Deep Blue
  "2280_II" = "#33A02C",  # Vibrant Green
  "2280_III" = "#E31A1C", # Rich Red
  "2526_III" = "#FF7F00", # Bright Orange
  "2526_IV" = "#B15928",  # Royal Purple
  "2526_I" = "#6A3D9A",   # Earthy Brown
  "grey" = "grey"         # Grey for unmatched
)


# Assign color for each bar
matched_otu_data$Color <- ifelse(is.na(matched_otu_data$Haplotype), "grey",
                                 matched_otu_data$Haplotype)

# Create the bar plot with haplotype-based coloring and labeled legend
bar_plot_haplo <- ggplot(matched_otu_data, aes(
  x = factor(Clean_Sample, levels = rev(matched_otu_data$Clean_Sample)), 
  y = Length, fill = Color)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  labs(y = "Adjusted Length", x = "", fill = "Haplotype") +
  theme_minimal() +
  scale_fill_manual(values = haplotype_colors) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  )

# Combine the tree plot and haplotype-colored bar plot side by side
combined_plot_haplo <- grid.arrange(tree_plot, bar_plot_haplo, ncol = 2, widths = c(1.5, 1))

# Save the updated combined figure to a PDF
pdf("tree_with_haplotype_barplot_with_labels.pdf", width = 12, height = 8)
grid.draw(combined_plot_haplo)
dev.off()