# Load required libraries
library(ape)  # For phylogenetic tree manipulation
library(ggtree)  # For tree visualization
library(ggplot2)
library(gridExtra)  # For combining plots
library(grid)
# Set working directory (adjust accordingly)
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/talia-backman-Ps1524_tailocin-520ad97/')

# Read the tree file
tree <- read.tree("/Users/cuijiajun/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/phylogeny_snp/CFML/anew2024/filteredfasta_std_94/filteredfasta_std_94.labelled_tree.newick")

# Read the OTU file
otu_data <- read.table("/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/extract_subgroups/modified_cds_53_lengths_sorted.txt", header = TRUE, sep = "\t")

# Remove the "Seq####_" prefix from OTU sample names
otu_data$Clean_Sample <- sub("_tape", "", otu_data$Sample)

# Replace missing lengths with 0 for unmatched samples
all_samples <- data.frame(Clean_Sample = tree$tip.label)  # Include all samples from the tree
matched_otu_data <- merge(all_samples, otu_data, by = "Clean_Sample", all.x = TRUE)  # Keep all tree labels
matched_otu_data$Length[is.na(matched_otu_data$Length)] <- 2525  # Replace missing lengths with 0
matched_otu_data$Length_Status <- ifelse(matched_otu_data$Length == 2525, "Missing", "Present")

# Adjust the Length column (subtract 2000, ensuring no negative values)
matched_otu_data$Length <- pmax(0, matched_otu_data$Length - 2000)



# Load required libraries
library(ape)  # For phylogenetic tree manipulation
library(ggtree)  # For tree visualization
library(ggplot2)
library(gridExtra)  # For combining plots
library(phangorn)

# Midpoint root the tree
tree <- midpoint(tree)

# Extract tip labels and identify modern (with "^p") and historical samples
# Identify modern (with "^p") and historical samples based on tip labels
modern_samples <- grepl("^p", tree$tip.label)
historical_samples <- !modern_samples

# Create a data frame to associate labels with types
label_data <- data.frame(
  Label = tree$tip.label,
  LabelType = ifelse(modern_samples, "modern", "historical")
)

# Define colors for labels
label_colors <- c("modern" = "green4", "historical" = "orange2")  # Modern as green, historical as orange
# Filter historical labels that appear in otu_data
#historical_labels_in_otu <- label_data$Label[
#  label_data$LabelType == "historical" & label_data$Label %in% otu_data$Clean_Sample
#]
historical_labels_in_otu <- label_data$Label[
 label_data$Label %in% otu_data$Clean_Sample
]
# Create the tree plot with labeled historical samples
tree_plot <- ggtree(tree) +
  theme_tree2() +
  geom_tiplab(size = 0, align = FALSE) +  # Remove textual tip labels
  geom_point2(aes(subset = isTip & label %in% label_data$Label, 
                  color = label_data$LabelType[match(label, label_data$Label)]), 
              size = 1.5) +  # Add colored dots for tips
  geom_text2(aes(subset = isTip & label %in% historical_labels_in_otu, 
                 label = label), 
             hjust = -0.2, size = 2.5) +  # Add text labels for selected historical samples
  scale_color_manual(values = label_colors) +  # Apply custom colors
  theme(
    legend.position = "left",  # Adjust legend position
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 10),  # Adjust legend text size
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    
  )




# Save the updated tree plot with labels as dots
pdf("mydata53_tree_midpoint_root_with_colored_labels.pdf", width = 6, height = 6)
print(tree_plot)
dev.off()
# Extract the order of tree tip labels from 
ggtree
tree_plot_obj <- ggtree(tree)  # Temporary object for label extraction
tree_label_order <- get_taxa_name(tree_plot_obj)  # Extract the actual plotting order

# Reorder matched_otu_data based on the plotting order of the tree
matched_otu_data <- matched_otu_data[match(tree_label_order, matched_otu_data$Clean_Sample), ]
row.names(matched_otu_data) <- seq(1, nrow(matched_otu_data))

# Create the bar plot
bar_plot <- ggplot(matched_otu_data, aes(x = factor(Clean_Sample, levels = rev(matched_otu_data$Clean_Sample)), 
                                         y = Length, fill = Length_Status)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +  # Flip coordinates to match tree orientation
  labs(y = "Adjusted Length", x = "") +
  theme_minimal() +  # Minimal theme for no background
  scale_fill_manual(values = c("Present" = "darkblue", "Missing" = "grey90")) +  # Different colors for presence and absence
  theme(
    axis.text.y = element_blank(),  # Remove y-axis text
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),  # Remove grid lines
    panel.background = element_blank(),  # Remove panel background
    plot.background = element_blank()  # Remove plot background
  )


# Combine the plots side by side
combined_plot <- grid.arrange(tree_plot, bar_plot, ncol = 2, widths = c(3, 1))

# Save the combined figure to a PDF
pdf("mydata_53tree_with_length_barplot_updated.pdf", width = 20, height = 8)
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
pdf("mydata53_tree_with_labels_only.pdf", width = 6, height = 6)  # Adjust size as needed
print(tree_plot_with_labels)
dev.off()
