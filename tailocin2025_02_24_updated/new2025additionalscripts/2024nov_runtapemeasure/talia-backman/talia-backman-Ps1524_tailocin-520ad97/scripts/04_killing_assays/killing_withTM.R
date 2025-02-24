# Load required libraries
library(ggplot2)
library(ape)
library(phytools)
library(ggtree)
library(dplyr)
library(viridis)
library(gridExtra) # For combining plots

# Set working directory
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/talia-backman-Ps1524_tailocin-520ad97/')

# Read in phylogeny tree
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")

# Read in killing assay data
dat2 <- read.csv("./data/killing_assay_data/3replicates_killing_assay_dat.csv")

# Subset tree to only the tester strains used in dat2
uneeded <- which(!(tree$tip.label %in% dat2$Tester_Strain))
tree <- drop.tip(tree, uneeded)

# Clean and process killing assay data
dat2 <- na.omit(dat2)
dat2 <- dat2[, c(1, 2, 7)] # Keep relevant columns
dat2$Final_qualitative[dat2$Final_qualitative == "Resistant"] <- 0
dat2$Final_qualitative[dat2$Final_qualitative == "Sensitive"] <- 1
dat2$Final_qualitative[dat2$Final_qualitative == "Partially Sensitive"] <- 0
dat2$Final_qualitative[dat2$Final_qualitative == "Mostly Sensitive"] <- 1
dat2$Final_qualitative <- as.numeric(dat2$Final_qualitative)

# Transform data to wide format
wide_dat2 <- reshape(dat2, idvar = "Tester_Strain", timevar = "Tailocin_donor_strain", direction = "wide")
rownames(wide_dat2) <- wide_dat2[, 1]
wide_dat2 <- wide_dat2[, -1] # Remove first column
colnames(wide_dat2) <- sub("Final_qualitative.", "", colnames(wide_dat2))

# Read in OTU length data
otu_data <- read.table("OTU52526and2403and2280.txt", header = TRUE, sep = "\t")
otu_data$Clean_Sample <- sub("Seq[0-9]+_", "", otu_data$Sample)

# Match length data to tree labels
all_samples <- data.frame(Clean_Sample = tree$tip.label)
matched_otu_data <- merge(all_samples, otu_data, by = "Clean_Sample", all.x = TRUE)
matched_otu_data$Length[is.na(matched_otu_data$Length)] <- 0 # Replace missing lengths with 0
matched_otu_data$Length <- pmax(0, matched_otu_data$Length - 2000) # Adjust lengths
rownames(matched_otu_data) <- matched_otu_data$Clean_Sample

# Reorder matched_otu_data to match tree tip labels
matched_otu_data <- matched_otu_data[match(tree$tip.label, matched_otu_data$Clean_Sample), ]

# Extract reordered length data
length_data <- matched_otu_data$Length

# Replace NAs with 0 for missing lengths
length_data[is.na(length_data)] <- 0

# Filter the matrix for only the p25.A12 column
dat_matrix2 <- as.matrix(wide_dat2)
dat_matrix2 <- dat_matrix2[, "p25.A12", drop = FALSE]
# Plot the phylogenetic tree
p <- ggtree(tree) + 
  geom_tiplab(size = 2) + 
  theme_tree2() + 
  theme(axis.line.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())


# Create the ggtree object with the heatmap
dat_matrix2 <- as.matrix(wide_dat2)
dat_matrix2 <- dat_matrix2[, "p25.A12", drop = FALSE]
p3 <- gheatmap(p, dat_matrix2, colnames = FALSE, legend_title = "Results", offset = 0.05, color = "black", 
               colnames_position = "top", font.size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  scale_fill_viridis(discrete = FALSE)
# Extract the reordered tip labels from the ggtree object
reordered_tips <- p3$data$label[order(p3$data$y, na.last = NA)]

# Reorder `matched_otu_data` based on the reordered tip labels
matched_otu_data <- matched_otu_data[match(reordered_tips, matched_otu_data$Clean_Sample), ]

# Extract reordered length data
length_data <- matched_otu_data$Length

# Replace NAs with 0 for missing lengths
length_data[is.na(length_data)] <- 0

# Create a data frame for the length bar plot, ensuring the order matches the reordered tree tips
length_df <- data.frame(
  Sample = factor(reordered_tips, levels = reordered_tips), # Match reordered tree tip order
  Length = length_data
)
# Create the bar plot for length distribution with an improved color palette
bar_plot <- ggplot(length_df, aes(x = Sample, y = Length, fill = Length)) +
  geom_bar(stat = "identity", width = 0.8) +
  coord_flip() +
  scale_fill_gradientn(
    colors = c( "#3d304f","#542c91", "#f0f06e"), # Improved gradient
    name = "Length"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),  # Align with tree
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank()
  ) +
  labs(y = "Length (adjusted)", x = "")

# Combine the tree + heatmap with the length bar plot
combined_plot <- grid.arrange(p3, bar_plot, ncol = 2, widths = c(2, 1))

# Save the combined plot as a PDF
pdf("killing_assays_p25A12_with_length_barplot_improved_colors.pdf", width = 12, height = 8)
grid.draw(combined_plot)
dev.off()