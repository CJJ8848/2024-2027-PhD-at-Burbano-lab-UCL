# Load required libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ape)

# Set file paths
tree_file_path <- "/Users/cuijiajun/Desktop/2023-2024 PhD ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/tailocin46/2024_oct_summaryfortalia/tree/coloredmuscle_aligned_33_7refsHTF.fasta.treefile"
depth_matrix_path <- "method2_afterassembly_gene_coverage_matrix.tsv"

# === Step 1: Load Tree and Modify Sample Names ===
tree <- ape::read.nexus(tree_file_path)

# Check if tree labels exist
if (is.null(tree$tip.label) || all(is.na(tree$tip.label))) {
  stop("Error: No valid tip labels found in the NEXUS tree. Check tree format.")
}

# Remove extra single quotes in tip labels
tree$tip.label <- gsub("'", "", tree$tip.label)

# **Convert `_` to `.` only for labels that start with "p"**
tree$tip.label[grepl("^p", tree$tip.label)] <- gsub("_", ".", tree$tip.label[grepl("^p", tree$tip.label)])

# Extract modified sample names
tree_sample_names <- tree$tip.label
print("Extracted sample names from tree:")
print(tree_sample_names)

# === Step 2: Load Depth Matrix and Subset for Tree Samples ===
depth_df <- read.table(depth_matrix_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Check if sample names exist in depth matrix
valid_samples <- tree_sample_names[tree_sample_names %in% colnames(depth_df)]

if (length(valid_samples) == 0) {
  stop("Error: No valid samples from the tree found in the depth matrix.")
}

# Subset depth matrix
filtered_depth_df <- depth_df[, valid_samples, drop=FALSE]

# Remove tagG_2 and tagH_2 if they exist
filtered_depth_df <- filtered_depth_df[!rownames(filtered_depth_df) %in% c("tagG_2", "tagH_2"), ]

# === Step 3: Order Samples Exactly as Requested ===
ordered_samples <- c(
  "p24.H2", "p8.H7", "p13.F3", "p23.A3", "p7.G11", "HTF_p7.G11",
  "p4.D2", "p12.G7", "p3.G9", "PL0080", "p13.D10", "HTF_p25.A12",
  "p8.E4", "PL0102", "PL0042", "p27.F2", "HTF_p5.D5", "PL0131",
  "109.NOR_1990", "p8.B9", "p6.B9", "p26.E7", "HTF_p26.D6", "HTF_p23.B8",
  "p3.F8", "p12.E2", "p26.B7", "p6.A10", "HTF_p25.C2", "HB0766",
  "p27.D6", "76.LTU_2009_S19", "p25.C11", "p22.D1", "p20.G9",
  "p13.C1", "HTF_p21.F9", "PL0240", "PL0059", "p8.C7"
)

# Ensure the order matches existing valid samples
ordered_samples <- ordered_samples[ordered_samples %in% colnames(filtered_depth_df)]
filtered_depth_df <- filtered_depth_df[, ordered_samples]

# === Step 4: Assign HTF Colors and Sample Annotations ===
htf_samples <- list(
  "HTF_p23.B8" = "1803", "HTF_p26.D6" = "1803", "HTF_p21.F9" = "1245",
  "HTF_p25.C2" = "1803", "HTF_p5.D5" = "1803", "HTF_p25.A12" = "1383", "HTF_p7.G11" = "1830"
)

htf_colors <- c(
  "1245" = "#FFD700",  # Yellow
  "1383" = "#1E90FF",  # Blue
  "1803" = "#FF0000",  # Red
  "1830" = "#FF69B4"   # Pink
)

# Identify all sample names in the dataset
sample_names <- colnames(filtered_depth_df)

# Create a dataframe for annotation
annotation_col <- data.frame(Group_Label = rep("Other", length(sample_names)), row.names = sample_names)

# Assign HTF length categories to matching samples
for (sample in names(htf_samples)) {
  if (sample %in% sample_names) {
    annotation_col[sample, "Group_Label"] <- htf_samples[[sample]]  # Assign HTF length
  }
}

# Identify "modern" (p-starting) and "historical" samples for non-HTF cases
annotation_col$Group_Label <- ifelse(
  annotation_col$Group_Label == "Other",  # If not assigned an HTF length
  ifelse(grepl("^p", rownames(annotation_col)), "modern", "historical"),  # Assign modern or historical
  annotation_col$Group_Label  # Otherwise, keep the HTF length
)

# Define a **unified** color scheme for both HTF length and historical/modern labels
group_colors <- c(
  "1245" = "#FFD700",  # Yellow
  "1383" = "#1E90FF",  # Blue
  "1803" = "#FF0000",  # Red
  "1830" = "#FF69B4",  # Pink
  "modern" = "#2C2B100C",  # Grey for modern samples
  "historical" = "#BEBEBE"  # Light grey for historical samples
)

# Convert to factor to ensure consistent ordering in pheatmap
valid_levels <- names(group_colors)  # Ensure levels match the color mapping
annotation_col$Group_Label <- factor(annotation_col$Group_Label, levels = valid_levels)

# Ensure annotation_colors matches the factor levels
annotation_colors <- list(Group_Label = group_colors)

# === Step 5: Adjust Margins and Generate Heatmap ===
if (ncol(filtered_depth_df) == 0) {
  stop("Error: No valid samples remain after filtering. Cannot generate heatmap.")
}

# Increase figure width dynamically and add **extra left margin**
fig_width <-  ncol(filtered_depth_df) * 0.5+ 2  # Increased extra margin

# ✅ Manually set the gene order
desired_gene_order <- c("WfgD", "rmlC_1", "tagG_1", "tagH_1", "epsE_2","spsA")

# ✅ Reorder filtered_depth_df to match the desired gene order
filtered_depth_df <- filtered_depth_df[desired_gene_order, , drop=FALSE]  # Drop=FALSE to keep matrix structure

# Generate heatmap with HTF and sample annotations
pdf("./plot/method2/m2Ordered_Heatmap_With_Colors.pdf", width=fig_width, height=8)

pheatmap(filtered_depth_df, 
         cluster_rows=FALSE,  # **No tree structure**
         cluster_cols=FALSE,  # **No clustering to preserve order**
         color = colorRampPalette(c("#FFFFFF", "#BEBEBE", "#2C2C2C"))(100),
         display_numbers = FALSE,
         main = "Heatmap of Coverage Proportion for Assembled Genes",
         fontsize_col = 16,  # Larger font for better visibility
         fontsize_row = 14,  
         angle_col = 90,
         border_color = "grey80",
         labels_col = ordered_samples,
         annotation_col = annotation_col["Group_Label"],  
         annotation_colors = annotation_colors,
         annotation_names_col = FALSE,  
         annotation_legend = TRUE)
dev.off()

print("✅ Heatmap successfully saved as Ordered_Heatmap_With_Colors.pdf!")