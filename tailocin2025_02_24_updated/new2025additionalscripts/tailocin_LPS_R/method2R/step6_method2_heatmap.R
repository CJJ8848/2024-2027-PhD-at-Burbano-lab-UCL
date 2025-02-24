# Load required libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

# Set working directory
setwd("/Users/cuijiajun/Desktop/others/tmphernan/2025/tailocin_LPS")

# Set file path
file_path <- "method2_afterassembly_gene_coverage_matrix.tsv"

# Load data
df <- read.table(file_path, header=TRUE, sep="\t", row.names=1, check.names=FALSE)

# Convert to matrix
data_matrix <- as.matrix(df)

# Replace NAs with 0 to avoid missing values in the heatmap
data_matrix[is.na(data_matrix)] <- 0

# Define HTF samples and their lengths
htf_samples <- list(
  "HTF_p23.B8" = "1803", "HTF_p26.D6" = "1803", "HTF_p21.F9" = "1245",
  "HTF_p25.C2" = "1803", "HTF_p5.D5" = "1803", "HTF_p25.A12" = "1383", "HTF_p7.G11" = "1830"
)

# Define colors for HTF lengths
htf_colors <- c(
  "1245" = "#FFD700",  # Yellow
  "1383" = "#1E90FF",  # Blue
  "1803" = "#FF0000",  # Red
  "1830" = "#FF69B4"   # Pink
)

# Identify all sample names in the dataset
sample_names <- colnames(data_matrix)

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
  "historical" = "#BEBEBE"  # Dark red for historical samples
)

# Convert to factor to ensure consistent ordering in pheatmap
valid_levels <- names(group_colors)  # Ensure levels match the color mapping
annotation_col$Group_Label <- factor(annotation_col$Group_Label, levels = valid_levels)

# Ensure annotation_colors matches the factor levels
annotation_colors <- list(Group_Label = group_colors)

# Set figure width dynamically
fig_width <- max(14, ncol(data_matrix) * 0.4)

# Define grayscale heatmap colors
heatmap_colors <- colorRampPalette(c("#FFFFFF", "#BEBEBE", "#2C2C2C"))(100)  # White → Light Grey → Dark Grey


# ✅ Manually set the gene order
desired_gene_order <- c("WfgD", "rmlC_1", "tagG_1", "tagH_1", "epsE_2", "tagG_2", "tagH_2", "spsA")

# ✅ Reorder data_matrix to match the desired gene order
data_matrix <- data_matrix[desired_gene_order, , drop=FALSE]  # Drop=FALSE to keep matrix structure

# === Save FULL heatmap ===
pdf("./plot/method2/m2Relative_Read_Depth_HTF_Fixed_Custom_Order.pdf", width=fig_width, height=8)
pheatmap(data_matrix, cluster_rows=FALSE, cluster_cols=TRUE,  
         color = heatmap_colors,  
         display_numbers = FALSE,  
         main = "Heatmap of Coverage Proportion for Assembled Genes",
         fontsize_col = 14,  
         fontsize_row = 14,  
         angle_col = 90,  
         border_color = "grey80",
         labels_col = sample_names,  
         annotation_col = annotation_col["Group_Label"],  # Use only 1 annotation row
         annotation_colors = annotation_colors,  # Use unified color scheme
         fontsize = 16,
         annotation_names_col = FALSE,  # Remove default column annotation names
         annotation_legend = TRUE)  # Show annotation legend
dev.off()

print("✅ Heatmap with **custom gene order** saved as Relative_Read_Depth_HTF_Fixed_Custom_Order.pdf!")

# === Remove tagG_2 and tagH_2 and Save Filtered Heatmap ===
filtered_matrix <- data_matrix[!rownames(data_matrix) %in% c("tagG_2", "tagH_2"), , drop=FALSE]

pdf("./plot/method2/m2Relative_Read_Depth_HTF_Fixed_No_tagG2_tagH2_Custom_Order.pdf", width=fig_width, height=8)
pheatmap(filtered_matrix, cluster_rows=FALSE, cluster_cols=TRUE,  
         color = heatmap_colors,  
         display_numbers = FALSE,  
         main = "Heatmap of Coverage Proportion for Assembled Genes (No tagG_2 & tagH_2)",  
         fontsize_col = 14,  
         fontsize_row = 14,  
         angle_col = 90,  
         border_color = "grey80",
         labels_col = sample_names,        
         annotation_col = annotation_col["Group_Label"],  # Use only 1 annotation row
         annotation_colors = annotation_colors,  # Use unified color scheme
         fontsize = 16,
         annotation_names_col = FALSE,  # Remove default column annotation names
         annotation_legend = TRUE)  # Show annotation legend
dev.off()

print("✅ Filtered heatmap with **custom gene order** saved as Relative_Read_Depth_HTF_Fixed_No_tagG2_tagH2_Custom_Order.pdf!")