# R script to create a point plot of Length of TM Highly Diverse Region vs Full Length of TM

# Load necessary libraries
library(ggplot2)

# Set working directory
setwd('/Users/cuijiajun/Desktop/others/tmphernan/talia-backman/comparetapemeasurelengthinmydata/')

# Read the data from the file
file_path <- "lengthdiff.txt"
data <- read.table(file_path, header = FALSE, fill = TRUE, col.names = c("Sample", "Length_TM_Diverse", "Sample_Dup", "Full_Length_TM"), na.strings = "")

# Remove duplicate sample column and rows with missing values
data$Sample_Dup <- NULL
data <- na.omit(data)

# Count the number of points for each (x, y) combination
data_count <- as.data.frame(table(data$Length_TM_Diverse, data$Full_Length_TM))
names(data_count) <- c("Length_TM_Diverse", "Full_Length_TM", "Count")
data_count$Length_TM_Diverse <- as.numeric(as.character(data_count$Length_TM_Diverse))
data_count$Full_Length_TM <- as.numeric(as.character(data_count$Full_Length_TM))

# Filter out points with a count of 0
data_count <- data_count[data_count$Count > 0, ]

# Create the point plot with customized labels, font size, and exact count labels
plot <- ggplot(data_count, aes(x = Length_TM_Diverse, y = Full_Length_TM, size = Count)) +
  geom_point(color = "#580", alpha = 0.7) +
  geom_text(aes(label = Count), vjust = -1, hjust = 0.5, size = 4) +  # Add exact count labels
  scale_x_continuous(breaks = c(590, 344, 467)) +
  scale_y_continuous(breaks = c(2526, 2280, 0)) +  # Include y = 0 in the y-axis
  labs(title = "Point Plot of Length of TM Highly Diverse Region vs Full Length of TM",
       x = "Length of TM Highly Diverse Region",
       y = "Full Length of TM") +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

# Print the plot
print(plot)

# Save the plot
ggsave('lengthr.png', plot)
