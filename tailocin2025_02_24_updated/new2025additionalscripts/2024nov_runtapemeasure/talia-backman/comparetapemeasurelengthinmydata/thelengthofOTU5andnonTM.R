# Load required libraries
library(ggplot2)

# Create the data frame with the provided data
data <- data.frame(
  Group = c(rep("OTU5", 8), rep("nonOTU5", 5)),
  Isolates = c(1, 445, 2, 898, 2, 1, 1, 1, 19, 9, 1, 1, 1),
  TM_Length = c(1797, 2280, 2403, 2526, 624, 765, 864, 315, 2142, 2157, 315, 2526, 315)
)

# Order the data by TM_Length
data <- data[order(data$TM_Length), ]

# Create the bar plot with bold and larger text
p <- ggplot(data, aes(x = factor(TM_Length, levels = unique(TM_Length)), y = Isolates, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("OTU5" = "#599", "nonOTU5" = "#958")) +
  labs(x = "Length of TM", y = "Number of Isolates", title = "Bar Plot of Number of Isolates by TM Length") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Bold and larger title
    axis.title.x = element_text(size = 16, face = "bold"),  # Bold and larger x-axis label
    axis.title.y = element_text(size = 16, face = "bold"),  # Bold and larger y-axis label
    axis.text = element_text(size = 14),  # Larger axis tick labels
    legend.title = element_text(size = 14, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 12)  # Larger legend text
  )

# Save the plot
ggsave('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/comparetapemeasurelengthinmydata/OTU5lengthTM.png', p, width = 10, height = 6)




#
# Install and load required libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  install.packages("scales")
}

library(dplyr)
library(ggplot2)
library(scales)

# Create the data frame with the provided data
data <- data.frame(
  Group = c(rep("OTU5", 8), rep("nonOTU5", 5)),
  Isolates = c(1, 445, 2, 898, 2, 1, 1, 1, 19, 9, 1, 1, 1),
  TM_Length = c(1797, 2280, 2403, 2526, 624, 765, 864, 315, 2142, 2157, 315, 2526, 315)
)

# Filter only the OTU5 group
otu5_data <- subset(data, Group == "OTU5")

# Mark specific lengths as "rare" if not 2526, 2280, or 2403
otu5_data$TM_Length <- ifelse(
  otu5_data$TM_Length %in% c(2526, 2280, 2403), 
  as.character(otu5_data$TM_Length), 
  "Rare"
)

# Order the data to maintain the custom levels
otu5_data$TM_Length <- factor(otu5_data$TM_Length, levels = unique(otu5_data$TM_Length))

# Calculate proportions for the y-axis
otu5_data <- otu5_data %>%
  group_by(TM_Length) %>%
  summarise(Total_Isolates = sum(Isolates)) %>%
  mutate(Proportion = Total_Isolates / sum(Total_Isolates))

# Create the bar plot with bold and larger text
p <- ggplot(otu5_data, aes(x = TM_Length, y = Proportion, fill = TM_Length)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("2526" = "darkblue", "2280" = "darkblue", "2403" = "darkblue", "Rare" = "#bbb")) +
  scale_y_continuous(labels = percent, limits = c(0, 1)) +  # Set y-axis as percentage with 100% limit
  labs(x = "Length of TM", y = "Proportion of Isolates", title = "Proportion of OTU5 Isolates by TM Length") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),  # Bold and larger title
    axis.title.x = element_text(size = 16, face = "bold"),  # Bold and larger x-axis label
    axis.title.y = element_text(size = 16, face = "bold"),  # Bold and larger y-axis label
    axis.text = element_text(size = 24),  # Larger axis tick labels
    legend.title = element_text(size = 14, face = "bold"),  # Bold legend title
    legend.text = element_text(size = 22)  # Larger legend text
  )

# Save the plot
ggsave('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/comparetapemeasurelengthinmydata/OTU5lengthTM_proportion.png', 
       p, width = 10, height = 6)