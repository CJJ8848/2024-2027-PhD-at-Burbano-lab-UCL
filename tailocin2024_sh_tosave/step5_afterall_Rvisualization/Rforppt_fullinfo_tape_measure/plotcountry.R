# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

setwd('/Users/cuijiajun/Desktop/2023-2024\ PhD\ ucl/2024_aMeta/wholepipeAt_Ps/2024_233_analysis/tailocin46/tree/fullinfo/fullinfo/')

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the country only data
country_data <- read.table("country_only.txt", header = FALSE, col.names = "country")
# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the country only data
# Create a column for the x-axis values
country_data$x <- 1


# Create a sequential index for each row to ensure all 68 entries are shown
country_data$index <- 1:nrow(country_data)
# Plot the data
ggplot(country_data, aes(x = x, y = factor(index, levels = rev(index)), color = country)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Country Distribution", x = "", y = "Country") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
ggsave('country.pdf')