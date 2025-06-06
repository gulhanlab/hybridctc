library(ggplot2)
library(dplyr)
library(patchwork)
#library(readxl)
data <- read.csv('/Users/ankit/Desktop/HybridCTC/Ezgi/Sample11/DNA_Sample11.csv')#, fileEncoding = "UTF-8")
data <- data[data$cluster_final != "BadQual", ]
data <- data %>% arrange(cluster_final, Sample)
# Reorder dna_id by cluster_final
data$Sample <- factor(data$Sample, levels = unique(data$Sample))

# Define color palette for categories
color_palette <- c("CTC" = "#7fc97f", "WBC" = "#beaed4", "Duplet" = "#fdc086", "WBC-like DP" = "#f0027f", "CTC-like DP" = "#386cb0", "BadQual" = "#ffff99", "TC_platelet" = "#bf5b17")

# Create the bar plots
plot1 <- ggplot(data, aes(x = Sample, y = ctc_fraction_mutation, fill = cluster_final)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "CTC Fraction Mutation", x = "DNA ID", y = "Fraction", fill = "Cluster")

plot2 <- ggplot(data, aes(x = Sample, y = wbc_fraction_mutation, fill = cluster_final)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "WBC Fraction Mutation", x = "DNA ID", y = "Fraction", fill = "Cluster")

# Arrange the two plots in a single row (2x1 grid)
combined_plot <- plot1 + plot2 + plot_layout(ncol = 1)
print(combined_plot)
