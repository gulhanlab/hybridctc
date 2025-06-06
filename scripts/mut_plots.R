library(ggplot2)
library(dplyr)
library(patchwork)

# Load and preprocess data
data <- read.csv('/Users/ankit/Desktop/HybridCTC/Ezgi/DNA_Blood_Ctrl_metadata.csv')#, fileEncoding = "UTF-8")
data <- data[data$cluster_final_new != "BadQual", ]
data <- data %>% 
  arrange(old.ident, cluster_final_new, dna_id)  # Arrange by old.ident, then cluster_final_new, then dna_id
data$dna_id <- factor(data$dna_id, levels = unique(data$dna_id))  # Preserve the ordering in dna_id

# Define color palettes
color_palette <- c(
  "CTC" = "#7fc97f", 
  "WBC" = "#beaed4", 
  "Duplet" = "#fdc086", 
  "WBC-like DP" = "#ffff99", 
  "CTC-like DP" = "#386cb0", 
  "Duplet-like DP" = "#f0027f" 
  #"TC_platelet" = "#bf5b17"
)
annotation_palette <- c(
  "Sample_7" = "#7fc97f", 
  "Sample_8" = "#beaed4",
  "Sample_9" = "#fdc086",
  "Sample_10" = "#666666",
  "Ctrl_Pair" = "#ffff99", 
  "Sample_4" = "#386cb0", 
  "Sample_5" = "#f0027f", 
  "Sample_6" = "#bf5b17",
  "Ctrl_Plt_DP" = "#a6cee3",
  "Sample_11" = "#1f78b4"
)

# Create the bar plots
plot1 <- ggplot(data, aes(x = dna_id, y = ctc_fraction_mutation, fill = cluster_final_new)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "CTC Fraction Mutation", x = "DNA ID", y = "Fraction", fill = "Cluster")

plot2 <- ggplot(data, aes(x = dna_id, y = wbc_fraction_mutation, fill = cluster_final_new)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = color_palette) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "WBC Fraction Mutation", x = "DNA ID", y = "Fraction", fill = "Cluster")

# Create the annotation tile
annotation_tile <- ggplot(data, aes(x = dna_id, y = 1, fill = old.ident)) +
  geom_tile() +
  scale_fill_manual(values = annotation_palette) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),  # Remove text from the annotation tile
    axis.ticks.x = element_blank()  # Remove ticks from the annotation tile
  ) +
  labs(fill = "Sample")

# Combine plots with the annotation tile at the bottom
combined_plot <- (plot1 / plot2 / annotation_tile) +  # Stack vertically
  plot_layout(heights = c(4, 4, 0.3))  # Adjust relative heights of plots

# Print the combined plot
print(combined_plot)