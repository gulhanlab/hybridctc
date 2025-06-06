library(ggplot2)
library(dplyr)
library(readxl)

# Load data
#meta <- read.delim('/Users/ankit/Desktop/HybridCTC/Ezgi/Sample10/DNA_Sample10.csv')
ave <- read.delim('~/Partners HealthCare Dropbox/Ankit Singh/GulhanLab/Data/hybrid_ctc_2025/dna/ref_files/average_ctc_cdp_cn.txt')
seg <- read.delim('/Users/ankit/Desktop/HybridCTC/Ezgi/Sample11/SegNorm.txt')

# Filter samples
# samples <- meta$Sample[which(meta$ginkgo_annot != "Bad Quality")]

# Choose a single sample for processing
sample_to_process <- 'EZ2073_gDNA_P1_Blood_DP'  # Example sample

# Define which scaled column to use (manually set to one of: "scaled_wbc", "scaled_ctc", "scaled_cdp")
chosen_cell_type <- "scaled_wbc"  # Change this to "scaled_ctc" or "scaled_cdp" as needed

vec_val <- numeric()
vec_id <- character()
vec_ctc <- numeric()
vec_cdp <- numeric()
vec_wbc <- numeric()

cn_ctc <- ave$cn_ctc
cn_cdp <- ave$cn_cdp
cn_wbc <- ave$cn_wbc

# Generate data for scaling calculation
for (val in seq(0, 20, by = 0.1)) {
  cn <- seg[, sample_to_process]
  med <- median(cn)
  cn[cn > med * 10] <- med * 10
  diff_ctc <- cn_ctc - cn * val
  diff_cdp <- cn_cdp - cn * val
  diff_wbc <- cn_wbc - cn * val
  
  vec_ctc <- c(vec_ctc, sum(abs(diff_ctc)))
  vec_cdp <- c(vec_cdp, sum(abs(diff_cdp)))
  vec_wbc <- c(vec_wbc, sum(abs(diff_wbc)))
  
  vec_val <- c(vec_val, val)
  vec_id <- c(vec_id, sample_to_process)
}

df_diff <- data.frame(
  id = vec_id, val = vec_val,
  diff_ctc = vec_ctc,
  diff_wbc = vec_wbc,
  diff_cdp = vec_cdp
)

# Compute the best scaling factors
scaling_ctc <- df_diff %>% filter(diff_ctc == min(diff_ctc)) %>% pull(val)
scaling_cdp <- df_diff %>% filter(diff_cdp == min(diff_cdp)) %>% pull(val)
scaling_wbc <- df_diff %>% filter(diff_wbc == min(diff_wbc)) %>% pull(val)

# Extract relevant segments
seg_filtered <- seg %>%
  select(CHR, START, END, all_of(sample_to_process))

# Compute scaled CN values
seg_filtered$scaled_ctc <- seg_filtered[[sample_to_process]] * scaling_ctc
seg_filtered$scaled_cdp <- seg_filtered[[sample_to_process]] * scaling_cdp
seg_filtered$scaled_wbc <- seg_filtered[[sample_to_process]] * scaling_wbc

# Create output dataframe with selected cell type
output_df <- seg_filtered %>%
  select(CHR, START, END, all_of(sample_to_process), all_of(chosen_cell_type))

# Cap the values in chosen_cell_type column to 10
output_df[[chosen_cell_type]] <- pmin(output_df[[chosen_cell_type]], 10)

# Rename the last column to match the sample name
colnames(output_df)[4] <- sample_to_process
colnames(output_df)[5] <- chosen_cell_type

# Save results as CSV file
output_file <- paste0("/Users/ankit/Desktop/HybridCTC/Ezgi/Sample11/scaled_cnv.csv")
write.csv(output_df, file = output_file, row.names = FALSE, quote = FALSE)