library(ggplot2)
library(dplyr)
library(readxl)

# Directory to save the plots
output_dir <- "/Users/ankit/Desktop/HybridCTC/Ezgi/Sample12/scaled_plots"

meta <- read.csv('/Users/ankit/Desktop/HybridCTC/Ezgi/Sample12/DNA_Sample12.csv')
ave <- read.delim('~/Partners HealthCare Dropbox/Ankit Singh/GulhanLab/Data/hybrid_ctc_2025/dna/ref_files/average_ctc_cdp_cn.txt')
seg <- read.delim('/Users/ankit/Desktop/HybridCTC/Ezgi/Sample12/SegNorm.txt')

#samples <- meta$Sample[which(meta$ginkgo_annot != "Bad Quality")]
samples <- meta$Sample

print(samples)
print(colnames(seg))

vec_val <- numeric()
vec_id <- character()
vec_ctc <- numeric()
vec_cdp <- numeric()
vec_wbc <- numeric()

cn_ctc <- ave$cn_ctc
cn_cdp <- ave$cn_cdp
cn_wbc <- ave$cn_wbc

for(id in samples){
  for(val in seq(0,20,by = 0.1)){
    cn <- seg[,id]
    med <- median(cn)
    cn[cn > med*10] <- med*10
    diff_ctc <- cn_ctc - cn*val
    diff_cdp <- cn_cdp - cn*val
    diff_wbc <- cn_wbc - cn*val
    
    vec_ctc <- c(vec_ctc, sum(abs(diff_ctc)))
    vec_cdp <- c(vec_cdp, sum(abs(diff_cdp)))
    vec_wbc <- c(vec_wbc, sum(abs(diff_wbc)))
    
    vec_val <- c(vec_val, val)
    vec_id <- c(vec_id, id)
    
  }
}

df_diff <- data.frame(
  id = vec_id, val = vec_val,
  diff_ctc = vec_ctc,
  diff_wbc = vec_wbc,
  diff_cdp = vec_cdp)

min_ctc_val <- by(df_diff, df_diff$id, function(x){ x[which.min(x[,'diff_ctc']),'val'] })
min_wbc_val <- by(df_diff, df_diff$id, function(x){ x[which.min(x[,'diff_wbc']),'val'] })
min_cdp_val <- by(df_diff, df_diff$id, function(x){ x[which.min(x[,'diff_cdp']),'val'] })

df_min <- data.frame(ctc = min_ctc_val, wbc = min_wbc_val, cdp = min_cdp_val)

cell_names <- samples

# Loop over cell names
for (cell_name in cell_names) {
  scaling_ctc <- df_min %>%
    filter(rownames(df_min) == cell_name) %>%
    select(ctc) %>%
    as.numeric()
  
  scaling_cdp <- df_min %>%
    filter(rownames(df_min) == cell_name) %>%
    select(cdp) %>%
    as.numeric()
  
  scaling_wbc <- df_min %>%
    filter(rownames(df_min) == cell_name) %>%
    select(wbc) %>%
    as.numeric()
  
  seg_filtered <- seg %>%
    select(CHR, START, END, all_of(cell_name))
  
  seg_filtered$scaled_ctc <- seg_filtered[[cell_name]] * scaling_ctc
  seg_filtered$scaled_cdp <- seg_filtered[[cell_name]] * scaling_cdp
  seg_filtered$scaled_wbc <- seg_filtered[[cell_name]] * scaling_wbc
  
  dots_data_ctc <- data.frame(
    x = 1:nrow(seg_filtered),
    y = seg_filtered$scaled_ctc
  )
  dots_data_ctc$y <- ifelse(dots_data_ctc$y > 10, 10, dots_data_ctc$y)
  
  dots_data_cdp <- data.frame(
    x = 1:nrow(seg_filtered),
    y = seg_filtered$scaled_cdp
  )
  dots_data_cdp$y <- ifelse(dots_data_cdp$y > 10, 10, dots_data_cdp$y)
  
  dots_data_wbc <- data.frame(
    x = 1:nrow(seg_filtered),
    y = seg_filtered$scaled_wbc
  )
  dots_data_wbc$y <- ifelse(dots_data_wbc$y > 10, 10, dots_data_wbc$y)
  
  ave$cn_ctc <- ifelse(ave$cn_ctc > 10, 10, ave$cn_ctc)
  ave$cn_cdp <- ifelse(ave$cn_cdp > 10, 10, ave$cn_cdp)
  ave$cn_wbc <- ifelse(ave$cn_wbc > 10, 10, ave$cn_wbc)
  
  # Generate plot
  plot <- ggplot(ave, aes(x = i)) +
    geom_line(aes(y = cn_ctc, color = "cn_ctc"), linewidth = 0.5) +
    geom_line(aes(y = cn_cdp, color = "cn_cdp"), linewidth = 0.5) +
    geom_line(aes(y = cn_wbc, color = "cn_wbc"), linewidth = 0.5) +
    geom_point(data = dots_data_wbc, aes(x = x, y = y), color = "blue", size = 0.01) +
    labs(
      y = "CN",
      color = "Legend"
    ) +
    scale_color_manual(values = c("cn_ctc" = "black", "cn_cdp" = "red", "cn_wbc" = "blue")) +
    theme_minimal() +
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5)
    )
  
  # Save plot as PDF
  ggsave(
    filename = file.path(output_dir, paste0("plot_", cell_name, "_WBC.pdf")),
    plot = plot,
    width = 12,
    height = 6,
    units = "in"
  )
}