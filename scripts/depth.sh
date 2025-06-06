#!/bin/bash
#SBATCH -J average_depth
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -p short
#SBATCH -o /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/average_depth.%J.out
#SBATCH -e /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/average_depth.%J.err
#SBATCH --mem-per-cpu=4G

# Load required modules
module load gcc samtools

# Output file
OUTPUT_FILE="/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2/average_depth.txt"
bam_path="/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2/04.rmdup"

# Write header
echo -e "Sample_Name\tAverage_Depth" > "$OUTPUT_FILE"

# Compute average depth for each BAM file
for bam_file in ${bam_path}/*.rmdup.sorted.bam; do
    sample_name=$(basename "$bam_file")  # Extract sample name
    genome_size=$(samtools view -H "$bam_file" | awk '$1 == "@SQ" {sum += substr($3, 4)} END {print sum}')
    
    avg_depth=$(samtools depth "$bam_file" | awk -v genome_size="$genome_size" '{ sum += $3 } END { print sum / genome_size }')
    
    echo -e "$sample_name\t$avg_depth" >> "$OUTPUT_FILE"
done
