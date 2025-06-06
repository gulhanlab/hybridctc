#!/bin/bash
#SBATCH -J genome_fraction
#SBATCH -n 1
#SBATCH -t 05:00:00
#SBATCH -p short
#SBATCH -o /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/genome_fraction.%J.out
#SBATCH -e /n/data1/hms/dbmi/gulhan/lab/ankit/slurm_output/genome_fraction.%J.err
#SBATCH --mem-per-cpu=4G

# Load required modules
module load gcc samtools

# Output file
OUTPUT_FILE="/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2/fraction_genome_covered.txt"
bam_path="/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2/04.rmdup"

# Write header
echo -e "Sample_Name\tFraction_Genome_Covered" > "$OUTPUT_FILE"

# Compute fraction of genome covered for each BAM file
for bam_file in ${bam_path}/*.rmdup.sorted.bam; do
    sample_name=$(basename "$bam_file")  # Extract sample name
    genome_size=$(samtools view -H "$bam_file" | awk '$1 == "@SQ" {sum += substr($3, 4)} END {print sum}')
    
    fraction_covered=$(samtools depth "$bam_file" | awk -v genome_size="$genome_size" '{ if ($3 > 0) covered++ } END { print covered / genome_size }')
    
    echo -e "$sample_name\t$fraction_covered" >> "$OUTPUT_FILE"
done