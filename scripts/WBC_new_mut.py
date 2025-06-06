import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pysam
import os

# Input folder containing BAM files
bam_folder = '/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2/04.rmdup'
exclusive_wbc_df_homozygous_filtered_snvs = pd.read_csv('/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/HybridCTC/exclusive_wbc_df_homozygous_filtered_snvs_new.csv')

# Output folder for results
output_folder = '/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/output_sample12/'

# Loop through all BAM files in the folder
for bam_file in os.listdir(bam_folder):
    if bam_file.endswith('.rmdup.sorted.bam'):  # Filter for BAM files
        bam_path = os.path.join(bam_folder, bam_file)
        bam = pysam.AlignmentFile(bam_path, 'rb')

        results = []

        # Loop through the regions from the filtered DataFrame
        for idx, row in exclusive_wbc_df_homozygous_filtered_snvs.iterrows():
            chromosome = row['CHROM']
            start = row['POS'] - 1  # Adjusting for 0-based indexing in BAM
            end = row['POS']
            alt_base = row['ALT']

            # Initialize counters
            alt_count = 0
            total_count = 0
            base_counts = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

            # Fetch reads from the BAM file
            for read in bam.fetch(chromosome, start, end):
                if not read.is_unmapped:
                    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                        if ref_pos == start:  # Match the position
                            if query_pos is not None:
                                read_base = read.query_sequence[query_pos].upper()
                                if read_base in base_counts:
                                    base_counts[read_base] += 1  # Increment base count
                                total_count += 1
                                if read_base == alt_base.upper():
                                    alt_count += 1

            # Calculate fraction
            fraction = round(alt_count / total_count, 4) if total_count > 0 else 0

            # Append results
            results.append({
                'Chromosome': chromosome,
                'Position': row['POS'],
                'Alt_Base': alt_base,
                'Alt_Read_Support': alt_count,
                'Total_Reads': total_count,
                'Fraction': fraction,
                'A_Count': base_counts['A'],
                'C_Count': base_counts['C'],
                'T_Count': base_counts['T'],
                'G_Count': base_counts['G']
            })

        # Close the BAM file
        bam.close()

        # Convert results to a DataFrame and save as CSV
        output_df = pd.DataFrame(results)
        output_file = os.path.join(output_folder, f"WBC_alt_read_support_{os.path.splitext(bam_file)[0]}.csv")
        output_df.to_csv(output_file, index=False)
        print(f"Processed {bam_file} and saved results to {output_file}")