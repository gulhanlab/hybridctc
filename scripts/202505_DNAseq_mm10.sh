#!/bin/bash
#SBATCH -J 202505_Sample12_scDNA
#SBATCH -n 8
#SBATCH -t 5-00:00:00
#SBATCH -q normal
#SBATCH -o /data/gulhanlab/ankit/HybridCTC/Sample12/slurm_output/202504_Sample12_scDNA_mm10.%J.out
#SBATCH -e /data/gulhanlab/ankit/HybridCTC/Sample12/slurm_output/202504_Sample12_scDNA_mm10.%J.err
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=asingh46@mgh.harvard.edu

bash /data/gulhanlab/data/HybridCTC/fromMin/02.scripts/scDNA_CNV_Sample12.sh /data/gulhanlab/data/HybridCTC/Sample12 mm10_EGFP_mCherry