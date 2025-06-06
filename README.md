# HybridCTC
Scripts for processing and analyzing DNA for HybridCTC cells

## Directory Setup
1.	[ERIS] Create directory for the current sample. I create folders in ‘/data/gulhanlab/data/HybridCTC’.
2.	[ERIS] Create required folders in the current sample’s directory to store outputs from different steps: 
```mkdir -p results 02.clean_data 03.align 04.rmdup 05.bigwig 06.bed 07.coverage```
3.	[ERIS] Create symlink for raw fastqs in current sample directory’s 01.RawData folder: navigate to the ‘01.RawData’ for current sample: 
```ln -s /data/rama/labMembers/ea865/HybridCells/23April2025_B16F10_Round3_Lung_P1P2/Round3_Lung_P1P2_DNA_MALBAC/01.RawData 01.RawData```
4.	[ERIS] Create the FilePrefix.txt file with all cell names in the 01.RawData folder: navigate to the ‘01.RawData’ for current sample:  
```ls -d * > /data/gulhanlab/data/HybridCTC/Sample12/FilePrefix.txt```


## Processing
### Now we are ready to begin processing. Refer to ‘scDNA_CNV.sh’ script for code. Below is description for each step:
1.	[ERIS] Trim Data: Edit the ‘Trim data’ (step1) step’s input file naming to match the cell naming convention obtained from sequencing. Following this, you can submit the job for this step using SLURM (refer 202504_DNAseq_mm10.sh). Outputs are stored in ‘02.clean_data’ folder that was set-up earlier in the current sample’s dir.
2.	[ERIS] Alignment: Following trimming, the script runs alignment using bwa-mem. Again, edit the input file names in the script following the naming convention of the cells. Outputs are stored in ‘03.align’ folder.
3.	[ERIS] Remove duplicates: After alignment, we remove duplicates using Picard library’s ‘MarkDuplicates’ command. Outputs are stored in ‘04.rmdup’ folder.
4.	[ERIS] BAM to Bed: We convert BAMs generated in previous step to bed file for running ginkgo CNV caller. Bedtools ‘bamToBed’ is used to accomplish this task. 
5.	[ERIS] Once bed files are generated, transfer the de-duplicated BAMs and bed file to 02 for further processing. Need to generate id_rsa keys on 02 and save it in ERIS to enable file transfer from ERIS to 02: 
```scp -r ./04.rmdup ./06.bed ans4371@transfer.rc.hms.harvard.edu:/n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2``` 
Approve Duo prompt and transfer should begin. I create an interactive job and run this transfer in that interactive job. 
6.	[02] Ginkgo: Navigate to dir: 
```/n/data1/hms/dbmi/gulhan/lab/ankit/scripts``` 
on 02. Run ginkgo using this command (example command for Ezgi’s Sample12): 
```./ginkgo_script.sh /n/data1/hms/dbmi/gulhan/lab/DATA/CTC/HybridCTC/23April2025_B16F10_Round3_Lung_P1P2``` 
as an interactive job. The script uses bed files to run ginkgo. The parameter in the script should be path to the folder with sub-folder containing the bed files. Usually requires 3GB memory and 4-5 hrs to process. Request resources for interactive job accordingly.
7.	[02] Avg depth and fraction of genome covered calculations: Navigate to dir: 
```/n/data1/hms/dbmi/gulhan/lab/ankit/scripts``` 
on 02. Modify the ‘depth.sh’ & ‘coverage.sh’ scripts to point to the correct input and output directories. Submit jobs to calculate depth and coverage using ‘sbatch depth.sh’ & ‘sbatch coverage.sh’ respectively. 
8.	[02] Calculate mutation support: Navigate to dir: 
```/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/mutation_calling/HybridCTC``` 
Run python scripts ‘CTC_new_mut.py’ & ‘WBC_new_mut.py’ to calculate CTC and WBC mutation support respectively. This will generate 2 files for each cell (CTC & WBC) with columns chr, position, Alt_Base,	Alt_Read_Support, Total_Reads, Fraction, A_Count, C_Count, T_Count, G_Count. 
9.	[Local] Calculate allele support and mutation fraction: Use R-script ‘muts.Rmd’ 1st code block to calculate CTC & WBC allele mutation fraction support. Calculate allele support raw numbers for CTC & WBC using 2nd  code block of ‘muts.Rmd’.
10.	[Local] Scaled Plots: Plot scaled plots using ‘scale_blood.R’ script for CDP, CTC, and WBC for each cell. Fill metadata table with sample names before running this so the script knows the cells to plot and these cell names should match the names in SegNorm file obtained as output from the ginkgo run.

## Analysis
Use the data and plots obtained in the processing step to populate the metadata table. 
1.	Use ginkgo’s CNV plots, scaled plots and mutation fraction numbers to derive consensus cell annotations.
2.	Discard cells with low coverage and depth as bad quality cells.

## Germline Plots
Using consensus cell annotations obtained from the metadata table, plot germline plots for the sample using ‘mut_plots_sample.R’ script. Integrate with plots for overall cells from previous samples (blood and/or lung) using ‘mut_plots.R’. Note that both these scripts ignore the ‘BadQual’ cells while plotting.

## Scaled SegNorm
Using consensus cell annotations obtained from the metadata table, scaled copy number values from the original SegNorm file obtained from ginkgo run using ‘scale_SegNorm.R’. 
