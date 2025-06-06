#/bin/bash
#Run on ERIStwo, DO NOT NEED TO ACTIVATE CONDA ENV
###########################
#using bwa to align

INPUT_DIR=$1
export PATH=/data/gulhanlab/data/HybridCTC/fromMin/00.software/TrimGalore-0.6.10:$PATH
module load parallel
module load R
module load samtools/1.11
module load bwa/0.7.17
module load cutadapt/1.18-foss-2018b-Python-3.6.6
module load bedtools2/2.31.1
module load picard/2.24.1
export PICARD_HOME=/apps/released/java-toolchain/java-1.8/picard/build/libs

# Activate Conda environment for bamCoverage
# source /data/gulhanlab/ankit/miniconda3/etc/profile.d/conda.sh
# conda activate /data/gulhanlab/ankit/msk_frags/.conda

if [ $2 == "hg38" ]; then
    REF_FILE=/data/gulhanlab/data/HybridCTC/fromMin/01.ref/Human/bowtie2_index/hg38.fa
    genomesize=2913022398
elif [ $2 == "hg19" ]; then
    REF_FILE=/data/gulhanlab/data/HybridCTC/fromMin/01.ref/Human/hg19_bowtie2/hg19.fa
    genomesize=2864785220
elif [ $2 == "mm10" ]; then
    REF_FILE=/data/gulhanlab/data/HybridCTC/fromMin/01.ref/Mouse/mm10_bwa/mm10_ucsc.fa
    genomesize=2652783500
elif [ $2 == "mm10_EGFP_mCherry" ]; then
    REF_FILE=/data/gulhanlab/references/Mouse/mm10_EGFP_mCherry_bwa/mm10_EGFP_mCherry.fa
    genomesize=2652783500
else
    echo "Please enter genome reference: hg38, hg19 or mm10"
    exit
fi

#1. Trim data
cat ${INPUT_DIR}/FilePrefix_Sample12_1.txt | parallel -j 5 --xapply \
   trim_galore --paired --fastqc --quality 20 --length 20 -j 8 -o ${INPUT_DIR}/02.clean_data ${INPUT_DIR}/01.RawData/{}/{}_CKDL250010596-1A_232WNKLT3_L6_1.fq.gz ${INPUT_DIR}/01.RawData/{}/{}_CKDL250010596-1A_232WNKLT3_L6_2.fq.gz


#2. Alignment

for i in `cat ${INPUT_DIR}/FilePrefix_Sample12_1.txt`
do
   bwa mem -t 20 $REF_FILE ${INPUT_DIR}/02.clean_data/${i}_CKDL250010596-1A_232WNKLT3_L6_1_val_1.fq.gz ${INPUT_DIR}/02.clean_data/${i}_CKDL250010596-1A_232WNKLT3_L6_2_val_2.fq.gz | \
       samtools view -q 20 -@ 20 -o ${INPUT_DIR}/03.align/${i}.bam
   samtools sort -@ 20 -O bam -o ${INPUT_DIR}/03.align/${i}.sorted.bam ${INPUT_DIR}/03.align/${i}.bam
   samtools index -@ 20 -b $INPUT_DIR/03.align/${i}.sorted.bam
   samtools flagstat $INPUT_DIR/03.align/${i}.sorted.bam > $INPUT_DIR/03.align/${i}.sorted.bam.stat
done

# for i in `cat ${INPUT_DIR}/FilePrefix_Sample12_2.txt`
# do
#    bwa mem -t 20 $REF_FILE ${INPUT_DIR}/02.clean_data/${i}_CKDL250010598-1A_232WNKLT3_L5_1_val_1.fq.gz ${INPUT_DIR}/02.clean_data/${i}_CKDL250010598-1A_232WNKLT3_L5_2_val_2.fq.gz | \
#        samtools view -q 20 -@ 20 -o ${INPUT_DIR}/03.align/${i}.bam
#    samtools sort -@ 20 -O bam -o ${INPUT_DIR}/03.align/${i}.sorted.bam ${INPUT_DIR}/03.align/${i}.bam
#    samtools index -@ 20 -b $INPUT_DIR/03.align/${i}.sorted.bam
#    samtools flagstat $INPUT_DIR/03.align/${i}.sorted.bam > $INPUT_DIR/03.align/${i}.sorted.bam.stat
# done

#3.Remove duplicates
cat ${INPUT_DIR}/FilePrefix.txt | parallel -j 5 --xapply \
   java -jar $PICARD_HOME/picard.jar MarkDuplicates REMOVE_DUPLICATES=TRUE I=$INPUT_DIR/03.align/{}.sorted.bam O=$INPUT_DIR/04.rmdup/{}.rmdup.bam M=$INPUT_DIR/04.rmdup/{}.rmdup.log TMP_DIR=$INPUT_DIR/tmp_files

for i in `cat ${INPUT_DIR}/FilePrefix.txt`
do
   samtools sort -O bam -@ 20 -o $INPUT_DIR/04.rmdup/${i}.rmdup.sorted.bam $INPUT_DIR/04.rmdup/${i}.rmdup.bam
   samtools index -@ 20 -b $INPUT_DIR/04.rmdup/${i}.rmdup.sorted.bam
   samtools flagstat $INPUT_DIR/04.rmdup/${i}.rmdup.sorted.bam > $INPUT_DIR/04.rmdup/${i}.rmdup.sorted.bam.stat
done

#4.Generate bw files
# cat ${INPUT_DIR}/FilePrefix.txt | parallel -j 5 --xapply \
#    bamCoverage --bam $INPUT_DIR/03.rmdup/{}.rmdup.sorted.bam -o $INPUT_DIR/04.bigwig/{}.100bp.rpkm.bw -p 20 --binSize 100 --normalizeUsing RPKM --effectiveGenomeSize ${genomesize}

#5.BamToBed
for i in `cat ${INPUT_DIR}/FilePrefix.txt`
do
   bamToBed -i ${INPUT_DIR}/04.rmdup/${i}.rmdup.sorted.bam > ${INPUT_DIR}/06.bed/${i}.bed
   samtools coverage -o ${INPUT_DIR}/07.coverage/${i}.coverage ${INPUT_DIR}/04.rmdup/${i}.rmdup.sorted.bam
done
gzip ${INPUT_DIR}/06.bed/*