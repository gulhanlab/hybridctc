#/bin/bash

INPUT_DIR=$1

#Run Ginkgo in local
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/n/data1/hms/dbmi/gulhan/lab/software/Anaconda3/lib
module load gcc R
ginkgo_script=/n/data1/hms/dbmi/gulhan/lab/ankit/scripts/ginkgo_CTC.sh
readLength=48 # This can be 48, 76, 101, 150
#size value can be avalible with 10000000, 5000000, 2500000, 1000000, 500000, 250000, 175000, 100000, 50000, 25000, 10000

for size in 500000 5000000
do
    mkdir -p ${INPUT_DIR}/07.ginkgo/bin_${size}.readLength_${readLength}
    $ginkgo_script ${INPUT_DIR}/06.bed ginkgo.${size}.${readLength} ${INPUT_DIR}/07.ginkgo/bin_${size}.readLength_${readLength} $size bwa $readLength mm10
done