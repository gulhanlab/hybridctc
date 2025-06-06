##
##  Min Yang's comments:
##  This script was adapted from Ben and I will use my own environment
##  /data/rama/labMembers/my1032/00.software/Anaconda3/envs/r_env
##  I add the choose_genome parameter in analyze.sh script

##  ginkgo CNV.
##  Based on ccd/projects/2017.guo/src/guo306.sh, which
##  stopped being runable when module load R/3.6.3 was taken away by ERIS due
##  to security concerns.
##
##  Run this script on grx23.partners.org, which is an ERIS2 machine.
##
##  Before running this script, put a non-commented-out version of the
##  following in your ~/.bashrc file and then start a new terminal:
##
##    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/data/rama/labMembers/my1032/00.software/Anaconda3/lib
##
##  Before running this script, execute the following command:
##    conda activate /data/rama/labMembers/my1032/00.software/Anaconda3/envs/r_env
##
##  ATTENTION: Using full path for "bedDir"

# bedDir=/data/rama/ccd/projects/2022.shihbo/output/shihbo288/data/12_ginkgo
# size=5000000
# runName=shihbo288.${size}
# oDir=/data/rama/ccd/projects/2022.shihbo/output/shihbo307/shihbo288/${size}
# aligner=bwa
# readLength=101          # default for web site.

bedDir=${1}
runName=${2}
oDir=${3}
size=${4}               # available values are 10000000, 5000000, 2500000,
                        # 1000000, 500000, 250000, 175000, 100000, 50000,
                        # 25000, 10000
aligner=${5}            # bowtie or bwa
readLength=${6}         # available values are: 150, 101, 76, 48
choose_genome=${7}      # available values are:hg19, mm10
gDir=/n/data1/hms/dbmi/gulhan/lab/software/ginkgo        #Change ginkgo folder location
uDir=$gDir/uploads/$runName

mkdir -p $uDir
mkdir -p $oDir

cd $uDir
echo $bedDir

for bedFile in $bedDir/*.bed.gz
do
    ln -s $bedFile
done

/bin/ls *.bed.gz > list
sed "s/500000_101/${size}_${readLength}/" \
    < $gDir/config.fromWebsite > $uDir/config

if [ $aligner = bwa ]
then
    sed "s/_bowtie/_bwa/" < $uDir/config > $uDir/config.temp
    mv $uDir/config.temp $uDir/config
fi

cd $gDir
./scripts/analyze_CTC.sh $runName ${choose_genome}
mkdir -p $oDir

cd $uDir
rm *.bed.gz
cp * $oDir

cd
rm -r $uDir

touch $oDir/$runName.done