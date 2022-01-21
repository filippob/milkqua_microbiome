## setting the enviornmnent
currpath=$(pwd)
datapath="data/subset"
targetdir="Analysis/raw_quality/subset"
core=8

module load python3/intel/2020
source activate milkqua

export PATH=/home/users/filippo.biscarini.est/.local/bin:$PATH
## Create analysis folders

echo "current path: $currpath"

targetpath="$currpath/$targetdir"
echo "target path: $targetpath"

if [ ! -d $targetpath ]; then
	mkdir $targetpath
fi


## FastQC
cd $currpath

# $HOME/software/FastQC/fastqc ${currpath}/${datapath}/*.fastq.gz -o $currpath/Analysis/$targetdir -t 8
## load fastqc
module load genetics/fastqc
fastqc ${currpath}/${datapath}/*.fastq.gz -o $targetpath -t 8

## MultiQC
cd $targetpath
multiqc .

echo "DONE!"
