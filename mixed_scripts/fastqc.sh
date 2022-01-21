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

fullpath="$currpath/$targetdir"
echo "target path: $fullpath"

if [ ! -d $fullpath ]; then
	mkdir $fullpath
fi


## FastQC
cd $currpath

# $HOME/software/FastQC/fastqc ${currpath}/${datapath}/*.fastq.gz -o $currpath/Analysis/$targetdir -t 8
## load fastqc
module load genetics/fastqc
fastqc ${currpath}/${datapath}/*.fastq.gz -o $fullpath -t 8

## MultiQC
cd $fullpath
multiqc .

echo "DONE!"
