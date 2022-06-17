## setting the enviornmnent
currpath=$(pwd)
datapath="data/mock-sequences"
targetdir="Analysis/mock_communities/0.fastqc"
core=8

## load the conda env with multiqc installed
module load python3/intel/2020
#source activate milkqua
source activate /gpfs/home/projects/MILKQUA/Conda/env-arriba

## MULTIQC can also be invoked from here (for some users who have multiqc installed in their locale)
## Exporting the path to bin by the way, also makes available additional software which may be useful
#export PATH=/home/users/luiz.dematos/.local/bin:$PATH
export PATH=/home/users/filippo.biscarini.est/.local/bin:$PATH

## Create analysis folders
echo "current path: $currpath"

targetpath="$currpath/$targetdir"
echo "target path: $targetpath"

if [ ! -d $targetpath ]; then
	mkdir -p $targetpath
	chmod g+rxw $targetpath
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
