## setting the enviornmnent
currpath=$(pwd)
datapath="data/210902_M04028_0139_000000000-JRGYP"
core=8

module load python3/intel/2020
source activate milkqua

export PATH=/home/users/filippo.biscarini.est/.local/bin:$PATH

if [ ! -d "$currpath/Analysis/raw_quality" ]; then
	mkdir $currpath/Analysis/raw_quality
fi

## Create analysis folders

echo $HOME
echo $currpath

## FastQC
cd $currpath

$HOME/software/FastQC/fastqc ${currpath}/${datapath}/*.fastq.gz -o $currpath/Analysis/raw_quality -t 8
cd $currpath/Analysis/raw_quality
multiqc .

echo "DONE!"
