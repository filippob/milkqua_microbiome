## setting the enviornmnent
currpath=$(pwd)
datapath="data/210902_M04028_0139_000000000-JRGYP_swab"
targetdir= "Analysis/raw_quality/swab"
core=8

module load python3/intel/2020
source activate milkqua

export PATH=/home/users/filippo.biscarini.est/.local/bin:$PATH

if [ ! -d `$currentpath/Analysis/$targetdir`]; then
	mkdir $currpath/Analysis/$targetdir
fi

## Create analysis folders

echo $HOME
echo $currpath

## FastQC
cd $currpath

$HOME/software/FastQC/fastqc ${currpath}/${datapath}/*.fastq.gz -o $currpath/Analysis/$targetdir -t 8
cd $currpath/Analysis/$targetdir
multiqc .

echo "DONE!"
