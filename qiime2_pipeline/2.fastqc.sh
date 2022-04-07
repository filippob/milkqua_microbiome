#!/bin/bash


#setting the environmnent
currpath=$(pwd)
project_home="$HOME/MILKQUA"
datapath="data/220225_M04028_0144_000000000-K6CMG_faeces"
outdir="Analysis/Faeces_analysis/0.fastqc"
core=8 

## Create analysis folders

echo $HOME
echo $currpath
export MPLCONFIGDIR=$project_home/R_MYPACKAGES

##launch conda env to allow multiqc job
module load python3/intel/2020
source activate /gpfs/home/projects/MILKQUA/Conda/env-arriba

## FastQC
cd $currpath

if [ ! -d $outdir ]; then
	mkdir -p $outdir
fi

/gpfs/software/Fastqc/FastQC-0.11.9/fastqc ${datapath}/*.fastq.gz -o ${outdir} -t 8
cd ${outdir}
multiqc .

echo "DONE!"
