#!/bin/bash

## script to run fastqc and multiqc on a subset of the samples
## the subset is defined by a numeric range via seq start end

## parameters
project_home=$HOME/MILKQUA
input_dir="220225_M04028_0144_000000000-K6CMG"
targetdir="Analysis/milkqua_stools/0.fastqc"
multiqc_sif="singularity_containers/multiqc-1.9.sif"
nstart=45  #number of first sample to include in the subset
nend=76    #number of last sample to include in the subset

cd $project_home

## create temporary folder
if [ ! -d ${project_home}/data/temp ]; then
	
	echo " - creating temporary folder"
	mkdir ${project_home}/data/temp

	## copying files of interest in the temporary folder
	cd ${project_home}/data
	## !! SPECIFY HERE THE SAMPLES TO SUBSET !!
	for i in `seq $nstart $nend`; do
        	echo "copying file ${input_dir}/${i}*.fastq.gz"
        	cp ${input_dir}/${i}*.fastq.gz temp/
	done;
	cd $project_home
fi

## FastQC

## load fastqc
#module load genetics/fastqc ## module with fastqc (but no multiqc)
## activate conda env with fastqc adn multiqc
module load python3/intel/2020
#source activate milkqua ## personal conda env
source activate /gpfs/home/projects/MILKQUA/Conda/env-arriba/

## running fastqc on single samples
if [ ! -d $targetdir ]; then
	echo " - creating trarget dir"
       	mkdir $targetdir
fi

if [ -z "$(ls -A $targetdir)" ]; then
	echo "(target folder is empty)"
	echo " - running fastqc on single samples"
        fastqc data/temp/*.fastq.gz -o $targetdir -t 8
        #singularity run ${project_home}/${multiqc_sif} fastqc data/temp/*.fastq.gz -o $targetdir -t 8
fi

## MultiQC
cd $targetdir
echo " - running multiqc on fastqc files"
#multiqc .
singularity run ${project_home}/${multiqc_sif} -f .

## removing temporary folder
cd $project_home
echo " - removing temporary folder"
rm -r data/temp

echo "done"



