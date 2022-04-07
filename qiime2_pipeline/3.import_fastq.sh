#!/bin/sh

## script that imports into Qiime2 already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## sequence still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)
## type:
## input-format:

## Importing data with different formats

project_home="$HOME/MILKQUA"
data_folder="data/220225_M04028_0144_000000000-K6CMG_faeces"
manifest="Config/manifest_faeces2.csv"
output_dir="Analysis/milkqua_stools" ##<project-name>_<sample-type>
qiime2_cont="/gpfs/software/Container/qiime2_2022.2.sif"
name="milkqua_stools"

##create import directory

if [ ! -d "${project_home}/${output_dir}/import2" ]; then
        mkdir -p ${project_home}/${output_dir}/import2
fi

echo "project folder is $project_home"
echo "output folder is $output_dir"

## importing data
echo " - running qiime to import data"

## Paired-end reads with quality (fastq files)

singularity run\
	 --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES\
	$qiime2_cont qiime tools import\
	--type 'SampleData[PairedEndSequencesWithQuality]'\
	--input-path $project_home/$manifest \
	--output-path $project_home/$output_dir/import2/${name}2.qza\
	--input-format PairedEndFastqManifestPhred33V2 

echo "DONE!"
