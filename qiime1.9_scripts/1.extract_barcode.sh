#!/bin/sh

## script that removes barcodes from already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## by default the first 8 bases are trimmed: otherwise, pass a manifest/mapping file with the barcode sequences

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="data/210902_M04028_0139_000000000-JRGYP_milk"
output_dir="Analysis/extract_barcode"
sing_container="${project_home}/Qiime1.9.sif"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
fi

## using the Singularity container
echo " - calling the singularity container"
singularity run ${sing_container} multiple_extract_barcodes.py --input_dir=${data_folder} --output_dir=${output_dir} --read1_indicator _R1 --read2_indicator _R2


echo "DONE!"

