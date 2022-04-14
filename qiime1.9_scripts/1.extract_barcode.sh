#!/bin/sh

## script that removes barcodes from already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## by default the first 8 bases are trimmed: otherwise, pass a manifest/mapping file with the barcode sequences

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="data/220225_M04028_0144_000000000-K6CMG"
output_dir="Analysis/milkqua_stools/qiime1.9/extract_barcode"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
temp_folder="temp/temp_fastq"

cd $currpath
echo "project folder is $project_home"

## copying files of interest to process
## e.g. K6CMG samples 45 - 76 is stools

echo " - creating temporary folder for sequence data ${temp_folder}"
if [ ! -d "${temp_folder}" ]; then
        mkdir -p ${temp_folder}
fi

echo " - copying relevant fastq files"
for i in {45..55}; 
do
	cp ${data_folder}/${i}*.fastq.gz ${project_home}/${temp_folder}
done


## make folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
fi

## using the Singularity container
echo " - calling the singularity container"
singularity run ${sing_container} multiple_extract_barcodes.py --input_dir=${temp_folder} --output_dir=${output_dir} --read1_indicator _R1 --read2_indicator _R2


echo "DONE!"

