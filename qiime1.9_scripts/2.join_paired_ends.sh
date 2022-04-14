#!/bin/sh

## script that joins paired end reads after the barcode has been removed
## input data format is fastq: R1 and R2 files (forward and reverse) from 1.extract_barcode.sh
## 

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
input_folder="Analysis/milkqua_stools/qiime1.9/extract_barcode"
output_dir="Analysis/milkqua_stools/qiime1.9/join_paired_ends"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
paramfile="Config/join.parameters"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
fi

## generating qiime_parameters file
#echo "join_paired_ends:pe_join_method SeqPrep" >| qiime_parameters_joining 

## using the Singularity container
echo " - calling the singularity container for read joining"
singularity run ${sing_container} multiple_join_paired_ends.py --input_dir=${input_folder} --output_dir=${output_dir} --include_input_dir_path --parameter_fp=$paramfile --read1_indicator _R1 --read2_indicator _R2

echo "DONE!"

