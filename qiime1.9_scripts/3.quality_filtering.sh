#!/bin/sh

## script that applies quality filters and assembles paired reads into a single Fasta file
## input data format is fastq files from 2.join_paired_ends.sh
## 

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
input_folder="Analysis/milkqua_skinswab/qiime1.9/join_paired_ends"
output_dir="Analysis/milkqua_skinswab/qiime1.9/filtered_fasta"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
paramfile="Config/filter.parameters"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "${output_dir}" ]; then
	mkdir -p ${output_dir}
	chmod g+rxw ${output_dir}
fi

## generating qiime_parameters file
#echo "split_libraries_fastq:max_bad_run_length 3 split_libraries_fastq:min_per_read_length_fraction 0.75 split_libraries_fastq:sequence_max_n 0 split_libraries_fastq:phred_quality_threshold 19" >| "$PWD/qiime_parameters_quality"

## using the Singularity container
echo " - calling the singularity container for quality filtering"
singularity run ${sing_container} multiple_split_libraries_fastq.py -p "${project_home}/$paramfile" --input_dir=$input_folder  --output_dir=$output_dir

echo "DONE!"

