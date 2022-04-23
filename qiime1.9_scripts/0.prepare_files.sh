#!/bin/sh

## script that prepare the data folder for multiqc and extract barcode etc.
## sample data may be a subset of one sequencing run or come from multiple sequencing runs (different folders)
## the file of metadata (mapping file) must also be changed accordingly (r scripts for this)
## !! beware of sample names: if they overlap (e.g. 1-32 10-44 from different runs), 
## they must be renamed (both the fastq file names and the sample names in the metadata file)

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder1="data/JRGYP"
data_folder2=""
output_dir="temp/temp_fastq"
#sing_container="${project_home}/Qiime1.9.sif"
#sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
temp_folder="temp/temp_fastq"
sample_start1=1 #first sample to use (in the sequence)
sample_end1=6 #last sample to use (in the sequence)
sample_start2=10 #first sample to use (in the sequence)
sample_end2=12 #last sample to use (in the sequence)

cd $currpath
echo "project folder is $project_home"

## copying files of interest to process
## e.g. K6CMG samples 45 - 76 is stools

echo " - creating temporary folder for sequence data ${temp_folder}"
if [ ! -d "${temp_folder}" ]; then
        mkdir -p ${temp_folder}
	chmod g+rwx ${temp_folder} 
fi

## clean the temporary folder (in case there were previous left overs)
echo " - cleaning the temp/temp_fastq folder"
cd ${temp_folder}
rm -f *
cd $currpath

echo " - copying relevant fastq files from data folder 1"
for i in $(seq ${sample_start1} ${sample_end1}); 
do
	echo "file ${data_folder1}/${i}_"
	cp ${data_folder1}/${i}_*.fastq.gz ${project_home}/${temp_folder}
done

if [ ! ${data_folder2} == "" ]; then
	echo " - copying relevant fastq files from data folder 2"
	for i in $(seq ${sample_start2} ${sample_end2});
	do
		echo "file ${data_folder2}/${i}_"
        	cp ${data_folder2}/${i}_*.fastq.gz ${project_home}/${temp_folder}
	done
fi

echo "DONE!"

