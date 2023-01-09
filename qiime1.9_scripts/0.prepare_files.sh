#!/bin/sh

## script that prepare the data folder for multiqc and extract barcode etc.
## sample data may be a subset of one sequencing run or come from multiple sequencing runs (different folders)
## the file of metadata (mapping file) must also be changed accordingly (r scripts for this)
## !! beware of sample names: if they overlap (e.g. 1-32 10-44 from different runs), 
## they must be renamed (both the fastq file names and the sample names in the metadata file)

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder1="${project_home}/data/mock_communities"
data_folder2=""
output_folder="Analysis/mock_communities/qiime1.9"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
temp_folder="temp"
sample_start1=1 #first sample to use (in the sequence)
sample_end1=3 #last sample to use (in the sequence)
#sample_start2=10 #first sample to use (in the sequence)
#sample_end2=12 #last sample to use (in the sequence)
prefix="mock" # prefix to remove from sample file names (if any: !! usually this is left empty !!)

echo "current folder $currpath"
cd $project_home

## copying files of interest to process
## e.g. K6CMG samples 45 - 76 is stools

echo " - creating temporary folder for sequence data ${temp_folder}"
if [ ! -d "${output_folder}/${temp_folder}" ]; then
        mkdir -p ${output_folder}/${temp_folder}
	chmod g+rwx ${output_folder}/${temp_folder} 
fi

## clean the temporary folder (in case there were previous left overs)
echo " - cleaning the temp/temp_fastq folder"
cd ${output_folder}/${temp_folder}
rm -f *
cd $project_home

echo " - copying relevant fastq files from data folder 1"
for i in $(seq ${sample_start1} ${sample_end1}); 
do
	echo "file ${data_folder1}/${i}_"
	cp ${data_folder1}/${prefix}_${i}_*.fastq.gz ${output_folder}/${temp_folder}
done

if [ ! ${data_folder2} == "" ]; then
	echo " - copying relevant fastq files from data folder 2"
	for i in $(seq ${sample_start2} ${sample_end2});
	do
		echo "file ${data_folder2}/${i}_"
        	cp ${data_folder2}/${i}_*.fastq.gz ${output_folder}/${temp_folder}
	done
fi

if [ ! $prefix == "" ]; then
	cd ${output_folder}/${temp_folder}
	for i in  *"${prefix}"*.gz;
	do 
		echo "renaming file $i"
		mv "$i" "${i#*_}";
	done
fi

echo "DONE!"

