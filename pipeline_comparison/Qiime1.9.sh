#!/bin/sh

## setting up the environment

currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder1="data/subset"
#data_folder2=""
output_dir="Analysis/milkqua_milk/qiime1.9"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
temp_folder="temp"
sample_start1=18 #first sample to use (in the sequence)
sample_end1=72 #last sample to use (in the sequence)
#sample_start2= #first sample to use (in the sequence)
#sample_end2= #last sample to use (in the sequence)
#prefix="" # prefix to remove from sample file names (if any: !! usually this is left empty !!)
joinparam="Config/join.parameters"
filterparam="Config/parameters/filter.parameters"
dbpath="Databases/SILVA_132_QIIME"
rpath="milkqua_microbiome/mixed_scripts"

cd $project_home
echo "we are in $currpath. Let's start by creating useful folders"

if [ ! -d "$output_dir/$temp_folder" ]; then
        mkdir -p $output_dir/$temp_folder
	chmod g+rwx $output_dir/$temp_folder
	cd $output_dir/$temp_folder
	rm *
	cd $project_home 
fi

if [ ! -d "${output_dir}/1.extract_barcode}" ]; then
	mkdir -p ${output_dir}/1.extract_barcode
	chmod g+rxw ${output_dir}/1.extract_barcode
fi

if [ ! -d "${output_dir}/2.join_paired_ends" ]; then
	mkdir -p ${output_dir}/2.join_paired_ends
	chmod g+rxw ${output_dir}/2.join_paired_ends
fi

if [ ! -d "${output_dir}/3.quality_filtering" ]; then
	mkdir -p ${output_dir}/3.quality_filtering
	chmod g+rxw ${output_dir}/3.quality_filtering
fi

if [ ! -d "${output_dir}/4.OTU_picking" ]; then
        mkdir -p ${output_dir}/4.OTU_picking
        chmod g+rxw ${output_dir}/4.OTU_picking
fi

if [ ! -d "${output_dir}/5.filter_OTUs" ]; then
        mkdir -p ${output_dir}/5.filter_OTUs
        chmod g+rxw ${output_dir}/5.filter_OTUs
fi

## script that prepare the data folder for multiqc and extract barcode etc. Sample data may be a subset of one sequencing run or come from multiple sequencing runs (different folders). The file of metadata (Mapping file) must also be changed accordingly (r scripts for this). 
## !! beware of sample names: if they overlap (e.g. 1-32 10-44 from different runs), they must be renamed (both the fastq file names and the sample names in the metadata file)

echo "we are in $currpath. Let's copy the files of interest to process"

cd $currpath

echo " - copying relevant fastq files from data folder 1"
for i in $(seq ${sample_start1} ${sample_end1}); 
do
	echo "file ${data_folder1}/${i}_"
	cp ${data_folder1}/${i}_*.fastq.gz $output_dir/${temp_folder}
done


##use the following line if you have to copy files from more than 1 folder
#echo " - copying relevant fastq files from data folder 1"

#if [ ! ${data_folder2} == "" ]; then
#	echo " - copying relevant fastq files from data folder 2"
#	for i in $(seq ${sample_start2} ${sample_end2});
#	do
#		echo "file ${data_folder2}/${i}_"
#        	cp ${data_folder2}/${i}_*.fastq.gz ${project_home}/${temp_folder}
#	done


echo "DONE!"

## script that removes barcodes from already demultiplexed data. Input data format is fastq: R1 and R2 files (forward and reverse). 
## By default the first 8 bases are trimmed: otherwise, pass a manifest/mapping file with the barcode sequences

cd $currpath
echo "we are in $currpath. Let's start barcodes removal"

## using the Singularity container
echo " - calling the singularity container"
singularity run ${sing_container} multiple_extract_barcodes.py \
	--input_dir=$output_dir/${temp_folder} \
	--output_dir=${output_dir}/1.extract_barcode \
	--read1_indicator _R1 --read2_indicator _R2

echo " - removing barcodes files"
cd $output_dir/1.extract_barcode
find . -name \*barcodes.fastq -type f -delete

cd $currpath

echo "DONE!"


## script that joins paired end reads after the barcode has been removed
## input data format is fastq: R1 and R2 files (forward and reverse) from 1.extract_barcode.sh

cd $currpath
echo "we are in $currpath. Let's start paired ends joining"

## generating qiime_parameters file
#cd $currpath/Config
#echo "join_paired_ends:pe_join_method SeqPrep" >| "join.parameters" 

cd $currpath 

echo " - calling the singularity container for read joining"
singularity run ${sing_container} multiple_join_paired_ends.py \
	--input_dir=${output_dir}/1.extract_barcode \
	--output_dir=${output_dir}/2.join_paired_ends \
	--include_input_dir_path \
	--parameter_fp=$joinparam \
	--read1_indicator=_R1 --read2_indicator=_R2

chmod g+rxw -R ${output_dir}/2.join_paired_ends

echo " - removing unassembled reads"

cd $output_dir/2.join_paired_ends
find . -name \*unassembled*.fastq.gz -type f -delete

# clean output folder from potential barcode subfolders

#echo " - cleaning the output folder (barcodes subfolders)"
#rm -r ${output_dir}/2.join_paired_ends/*barcodes

echo "DONE!"

echo "running the report file count_reads"

cd $currpath
sh $rpath/count_reads.sh

echo "done!"

## script that applies quality filters and assembles paired reads into a single Fasta file
## input data format is fastq files from 2.join_paired_ends.sh

cd $currpath
echo "we are in $currpath. Lets start quality filtering"


## generating qiime_parameters file

#cd $currpath/Config
#echo "split_libraries_fastq:max_bad_run_length 3 split_libraries_fastq:min_per_read_length_fraction 0.75 split_libraries_fastq:sequence_max_n 0 split_libraries_fastq:phred_quality_threshold 19" >| "filter.parameters"

cd $currpath

## using the Singularity container
echo " - calling the singularity container for quality filtering"

singularity run ${sing_container} multiple_split_libraries_fastq.py \
	-p $filterparam \
	-i $output_dir/2.join_paired_ends \
	--read_indicator=_R1 --include_input_dir_path \
	--output_dir=$output_dir/3.quality_filtering

chmod g+rxw -R $output_dir/3.quality_filtering

echo "DONE!"

echo "running the report file reads_after_filter"

cd $currpath
sh $rpath/reads_after_filter.sh

echo "done!"

#!/bin/bash
## Script that determines OTUS. OTUs were determined by aligning quality-filtered reads against the QIIME-compatible SILVA reference FASTA file, release 123, with minimum 97% clustering (https://www.arb-silva.de/download/archive/qiime/)

###NOTE: database can be changed

cd $currpath
echo "we are in $currpath. Lets start OTU picking with $dbpath"

## using the Singularity container

echo " - calling the singularity container for OTU picking"
singularity run ${sing_container} \
        pick_closed_reference_otus.py \
	--input_fp ${output_dir}/3.quality_filtering/seqs.fna \
        --output_dir ${output_dir}/4.OTU_picking \
	--reference_fp $dbpath/rep_set/rep_set_all/97/silva132_97.fna \
        --taxonomy_fp $dbpath/taxonomy/taxonomy_all/97/raw_taxonomy.txt \
        --force

chmod g+xrw -R $output_dir/4.OTU_picking

echo " - converting biom file to tsv file"

cd $output_dir/4.OTU_picking
singularity run ${sing_container} biom convert -i otu_table.biom -o otu_table.txt --to-tsv --header-key taxonomy

echo "DONE!"

###Script that performs filtering on by total count across samples greater than 15 of the number of OTUs in at least 2 samples

cd $currpath

## using the Singularity container

echo " - calling the singularity container"

singularity run ${sing_container} filter_otus_from_otu_table.py \
        -i ${output_dir}/4.OTU_picking/otu_table.biom \
        -n 15 \
        -s 2 \
        -o ${output_dir}/5.filter_OTUs/otu_table_filtered.biom

echo "DONE!"

echo " - converting biom file to tsv file"

cd $output_dir/5.filter_OTUs

singularity run ${sing_container} biom convert -i otu_table_filtered.biom -o otu_table_filtered.txt --to-tsv --header-key taxonomy

echo "DONE!"

echo "Excellent, your analysis with QIIME 1.9 are completed. REMEMBER TO RUN PROCESSING STATS LOCALLY. Have a nice day!"
