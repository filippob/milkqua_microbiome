#!/bin/sh

## script that subsamples demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse) already imported into Qiime2
## sequences still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/import2"
output_dir="Analysis/milkqua_stools"
name="milkqua_stools2"

echo "project folder is $project_home"

cd $project_home/$data_folder

## creating folder if not existing

if [ ! -d "${project_home}/${output_dir}/subsample2" ]; then
	mkdir -p ${project_home}/${output_dir}/subsample2
fi

## running subsampling with 0.99 

echo " - running qiime for subsampling"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif \
	qiime demux subsample-paired \
	--i-sequences $project_home/$data_folder/${name}.qza \
	--p-fraction 0.99 \
	--verbose \
	--o-subsampled-sequences $project_home/$output_dir/subsample2/${name}_subsample.qza

echo " - running qiime for summary"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif \
	qiime demux summarize \
	--i-data $project_home/$output_dir/subsample2/${name}_subsample.qza \
	--o-visualization $project_home/$output_dir/subsample2/${name}_subsample.qzv

echo "DONE!"

