#!/bin/sh

## script that trims and truncates reads (quality filtering in the DADA2 framework: integrated quality filtering + read merging)
## with DADA2 the approach is OTU binning and then taxonomy classification (no clustering)

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/filter2"
output_dir="Analysis/milkqua_stools"
trim_forward=10
trim_reverse=10
trunc_forward=200
trunc_reverse=200
name="milkqua_stools2"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "${project_home}/${output_dir}/denoise2" ]; then
        mkdir -p ${project_home}/${output_dir}/denoise2
fi

echo " - running qiime for denoising"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime dada2 denoise-paired \
	--i-demultiplexed-seqs $project_home/$data_folder/${name}_filtered.qza \
	--p-trim-left-f $trim_forward \
	--p-trim-left-r $trim_reverse \
	--p-trunc-len-f $trunc_forward \
	--p-trunc-len-r $trunc_reverse \
	--o-table $project_home/$output_dir/denoise2/${name}_table.qza \
	--o-representative-sequences $project_home/$output_dir/denoise2/${name}_rep-seqs.qza \
	--o-denoising-stats $project_home/$output_dir/denoise2/${name}_denoising-stats.qza \
	--p-n-threads $SLURM_CPUS_PER_TASK

echo "DONE!"

