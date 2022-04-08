#!/bin/sh

## script that trims and truncates reads (quality filtering in the DADA2 framework: integrated quality filtering + read merging)
## with DADA2 the approach is OTU binning and then taxonomy classification (no clustering)

project_home="$HOME/MILKQUA"
data_folder="220225_M04028_0144_000000000-K6CMG"
manifest="Config/manifest_faeces.csv"
outdir="Analysis/milkqua_stools"
trim_forward=10
trim_reverse=10
trunc_forward=200
trunc_reverse=200
sif="/gpfs/software/Container/qiime2_2022.2.sif"
use_singularity=true

echo "project folder is $project_home"

if [ ! -d "${outdir}/denoise" ]; then
	mkdir -p ${outdir}/denoise
fi

if [ "${use_singularity}" = true ]; then

else
	echo "CONDA ENV (chosen qiime distribution)"
	## exporting conda to PATH then activating the Qiime2 Conda env
	echo " - activating conda env"
	module load python3/intel/2020
	source activate qiime2-2019.10

	echo " - running qiime for denoising"
	qiime dada2 denoise-paired --i-demultiplexed-seqs $project_home/Analysis/${data_folder}_filtered.qza \
	--p-trim-left-f $trim_forward --p-trim-left-r $trim_reverse --p-trunc-len-f $trunc_forward --p-trunc-len-r $trunc_reverse \
	--o-table $project_home/Analysis/${data_folder}_table.qza --o-representative-sequences $project_home/Analysis/${data_folder}_rep-seqs.qza \
	--o-denoising-stats $project_home/Analysis/${data_folder}_denoising-stats.qza

fi

echo " - changing permissions to the output folder (r+x)"
chmod -R a+rx ${project_home}/${outdir}/subsample

echo "DONE!"

