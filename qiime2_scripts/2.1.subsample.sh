#!/bin/sh

## script that subsamples demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse) already imported into Qiime2
## sequences still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)

project_home="$HOME/MILKQUA"
data_folder="220225_M04028_0144_000000000-K6CMG"
manifest="Config/manifest_faeces.csv"
outdir="Analysis/milkqua_stools"
fraction=0.33
use_singularity=true #if false the conda env is used

if [ ! -d "${outdir}/subsample" ]; then
	mkdir -p ${outdir}/subsample
fi

echo "project folder is $project_home"

if [ "${use_singularity}" = true ]; then
	
	echo "SINGULARITY (chosen qiime distribution)"
	echo " - running qiime for subsampling"	
	singularity run /gpfs/software/Container/qiime2_2022.2.sif \
	qiime demux subsample-paired --i-sequences $project_home/$outdir/import/${data_folder}.qza \
	 --p-fraction $fraction --verbose --o-subsampled-sequences $project_home/$outdir/subsample/${data_folder}_subsample.qza

	echo " - running qiime for summary"
	singularity run /gpfs/software/Container/qiime2_2022.2.sif \
	qiime demux summarize --i-data $project_home/$outdir/subsample/${data_folder}_subsample.qza --o-visualization $project_home/$outdir/subsample/${data_folder}_subsample.qzv

else
	echo "CONDA ENV (chosen qiime distribution)"
	## exporting conda to PATH then activating the Qiime2 Conda env
	echo " - activating conda env"
	module load python3/intel/2020
	source activate qiime2-2019.10

	echo " - running qiime for subsampling"
	qiime demux subsample-paired --i-sequences $project_home/$outdir/import/${data_folder}.qza \
	--p-fraction $fraction --verbose --o-subsampled-sequences $project_home/$outdir/subsample/${data_folder}_subsample.qza

	echo " - running qiime for summary"
	qiime demux summarize --i-data $project_home/$outdir/subsample/${data_folder}_subsample.qza --o-visualization $project_home/$outdir/subsample/${data_folder}_subsample.qzv

fi

echo " - changing permissions to the output folder (r+x)"
chmod -R a+rx ${project_home}/${outdir}/subsample

echo "DONE!"

