#!/bin/sh

## script that subsamples demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse) already imported into Qiime2
## sequences still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP_milk"
manifest="Config/manifest_milk.csv"
outdir="Analysis/milkqua_milk_subset"

if [ ! -d "${outdir}/subsample" ]; then
	mkdir -p ${outdir}/subsample
fi

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime for subsampling"
qiime demux subsample-paired --i-sequences $project_home/$outdir/import/$data_folder.qza --p-fraction 0.3 --verbose --o-subsampled-sequences $project_home/$outdir/subsample/${data_folder}_subsample.qza

echo " - running qiime for summary"
qiime demux summarize --i-data $project_home/$outdir/subsample/${data_folder}_subsample.qza --o-visualization $project_home/$outdir/subsample/${data_folder}_subsample.qzv

echo "DONE!"

