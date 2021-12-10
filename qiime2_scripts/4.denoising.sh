#!/bin/sh

## script that trims and truncates reads (quality filtering in the DADA2 framework: integrated quality filtering + read merging)
## with DADA2 the approach is OTU binning and then taxonomy classification (no clustering)

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"
trim_forward=10
trim_reverse=10
trunc_forward=200
trunc_reverse=200

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime for denoising"
qiime dada2 denoise-paired --i-demultiplexed-seqs $project_home/Analysis/${data_folder}_filtered.qza --p-trim-left-f $trim_forward --p-trim-left-r $trim_reverse --p-trunc-len-f $trunc_forward --p-trunc-len-r $trunc_reverse --o-table $project_home/Analysis/${data_folder}_table.qza --o-representative-sequences $project_home/Analysis/${data_folder}_rep-seqs.qza --o-denoising-stats $project_home/Analysis/${data_folder}_denoising-stats.qza

echo "DONE!"

