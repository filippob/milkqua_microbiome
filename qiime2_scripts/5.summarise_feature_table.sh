#!/bin/sh

## script that summarises and visualises the results from DADA2 (feature table, representative sequences, quality filtering statistics)

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"
manifest="Config/manifest_milk.csv"
outdir="Analysis/milkqua_milk_subset"

if [ ! -d "${outdir}/filter" ]; then
	mkdir -p ${outdir}/filter
fi

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime to summarise the feature table"
echo "(input from denoise/)"
qiime feature-table summarize --i-table $project_home/$outdir/denoise/${data_folder}_table.qza --o-visualization $project_home/$outdir/denoise/${data_folder}_table.qzv --m-sample-metadata-file $project_home/Config/mapping_file.tsv
qiime feature-table tabulate-seqs --i-data $project_home/$outdir/denoise/${data_folder}_rep-seqs.qza --o-visualization $project_home/$outdir/denoise/${data_folder}_rep-seqs.qzv
qiime metadata tabulate --m-input-file $project_home/$outdir/denoise/${data_folder}_denoising-stats.qza --o-visualization $project_home/$outdir/denoise/${data_folder}_denoising-stats.qzv

echo "DONE!"

