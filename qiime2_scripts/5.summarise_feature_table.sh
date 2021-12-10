#!/bin/sh

## script that trims and truncates reads (a little quality filtering)

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime to summarise the feature table"

qiime feature-table summarize --i-table $project_home/Analysis/${data_folder}_table.qza --o-visualization $project_home/Analysis/${data_folder}_table.qzv --m-sample-metadata-file $project_home/Config/mapping_file.tsv

qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv

echo "DONE!"

