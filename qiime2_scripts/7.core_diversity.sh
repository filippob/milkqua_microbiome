#!/bin/sh

## script that calculates ore diversity metrics

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime to calculate core diversity metrics"
qiime diversity core-metrics-phylogenetic --i-phylogeny $project_home/Analysis/${data_folder}_rooted-tree.qza --i-table $project_home/Analysis/${data_folder}_table.qza --p-sampling-depth 110 --m-metadata-file $project_home/Config/mapping_file.tsv --output-dir $project_home/Analysis/${data_folder}_core-metrics-results

echo "DONE!"

