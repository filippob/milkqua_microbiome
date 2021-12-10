#!/bin/sh

## script that runs the phylogenetic analysis with qiime2

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime for the phylogenetic analysis"
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences $project_home/Analysis/${data_folder}_rep-seqs.qza --o-alignment $project_home/Analysis/${data_folder}_aligned-rep-seqs.qza --o-masked-alignment  $project_home/Analysis/${data_folder}_masked-aligned-rep-seqs.qza --o-tree  $project_home/Analysis/${data_folder}_unrooted-tree.qza --o-rooted-tree  $project_home/Analysis/${data_folder}_rooted-tree.qza

echo "DONE!"

