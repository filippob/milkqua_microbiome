#!/bin/sh

## script that summarises and visualises the results from DADA2 (feature table, representative sequences, quality filtering statistics)

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/denoise2"
output_dir="Analysis/milkqua_stools"
name="milkqua_stools2"
manifest="Config/manifest_faeces2.csv"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "${project_home}/${outdir}/feature_table2" ]; then
        mkdir -p ${project_home}/${outdir}/feature_table2
fi

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime feature-table summarize \
	--i-table $project_home/$data_folder/${name}_table.qza \
	--o-visualization $project_home/$output_dir/feature_table2/${name}_table.qzv \
	--m-sample-metadata-file $project_home/$manifest

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime feature-table tabulate-seqs \
	--i-data $project_home/$data_folder/${name}_rep-seqs.qza \
	--o-visualization $project_home/$output_dir/feature_table2/${name}_rep-seqs.qzv

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime metadata tabulate \
	--m-input-file $project_home/$data_folder/${name}_denoising-stats.qza \
	--o-visualization $project_home/$output_dir/feature_table2/${name}_denoising-stats.qzv

echo "DONE!"
