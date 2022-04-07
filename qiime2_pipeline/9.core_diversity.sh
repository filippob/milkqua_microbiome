#!/bin/bash

##Script that calculates core diversity metrics

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/denoise2"
input_dir="Analysis/milkqua_stools/phylogeny2"
output_dir="Analysis/milkqua_stools"
name="milkqua_stools2"
manifest="Config/manifest_faeces2.csv"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "${project_home}/${outdir}/core_metrics2" ]; then
        mkdir -p ${project_home}/${outdir}/core_metrics2
fi

echo " - running qiime to calculate core diversity metrics"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime diversity core-metrics-phylogenetic \
	--i-phylogeny $project_home/$input_dir/${name}_rooted-tree.qza \
	--i-table $project_home/$data_folder/${name}_table.qza \
	--p-sampling-depth 110 \
	--m-metadata-file $project_home/$manifest \
	--output-dir $project_home/$output_dir/core_metrics2/${name}_core-metrics-results


echo "DONE!!"
