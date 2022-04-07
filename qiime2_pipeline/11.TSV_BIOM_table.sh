#!/bin/sh

## script that convert table.qza in OTU table

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools"
name="milkqua_stools2"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "${project_home}/${data_folder}/collapsed_biom_table" ]; then
        mkdir -p ${project_home}/${output_dir}/collapsed_biom_table
fi


echo " - running qiime for collapsing table and taxonomy"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        /gpfs/software/Container/qiime2_2022.2.sif qiime taxa collapse \
        --i-table $project_home/$data_folder/denoise2/${name}_table.qza \
        --i-taxonomy $project_home/$data_folder/taxonomy2/${name}_taxonomy.qza \
        --p-level 6 \
        --output-dir $project_home/$data_folder/collapsed_biom_table

echo "DONE!"

echo " - running qiime for export collapsed table into biom"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        /gpfs/software/Container/qiime2_2022.2.sif qiime tools export \
	--input-path $project_home/$data_folder/collapsed_biom_table/collapsed_table.qza \
	--output-path $project_home/$data_folder/collapsed_biom_table/converted_to_biom

echo "- running qiime to convert biom into tsv"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        /gpfs/software/Container/qiime2_2022.2.sif biom convert \
	-i $project_home/$data_folder/collapsed_biom_table/converted_to_biom/feature-table.biom \
	-o $project_home/$data_folder/collapsed_biom_table/feature-table.tsv \
	--to-tsv


echo "DONE!"

