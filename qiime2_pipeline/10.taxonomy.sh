#!/bin/sh

## script that obtain taxonomy analysis

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/denoise2"
output_dir="Analysis/milkqua_stools"
name="milkqua_stools2"
manifest="Config/manifest_faeces2.csv"
database="OnTest_Chiara/Databases/Silva_138_QIIME2"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "$project_home/$output_dir/taxonomy2" ]; then
        mkdir -p $project_home/$output_dir/taxonomy2
fi

echo "project folder is $project_home"
echo "- running qiime for taxonomic analysis"

export TMPDIR=/home/users/chiara.gini/MILKQUA/scratch

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime feature-classifier classify-sklearn \
  	--i-classifier $project_home/$database/silva-138-99-515-806-nb-classifier.qza \
  	--i-reads $project_home/$data_folder/${name}_rep-seqs.qza \
  	--o-classification $project_home/$output_dir/taxonomy2/${name}_taxonomy.qza\

echo "- running qiime to summarize taxonomic analysis"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif  qiime metadata tabulate \
	--m-input-file $project_home/$output_dir/taxonomy2/${name}_taxonomy.qza \
  	--o-visualization $project_home/$output_dir/taxonomy2/${name}_taxonomy.qzv\

echo "DONE!!_"
