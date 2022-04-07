#!/bin/sh

## script that prepares emperor plot

project_home="$HOME/MILKQUA"
data_folder="$project_home/Analysis_subset"
metadata="$project_home/Config/mapping_file.tsv"
database="$project_home/OnTest_Chiara/Databases"

echo "project folder is $project_home"
echo "- running qiime for emperor plot"

singularity run --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES  /gpfs/software/Container/qiime2_2022.2.sif qiime emperor plot \
  --i-pcoa $data_folder/core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file $metadata \
  --p-custom-axes treatment \
  --o-visualization $data_folder/core-metrics-results/unweighted-unifrac-emperor-treatment.qzv

singularity run --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES  /gpfs/software/Container/qiime2_2022.2.sif qiime emperor plot \
  --i-pcoa $data_folder/core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file $metadata \
  --p-custom-axes timepoint \
  --o-visualization $data_folder/core-metrics-results/bray-curtis-emperor-timepoint.qzv

echo "DONE!!"
