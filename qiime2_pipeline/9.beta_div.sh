#!/bin/sh

## script that calculates beta diversity

project_home="$HOME/MILKQUA"
data_folder="$project_home/Analysis_subset"
metadata="$project_home/Config/mapping_file.tsv"
database="$project_home/OnTest_Chiara/Databases"

echo "project folder is $project_home"
echo "- running qiime for beta diversity"

singularity run --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES /gpfs/software/Container/qiime2_2022.2.sif qiime diversity beta-group-significance \
  --i-distance-matrix $data_folder/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $metadata \
  --o-visualization $data_folder/core-metrics-results/unweighted-unifrac-treatment-group-significance.qzv \
  --p-pairwise \
  --m-metadata-column treatment

singularity run --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES /gpfs/software/Container/qiime2_2022.2.sif qiime diversity beta-group-significance \
  --i-distance-matrix $data_folder/core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $metadata \
  --o-visualization $data_folder/core-metrics-results/unweighted-unifrac-timepoint-group-significance.qzv \
  --p-pairwise --m-metadata-column timepoint

echo "DONE!!"
