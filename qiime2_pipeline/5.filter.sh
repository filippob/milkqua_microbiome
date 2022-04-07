#!/bin/sh

## script that filters out samples with few reads (less than 100 in this example)
## output saved in a qiime artifact that contains two files per sample: R1 and R2 (similar to the startign point from miSeq)

currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/subsample2"
output_dir="Analysis/milkqua_stools/filter2"
min_seqs=100
manifest="Config/manifest_faeces2.csv"
name="milkqua_stools2"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "${project_home}/${output_dir}/filter2" ]; then
        mkdir -p ${project_home}/${output_dir}/filter2
fi

echo " - running qiime for exporting data (to get information on number of reads per sample)"

singularity run\
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime tools export \
	--input-path $project_home/$data_folder/${name}_subsample.qzv\
	--output-path $project_home/$output_dir/demux-subsample_temp/

sed -i 's/,/\t/g' $project_home/$output_dir/demux-subsample_temp/per-sample-fastq-counts.tsv

tail -n +2 $project_home/$output_dir/demux-subsample_temp/per-sample-fastq-counts.tsv > $project_home/$output_dir/temp

echo -e "sample-id\tforward sequence count\treverse sequence count" | cat - $project_home/$output_dir/temp > $project_home/$output_dir/demux-subsample_temp/per-sample-fastq-counts.tsv

echo " - running qiime for filtering samples"
echo " - removing samples with fewer than $min_seqs reads"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime demux filter-samples \
	--i-demux $project_home/$data_folder/${name}_subsample.qza \
	--m-metadata-file $project_home/$output_dir/demux-subsample_temp/per-sample-fastq-counts.tsv \
	--p-where "CAST([forward sequence count] AS INT) > $min_seqs" \
	--o-filtered-demux $project_home/$output_dir/${name}_filtered.qza

echo " - removing temporary folders"

rm -r $project_home/$output_dir/demux-subsample_temp
rm -r $project_home/$output_dir/temp

echo "DONE!"

