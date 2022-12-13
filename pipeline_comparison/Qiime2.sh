#!/bin/sh

## setting up the environment

currpath=$(pwd)
project_home="$HOME/MILKQUA"
input_folder="data/milk_subset"
manifest="Config/manifests/manifest_milk.tsv"
output_dir="Analysis/testing_q2_pipe/qiime2"
qiime2_cont="/gpfs/software/Container/qiime2_2022.2.sif"
name="milkqua_milk"
fraction=0.05
min_seqs=100
trim_forward=10
trim_reverse=10
trunc_forward=200
trunc_reverse=200
name="milkqua_milk"
database="Databases/Silva_138_QIIME2"
SLURM_CPUS_PER_TASK=24

cd $currpath
echo "we are in $currpath. Let's start by creating useful folders"

if [ ! -d "${output_dir}/1.import_fastq" ]; then
        mkdir -p ${output_dir}/1.import_fastq
        chmod g+rwx $output_dir/1.import_fastq
fi

if [ ! -d "${output_dir}/2.subsample" ]; then
        mkdir -p ${output_dir}/2.subsample
        chmod g+rwx $output_dir/2.subsample
fi

if [ ! -d "${output_dir}/3.filter" ]; then
        mkdir -p ${output_dir}/3.filter
        chmod g+rwx $output_dir/3.filter
fi

if [ ! -d "${output_dir}/4.denoise" ]; then
        mkdir -p ${output_dir}/4.denoise
        chmod g+rwx $output_dir/4.denoise
fi

if [ ! -d "${output_dir}/5.taxonomy" ]; then
        mkdir -p ${output_dir}/5.taxonomy
        chmod g+rwx $output_dir/5.taxonomy
fi

## importing data
echo " - running qiime to import data"

cd $currpath 

## Paired-end reads with quality (fastq files)

singularity run\
	 --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES\
	$qiime2_cont qiime tools import\
	--type 'SampleData[PairedEndSequencesWithQuality]'\
	--input-path $manifest \
	--output-path $output_dir/1.import_fastq/${name}.qza\
	--input-format PairedEndFastqManifestPhred33V2 

echo "DONE!"

## script that subsamples demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse) already imported into Qiime2
## sequences still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)

# running subsampling with $fraction

echo " - running qiime for subsampling"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont \
	qiime demux subsample-paired \
	--i-sequences $output_dir/1.import_fastq/${name}.qza \
	--p-fraction $fraction \
	--o-subsampled-sequences $output_dir/2.subsample/${name}_subsample$fraction.qza

echo " - running qiime for summary"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont \
	qiime demux summarize \
	--i-data $output_dir/2.subsample/${name}_subsample$fraction.qza \
	--o-visualization $output_dir/2.subsample/${name}_subsample$fraction.qzv	

echo "DONE!"

## script that filters out samples with few reads (less than 100 in this example)
## output saved in a qiime artifact that contains two files per sample: R1 and R2 (similar to the startign point from miSeq)

echo " - running qiime for exporting data (to get information on number of reads per sample)"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime tools export \
	--input-path $output_dir/2.subsample/${name}_subsample$fraction.qzv \
	--output-path $output_dir/3.filter/demux-subsample_temp$fraction/

sed -i 's/,/\t/g' $output_dir/3.filter/demux-subsample_temp$fraction/per-sample-fastq-counts.tsv

tail -n +2 $output_dir/3.filter/demux-subsample_temp$fraction/per-sample-fastq-counts.tsv > $output_dir/3.filter/temp

echo -e "sample-id\tforward sequence count\treverse sequence count" | cat - $output_dir/3.filter/temp > $output_dir/3.filter/demux-subsample_temp$fraction/per-sample-fastq-counts.tsv

echo " - running qiime for filtering samples. Removing samples with fewer than $min_seqs reads"

singularity run \
 	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime demux filter-samples \
	--i-demux $output_dir/2.subsample/${name}_subsample$fraction.qza \
	--m-metadata-file $output_dir/3.filter/demux-subsample_temp$fraction/per-sample-fastq-counts.tsv \
	--p-where "CAST([forward sequence count] AS INT) > $min_seqs" \
	--o-filtered-demux $output_dir/3.filter/${name}_filtered.qza

echo " - removing temporary folders"

rm -r $output_dir/3.filter/demux-subsample_temp$fraction
rm -r $output_dir/3.filter/temp

echo "DONE!"

## script that trims and truncates reads (quality filtering in the DADA2 framework: integrated quality filtering + read merging)
## with DADA2 the approach is OTU binning and then taxonomy classification (no clustering)

echo " - running qiime for denoising"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime dada2 denoise-paired \
	--i-demultiplexed-seqs $output_dir/3.filter/${name}_filtered.qza \
	--p-trim-left-f $trim_forward \
	--p-trim-left-r $trim_reverse \
	--p-trunc-len-f $trunc_forward \
	--p-trunc-len-r $trunc_reverse \
	--o-table $output_dir/4.denoise/${name}_table.qza \
	--o-representative-sequences $output_dir/4.denoise/${name}_rep-seqs.qza \
	--o-denoising-stats $output_dir/4.denoise/${name}_denoising-stats.qza \
	--p-n-threads $SLURM_CPUS_PER_TASK

echo "DONE!"

## script that summarises and visualises the results from DADA2 (feature table, representative sequences, quality filtering statistics)

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime feature-table summarize \
	--i-table $output_dir/4.denoise/${name}_table.qza \
	--o-visualization $output_dir/4.denoise/${name}_table.qzv \
	--m-sample-metadata-file $manifest

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime feature-table tabulate-seqs \
	--i-data $output_dir/4.denoise/${name}_rep-seqs.qza \
	--o-visualization $output_dir/4.denoise/${name}_rep-seqs.qzv

echo " - running qiime to summarise the feature table"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime metadata tabulate \
	--m-input-file $output_dir/4.denoise/${name}_denoising-stats.qza \
	--o-visualization $output_dir/4.denoise/${name}_denoising-stats.qzv

echo "DONE!"

## script that obtain taxonomy analysis

echo "- running qiime for taxonomic analysis"

export TMPDIR=/home/users/chiara.gini/MILKQUA/scratch

cd $currpath

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime feature-classifier classify-sklearn \
  	--i-classifier $database/silva-138-99-515-806-nb-classifier.qza \
  	--i-reads $output_dir/4.denoise/${name}_rep-seqs.qza  \
  	--o-classification $output_dir/5.taxonomy/${name}_taxonomy.qza

echo "- running qiime to summarize taxonomic analysis"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	$qiime2_cont qiime metadata tabulate \
	--m-input-file $output_dir/5.taxonomy/${name}_taxonomy.qza \
  	--o-visualization $output_dir/5.taxonomy/${name}_taxonomy.qzv

echo "DONE!!"

## script that convert table.qza in OTU table

echo " - running qiime for collapsing table and taxonomy"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        $qiime2_cont qiime taxa collapse \
        --i-table $output_dir/4.denoise/${name}_table.qza \
        --i-taxonomy $output_dir/5.taxonomy/${name}_taxonomy.qza \
        --p-level 6 \
        --output-dir $output_dir/5.taxonomy/collapsed_biom_table

echo "DONE!"

echo " - running qiime for export collapsed table into biom"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        $qiime2_cont qiime tools export \
	--input-path $output_dir/5.taxonomy/collapsed_biom_table/collapsed_table.qza \
	--output-path $output_dir/5.taxonomy/collapsed_biom_table/converted_to_biom

echo "- running qiime to convert biom into tsv"

singularity run \
        --env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
      	$qiime2_cont biom convert \
	-i $output_dir/5.taxonomy/collapsed_biom_table/converted_to_biom/feature-table.biom \
	-o $output_dir/5.taxonomy/collapsed_biom_table/feature-table$fraction.tsv \
	--to-tsv


echo "DONE!"

echo "Excellent, your analysis with QIIME 2 are completed. Have a nice day!"
