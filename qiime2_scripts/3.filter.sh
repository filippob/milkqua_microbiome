#!/bin/sh

## script that filters out samples with few reads (less than 100 in this example)
## output saved in a qiime artifact that contains two files per sample: R1 and R2 (similar to the startign point from miSeq)

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP_milk"
min_seqs=100
manifest="Config/manifest_milk.csv"
outdir="Analysis/milkqua_milk_subset"

if [ ! -d "${outdir}/filter" ]; then
	mkdir -p ${outdir}/filter
fi

echo "project folder is $project_home"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

echo " - running qiime for exporting data (to get information on number of reads per sample)"
echo "(in the subfolder 'subsamples')"
qiime tools export --input-path $project_home/$outdir/subsample/${data_folder}_subsample.qzv --output-path $project_home/$outdir/subsample/demux-subsample_temp/
sed -i 's/,/\t/g' $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.csv
tail -n +2 $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.csv > $project_home/$outdir/subsample/temp
echo -e "sample-id\tforward sequence count" | cat - $project_home/$outdir/subsample/temp > $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.csv
rm $project_home/$outdir/subsample/temp

echo " - running qiime for filtering samples"
echo " - removing samples with fewer than $min_seqs reads"
echo "(in the subfolder 'filter')"
qiime demux filter-samples --i-demux $project_home/$outdir/subsample/${data_folder}_subsample.qza --m-metadata-file $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.csv --p-where "CAST([forward sequence count] AS INT) > $min_seqs" --o-filtered-demux $project_home/$outdir/filter/${data_folder}_filtered.qza

echo " - removing temporary folders"
#rm -r $project_home/$outdir/subsample/demux-subsample_temp

echo "DONE!"

