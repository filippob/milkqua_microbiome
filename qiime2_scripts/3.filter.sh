#!/bin/sh

## script that filters out samples with few reads (less than 100 in this example)
## output saved in a qiime artifact that contains two files per sample: R1 and R2 (similar to the startign point from miSeq)

project_home="$HOME/MILKQUA"
data_folder="220225_M04028_0144_000000000-K6CMG"
min_seqs=100
manifest="Config/manifest_faeces.csv"
outdir="Analysis/milkqua_stools"
sif="/gpfs/software/Container/qiime2_2022.2.sif"
use_singularity=true

if [ ! -d "${outdir}/filter" ]; then
	mkdir -p ${outdir}/filter
fi

echo "project folder is $project_home"

if [ "${use_singularity}" = true ]; then
	echo "SINGULARITY (chosen qiime distribution)"
	echo " - running qiime for exporting data (to get information on number of reads per sample)"
	echo "(in the subfolder 'subsamples')"
	singularity run $sif \
	qiime tools export --input-path $project_home/$outdir/subsample/${data_folder}_subsample.qzv --output-path $project_home/$outdir/subsample/demux-subsample_temp/
	sed -i 's/,/\t/g' $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv
        tail -n +2 $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv > $project_home/$outdir/subsample/temp
        echo -e "sample-id\tforward sequence count\treverse sequence count" | cat - $project_home/$outdir/subsample/temp > $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv
        rm $project_home/$outdir/subsample/temp
	
	echo " - running qiime for filtering samples"
        echo " - removing samples with fewer than $min_seqs reads"
        echo "(in the subfolder 'filter')"
        singularity run $sif \
	qiime demux filter-samples --i-demux $project_home/$outdir/subsample/${data_folder}_subsample.qza \
        --m-metadata-file $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv \
        --p-where "CAST([forward sequence count] AS INT) > $min_seqs" --o-filtered-demux $project_home/$outdir/filter/${data_folder}_filtered.qza
else
	## exporting conda to PATH then activating the Qiime2 Conda env
	echo "CONDA ENV (chosen qiime distribution)"
	echo " - activating conda env"
	module load python3/intel/2020
	source activate qiime2-2019.10

	echo " - running qiime for exporting data (to get information on number of reads per sample)"
	echo "(in the subfolder 'subsamples')"
	qiime tools export --input-path $project_home/$outdir/subsample/${data_folder}_subsample.qzv --output-path $project_home/$outdir/subsample/demux-subsample_temp/
	sed -i 's/,/\t/g' $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv
	tail -n +2 $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv > $project_home/$outdir/subsample/temp
	echo -e "sample-id\tforward sequence count\treverse sequence count" | cat - $project_home/$outdir/subsample/temp > $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv
	rm $project_home/$outdir/subsample/temp

	echo " - running qiime for filtering samples"
	echo " - removing samples with fewer than $min_seqs reads"
	echo "(in the subfolder 'filter')"
	qiime demux filter-samples --i-demux $project_home/$outdir/subsample/${data_folder}_subsample.qza \
	--m-metadata-file $project_home/$outdir/subsample/demux-subsample_temp/per-sample-fastq-counts.tsv \
	--p-where "CAST([forward sequence count] AS INT) > $min_seqs" --o-filtered-demux $project_home/$outdir/filter/${data_folder}_filtered.qza

	echo " - removing temporary folders"
fi

echo " - removing temporary folders"
rm -r $project_home/$outdir/subsample/demux-subsample_temp

echo " - changing permissions to the output folder (r+x)"
chmod -R a+rx ${project_home}/${outdir}/filter

echo "DONE!"

