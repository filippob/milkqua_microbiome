#!/bin/sh

## script that imports into Qiime2 already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## sequence still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)
## type:
## input-format:

## Importing data with different formats

project_home="$HOME/MILKQUA"
data_folder="220225_M04028_0144_000000000-K6CMG"
manifest="Config/manifest_faeces.csv"
outdir="Analysis/milkqua_stools" ## <project-name>_<sample-type>

if [ ! -d "${outdir}/import" ]; then
	mkdir -p ${outdir}/import
fi

echo "project folder is $project_home"
echo "output folder is $outdir"

## exporting conda to PATH then activating the Qiime2 Conda env
echo " - activating conda env"
module load python3/intel/2020
source activate qiime2-2019.10

## importing data
echo " - running qiime to import data"

## CASAVA foramt
#qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path $HOME/MILKQUA/data/subset \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --output-path $HOME/MILKQUA/Analysis/demux-paired-end.qza

## Paired-end reads with quality (fastq files)
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $project_home/$manifest \
  --output-path $project_home/$outdir/import/${data_folder}.qza \
  --input-format PairedEndFastqManifestPhred33V2

echo " - changing permissions to the output folder (r+x)"
chmod -R a+rx ${project_home}/${outdir}/import

echo "DONE!"

