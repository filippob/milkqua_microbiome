#!/bin/sh

## script that imports into Qiime2 already demultiplexed data
## input data format is fastq: R1 and R2 files (forward and reverse)
## sequence still contain the barcode (to be trimmed later)
## (the function `qiime demux emp-paired` demultiplexes multiplexed data and leave the barcode intact)
## type:
## input-format:

## Importing data with different formats

project_home="$HOME/MILKQUA"
data_folder="210902_M04028_0139_000000000-JRGYP"

echo "project folder is $project_home"

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
  --input-path $project_home/Config/manifest.csv \
  --output-path $project_home/Analysis/$data_folder.qza \
  --input-format PairedEndFastqManifestPhred33V2

echo "DONE!"

