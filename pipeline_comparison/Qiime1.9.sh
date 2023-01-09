

#!/bin/sh

## setting up the environment

## input/output
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder1="data/mock_communities"
#data_folder2=""
output_dir="Analysis/mock_communities/qiime1.9"
temp_folder="temp"
sample_start1=1 #first sample to use (in the sequence)
sample_end1=3 #last sample to use (in the sequence)
#sample_start2= #first sample to use (in the sequence)
#sample_end2= #last sample to use (in the sequence)

## software
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
repo_folder="${project_home}/milkqua_microbiome"
pipeline="${repo_folder}/qiime1.9_scripts"
rpath="${repo_folder}/mixed_scripts"

## parameters
#prefix="" # prefix to remove from sample file names (if any: !! usually this is left empty !!)
joinparam="Config/join.parameters"
filterparam="Config/parameters/filter.parameters"
dbpath="Databases/SILVA_132_QIIME"

## sequencing file name format
prefix=1 ## 1 if present (e.g. 'mock_1', 'mock_2' etc.); 0 if absent
interfix=0 ## ## 1 if present (e.g. '_S1_', '_S2_' etc.); 0 if absent
field=1 ## which field to retain after splitting by '_' (usually 1, for mock communities it is 2 --> mock_1, mock_2 etc., we want to retain 1, 2 etc. not mock, mock, mock ...
#r1="_L001_R1_001"
#r2="_L001_R2_001"
r1="R1"
r2="R2

cd $project_home
echo " - step 0: prepare files"
bash $pipeline/0.prepare_files.sh
echo "files have been prepared"

cd $project_home
echo " - step 1: extract barcodes"
bash $pipeline/1.extract_barcode.sh
echo "barcodes have been extracted"

cd $project_home
echo " - step 2: join paired-end reads"
bash $pipeline/2.join_paired_ends.sh
echo "reads have been joined"

cd $project_home
echo " - step 3: quality filtering"
bash $pipeline/3.quality_filtering.sh
echo "reads filtered and assembled into a single fasta file"

cd $project_home
echo " - step 4: closed-reference OTU picking"
bash $pipeline/4.OTU_picking.sh $dbpath
echo "OTU table is ready"

echo "DONE!"



