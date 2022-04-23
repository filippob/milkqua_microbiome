#!/bin/bash
## Script that determines OTUS. OTUs were determined by aligning quality-filtered reads against the QIIME-compatible SILVA reference FASTA file, release 123, with minimum 97% clustering (https://www.arb-silva.de/download/archive/qiime/)

###NOTE: database can be changed

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="Analysis/prova_qiime1.9/3.quality_filtering"
output_dir="Analysis/prova_qiime1.9/4.OTU_picking"
dbpath="Databases/SILVA_132/SILVA_132_QIIME_release"
#sing_container="${project_home}/Qiime1.9.sif"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist

if [ ! -d "$project_home/${output_dir}" ]; then
        mkdir -p $project_home/${output_dir}
        chmod g+rxw $project_home/${output_dir}
fi

## using the Singularity container

echo " - calling the singularity container for OTU picking"
singularity run ${sing_container} \
        pick_closed_reference_otus.py \
	--input_fp $project_home/${data_folder}/seqs.fna \
        --output_dir $project_home/${output_dir} \
	--reference_fp $project_home/$dbpath/rep_set/rep_set_all/97/silva132_97.fna \
        --taxonomy_fp $project_home/$dbpath/taxonomy/taxonomy_all/97/raw_taxonomy.txt \
        --force

chmod g+xrw -R $output_dir

echo " - converting biom file to tsv file"
cd $output_dir
singularity run ${sing_container} biom convert -i otu_table.biom -o otu_table.txt --to-tsv --header-key taxonomy

echo "DONE!"
