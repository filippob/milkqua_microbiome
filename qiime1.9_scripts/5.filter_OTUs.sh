###Script that performs filtering on by total count across samples greater than 15 of the number of OTUs in at least 2 samples

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_skinswab/qiime1.9/4.OTU_picking"
output_dir="Analysis/milkqua_skinswab/qiime1.9"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "$project_home/${output_dir}/5.filter_OTUs" ]; then
        mkdir -p $project_home/${output_dir}/5.filter_OTUs
        chmod g+rxw ${output_dir}
fi

## using the Singularity container

echo " - calling the singularity container"

singularity run ${sing_container} filter_otus_from_otu_table.py \
        -i $project_home/${data_folder}/otu_table.biom \
        -n 15 \
        -s 2 \
        -o $project_home/${output_dir}/5.filter_OTUs/otu_table_filtered.biom

echo "DONE!"
