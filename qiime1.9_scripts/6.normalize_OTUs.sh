###Script that performs OTU counts normalizaiton

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_skinswab/qiime1.9/5.filter_OTUs"
output_dir="Analysis/milkqua_skinswab/qiime1.9"
sing_container="/gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif"
sing_container2="/home/users/chiara.gini/MILKQUA/singularity_containers/Qiime1.9.sif"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist
if [ ! -d "$project_home/${output_dir}/6.normalize_OTUs" ]; then
        mkdir -p $project_home/${output_dir}/6.normalize_OTUs
        chmod g+rxw ${output_dir}
fi

## using the Singularity container

echo " - calling the singularity container"

singularity run \
        --env  R_LIBS=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
        ${sing_container} normalize_table.py \
        -i $project_home/${data_folder}/otu_table_filtered.biom \
        -o $project_home/${output_dir}/6.normalize_OTUs/CSS_normalized_otu_table.biom

echo "DONE!"    
