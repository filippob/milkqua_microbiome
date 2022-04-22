## Script that determines OTUS. OTUs were determined by aligning quality-filtered reads against the QIIME-compatible SILVA reference FASTA file, release 123, with minimum 97% clustering (https://www.arb-silva.de/download/archive/qiime/)

###NOTE: database can be changed

## setting up the environment
currpath=$(pwd)
project_home="/gpfs/home/users/chiara.gini/MILKQUA"
data_folder="Analysis/milkqua_skinswab/qiime1.9/3.quality_filtering"
output_dir="Analysis/milkqua_skinswab/qiime1.9"
param_file="temp/OnTest_Chiara/Databases/Silva_132/SILVA_132_QIIME_release"

cd $currpath
echo "project folder is $project_home"

## make folder if it does not exist

if [ ! -d "$project_home/${output_dir}/4.OTU_picking" ]; then
        mkdir -p $project_home/${output_dir}/4.OTU_picking
fi

## using the Singularity container

echo " - calling the singularity container"

singularity run /gpfs/software/Container/qiime_docker:fischuu-qiime-1.9.1.sif \
        pick_closed_reference_otus.py \
        --reference_fp $project_home/${param_file}/rep_set/rep_set_all/97/silva132_97.fna \
        --taxonomy_fp $project_home/${param_file}/taxonomy/taxonomy_all/97/raw_taxonomy.txt \
        --force --parallel --jobs_to_start=32 \
        --input_fp=$project_home/${data_folder}/seqs.fna \
        --output_dir=$project_home/${output_dir}/4.OTU_picking

echo "DONE!"
