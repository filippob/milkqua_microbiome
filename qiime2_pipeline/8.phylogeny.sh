#!/bin/sh

## script that runs the phylogenetic analysis with qiime2

project_home="$HOME/MILKQUA"
data_folder="Analysis/milkqua_stools/denoise2"
output_dir="Analysis/milkqua_stools"
name="milkqua_stools2"
manifest="Config/manifest_faeces2.csv"

echo "project folder is $project_home"

## creating folder if not existing

if [ ! -d "$project_home/$output_dir/phylogeny2" ]; then
        mkdir -p $project_home/$output_dir/phylogeny2
fi

echo " - running qiime for the phylogenetic analysis"

singularity run \
	--env MPLCONFIGDIR=/gpfs/home/users/chiara.gini/MILKQUA/R_MYPACKAGES \
	/gpfs/software/Container/qiime2_2022.2.sif qiime phylogeny align-to-tree-mafft-fasttree \
	--i-sequences $project_home/$data_folder/${name}_rep-seqs.qza \
	--o-alignment $project_home/$output_dir/phylogeny2/${name}_aligned-rep-seqs.qza \
	--o-masked-alignment  $project_home/$output_dir/phylogeny2/${name}_masked-aligned-rep-seqs.qza \
	--o-tree  $project_home/$output_dir/phylogeny2/${name}_unrooted-tree.qza \
	--o-rooted-tree  $project_home/$output_dir/phylogeny2/${name}_rooted-tree.qza\
	--p-n-threads $SLURM_CPUS_PER_TASK

echo "DONE!"
