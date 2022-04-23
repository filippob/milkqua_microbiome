#!/bin/bash
#SBATCH --job-name=run-join-MILKQUA
#SBATCH --account=MILKQUA
#SBATCH --get-user-env
#SBATCH --partition=light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:0:0
##SBATCH --output=/home/users/chiara.gini/MILKQUA/%x-%j.%N.out
##SBATCH --error=/home/users/chiara.gini/MILKQUA/%x-%j.%N.err
#SBATCH --output=/home/users/filippo.biscarini.est/MILKQUA/log/%x-%j.%N.out
#SBATCH --error=/home/users/filippo.biscarini.est/MILKQUA/log/%x-%j.%N.err
#SBATCH --mail-type=all
#SBATCH --mail-user=filippo.biscarini@gmail.com

bash $HOME/MILKQUA/milkqua_microbiome/qiime1.9_scripts/2.join_paired_ends.sh
