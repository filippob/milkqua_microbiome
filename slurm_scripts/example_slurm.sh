#!/bin/bash
#SBATCH --job-name=run-MILKQUA
#SBATCH --account=MILKQUA
#SBATCH --get-user-env
#SBATCH --partition=light
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:0:0
#SBATCH --output=/home/users/chiara.gini/MILKQUA/%x-%j.%N.out
#SBATCH --error=/home/users/chiara.gini/MILKQUA/%x-%j.%N.err
#SBATCH --mail-type=all
#SBATCH --mail-user=chiara.gini@unimi.it

sh /home/users/chiara.gini/MILKQUA/OnTest_Chiara/scripts_Qiime2/12.taxonomy.sh
                    
~                                                                                                    
~                                                                                                    
~                                                                                                    
~     
