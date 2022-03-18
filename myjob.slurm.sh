#!/bin/bash

#SBATCH -o /gpfs/home/users/chiara.gini/MILKQUA/Analysis_subset/4.OTU_picking/test/myjob.%j.%N.out
#SBATCH -D /gpfs/home/users/chiara.gini
#SBATCH -J Prova
#SBATCH --get-user-env
#SBATCH -p light
#SBATCH --nodes=1
#SBATCH -c 1 # Number of cores per task
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-type=all
#SBATCH --mail-user=chiara.gini@unimi.it
#SBATCH --time=0:10:00
#SBATCH --account=MILKQUA

cd /gpfs/home/users/chiara.gini/test
/home/users/chiara.gini/MILKQUA/OnTest_Chiara/scripts_Qiime1.9/4.OTU_picking.sh
~                                                                                                    
~                                                                                                    
~                                                                                                    
~                                                                                                    
~     
