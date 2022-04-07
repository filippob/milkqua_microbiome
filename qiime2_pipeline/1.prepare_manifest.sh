#!/bin/bash

project_home="$HOME/MILKQUA"
r_pckgs="${project_home}/r_packages"
rscript="milkqua_microbiome/r_scripts/create_manifest.r"

## setting up R
module load R/3.6.0
export R_LIBS=$r_pckgs

echo " - running Rscript"
Rscript --vanilla ${project_home}/$rscript

echo "DONE!!"
