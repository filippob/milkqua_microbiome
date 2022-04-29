#!/bin/sh

## script to donwload metadata from Google Drive spreadsheets
## python 3 is needed with installed dependencies
## !! IMPORTANT: before running, change parameters in the python script !!
proj_folder="$HOME/MILKQUA"

echo " - loading Python 3"
module load python3/intel/2020

echo " - moving to the project home folder"
cd $proj_folder

echo " - downloading the metadata"
python milkqua_microbiome/python_scripts/get_google_spreadsheet.py

echo "DONE!"

