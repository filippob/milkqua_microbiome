#!/bin/bash

## script to download the reference sequence databases for metataxonomics
## e.g. SILVA, GreenGenes, Unite etc. (currently only SILVA db releases)

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"
target="${project_home}/Databases"
db_url="https://www.arb-silva.de/fileadmin/silva_databases/qiime/Silva_132_release.zip" 
folder_name="SILVA_132" ## name of folder where to unzip the downloaded database

cd $currpath

## if file exists, remove it to avoid appending to old results
if [ ! -d $target ]; then
	mkdir -p $target
	chmod g+rwx $target
fi

echo " - download database" 
cd $target
wget -O ${db_url}

## unzip the database archive
unzip $(basename -- ${db_url}) -d $folder_name
chmod g+rxw -R $folder_name

echo "DONE!!"

