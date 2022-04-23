#!/bin/bash

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"

target="${project_home}/Analysis/milkqua_skinswab/qiime1.9/3.quality_filtering"
output="$target/reads_after_filter.tsv" 

## if file exists, remove it to avoid appending to old results
if [ -f $output ]; then
	rm $output
fi

touch $output
chmod a+rxw $output

## count reads per sampleq
grep 'Sequence read filepath:' ${target}/split_library_log.txt >> $output
grep 'Total number of input sequences:' ${target}/split_library_log.txt >> $output
grep 'Total number seqs written' ${target}/split_library_log.txt >> $output

## fix permissions to output file
chmod g+xrw $output

echo "DONE!!"

