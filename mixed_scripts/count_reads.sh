#!/bin/bash

## setting up the environment
currpath=$(pwd)
project_home="$HOME/MILKQUA"

target="${project_home}/Analysis/prova_qiime1.9/2.join_paired_ends"
output="$target/readsPerSample.tsv" 

## if file exists, remove it to avoid appending to old results
if [ -f $output ]; then
	rm $output
fi

touch $output
chmod a+rxw $output

## count reads per sample
for F in $target/*reads/*assembled.fastq.gz ; do
	fname="$(basename -- $F)"
	echo $fname 
	nrows=$(zcat -f < "$F" | wc -l)
	echo "N. lines: $nrows"
	nreads=$((nrows/4))
	echo "$F $nreads" >> $output
done;

## fix permissions to output file
#sed -i "s/$target//g" $output
sed -i 's,'"$target"',,' "$output"
sed -i 's/reads\/seqprep_assembled.fastq.gz//g' "$output"

echo "DONE!!"

