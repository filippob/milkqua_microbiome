#!/bin/sh

# setting up the environment

currpath=$(pwd)
project_home="$HOME/MILKQUA"
data_folder1="data/subset" ## in some experiments you may have multiple data folders (hence the numbering: 1, 2 etc.)
output_dir="Analysis/milkqua_milk/micca"
temp_folder="temp"
sample_start1=18 #first sample to use (in the sequence)
sample_end1=72 #last sample to use (in the sequence)
fwd_primer="CCTACGGGNGGCWGCAG" #(adapter: TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG)
rev_primer="GACTACHVGGGTATCTAATCC" #(adapter: GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG)
sing_container="singularity_containers/Micca.sif"
sickle_exe="singularity_containers/sickle"
core=8
dbpath="Databases/SILVA_132_QIIME_release"
Q=20 ## Phred score threshold

cd $project_home

echo "##################################"
echo "BEGINNING THE ANALYSIS WITH MICCA"
echo "##################################"
echo "we are in $(pwd). Let's start by creating useful folders"

if [ ! -d "$output_dir/$temp_folder" ]; then
        mkdir -p $output_dir/$temp_folder
	chmod g+rwx $output_dir/$temp_folder 
fi

if [ ! -d "$output_dir/1.cutadapt" ]; then
        mkdir -p $output_dir/1.cutadapt
	chmod g+rwx $output_dir/1.cutadapt
fi

if [ ! -d "$output_dir/2.trimming" ]; then
        mkdir -p $output_dir/2.trimming
	chmod g+rwx $output_dir/2.trimming 
fi

if [ ! -d "$output_dir/3.mergepairs" ]; then
        mkdir -p $output_dir/3.mergepairs
	chmod g+rwx $output_dir/3.mergepairs 
fi

if [ ! -d "$output_dir/4.filter" ]; then
        mkdir -p $output_dir/4.filter
	chmod g+rwx $output_dir/4.filter 
fi

if [ ! -d "$output_dir/5.otu_picking" ]; then
        mkdir -p $output_dir/5.otu_picking
	chmod g+rwx $output_dir/5.otu_picking 
fi

if [ ! -d "$output_dir/6.classify" ]; then
        mkdir -p $output_dir/6.classify
	chmod g+rwx $output_dir/6.classify 
fi

if [ ! -d "$output_dir/Quality_control" ]; then
        mkdir -p $output_dir/Quality_control
	chmod g+rwx $output_dir/Quality_control 
fi

cd $project_home

####################################
### STEP 0: COPYING AND RENAMING ###
####################################
## Script to rename samples and count input reads

echo "we are in $(pwd). Let's rename the samples"

este=".fastq.gz"

echo " - copying files from seq data folder to temp folder"
for sample in ${data_folder1}/*.fastq.gz
do
  echo "copying to a temporary folder the file $sample"
  cp $sample $output_dir/$temp_folder
  echo -e "$sample\t-->\t$sample"_"$read$este" >> $output_dir/$temp_folder/log_renamer.txt
done

## Count reads
echo " - counting initial reads from sequencing"
cd $output_dir/$temp_folder
for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_raw.txt
        echo -n " " >> seq_count_16S_raw.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_raw.txt
done

echo "Copying and renaming: DONE"

########################################
### STEP 1: REMOVING PRIMERS & CO.   ###
########################################
cd $project_home
echo "we are in $(pwd). Let's create single_names file and then remove adapters"
## Remove adapters. Create a file of names that will be used for looping. Only file/sample name, remove extension and R1/R2
echo " - removing primer and adapters (and barcodes)"
cd $output_dir/$temp_folder

for i in *.fastq.gz
do
	echo "$i" | cut -d "_" -f1 >> names.txt
	sed 'n; d' names.txt > names_single.txt
done

# remove primers with cutadapt
# Primers (Sequences from Pindo and FMACH): forward: CCTACGGGNGGCWGCAG, reverse: GACTACNVGGGTWTCTAATCC

cd $project_home

echo "Running cutadapt by calling singularity container"

while read -r line;
do
	echo "removing adapters from file $line";

	singularity run singularity_containers/Micca.sif cutadapt -g Forward=$fwd_primer -G Reverse=$rev_primer --discard-untrimmed --pair-filter=any -o ${output_dir}/1.cutadapt/${line}_R1_cutadapt.fastq.gz -p ${output_dir}/1.cutadapt/${line}_R2_cutadapt.fastq.gz ${output_dir}/temp/${line}_S${line}_L001_R1_001.fastq.gz ${output_dir}/temp/${line}_S${line}_L001_R2_001.fastq.gz >> ${output_dir}/Quality_control/cutadapt_report.txt

done <  ${output_dir}/$temp_folder/names_single.txt

# --discard-untrimmed, --trimmed-only
#                        Discard reads that do not contain an adapter.

#--pair-filter=(any|both|first)
#                        Which of the reads in a paired-end read have to match
#                        the filtering criterion in order for the pair to be
#                        filtered. Default: any

echo "Cutadapt: DONE"

##########################################
### STEP 2: TRIMMING READS FOR QUALITY ###
##########################################
## TRIMMING. # trim low quality part. Q = 20 inspect quality, eventually for 16S Q can be set to 25
echo " - trimming for quality parameters"
echo "we are in $(pwd). Let's do the trimming for quality. Q is set at $Q" 
cd $project_home

while read -r line;
do
	#echo "Running sickle on file "$line""
	echo "Running sickle on file "$line"" >> ${output_dir}/Quality_control/stats_trim.txt
	echo "Running sickle by calling singularity container on file ${line}"

	singularity run $sickle_exe sickle pe \
	-f ${output_dir}/1.cutadapt/${line}_R1_cutadapt.fastq.gz\
	-r ${output_dir}/1.cutadapt/${line}_R2_cutadapt.fastq.gz \
	-o ${output_dir}/2.trimming/${line}_trimmed_R1.fastq.gz \
	-p ${output_dir}/2.trimming/${line}_trimmed_R2.fastq.gz \
	-s ${output_dir}/2.trimming/${line}_singles.gz \
	-t sanger -q $Q -g 1>> ${output_dir}/Quality_control/stats_trim.txt

done < ${output_dir}/$temp_folder/names_single.txt

echo "TRIMMING: DONE"

cd $project_home

# Count trimmed data
echo " - counting data after trimming"
cd ${output_dir}/2.trimming

echo " - counting sequences "
for i in *.fastq.gz
do
        echo -n $i >> seq_count_16S_QC.txt
        echo -n " " >> seq_count_16S_QC.txt
        echo $(zcat $i | wc -l) / 4 | bc >> seq_count_16S_QC.txt
done
echo "Counting after trimming: Done"
cd $project_home

####################################
### STEP 3: JOINING READS        ###
####################################
## Join reads (MICCA)
echo " - joining R1 and R2 reads files"
echo "we are in $(pwd). Let's join reads!"
 
cd ${output_dir}/2.trimming

# remove singles reads from sickle
echo " -  removing single reads from sickle"
rm *.singles.gz

echo " - uncompressing trimmed fastq files "
cd $project_home
gunzip $output_dir/2.trimming/*.fastq.gz

echo " - joining reads"
cd $project_home

singularity run $sing_container micca mergepairs \
	-i ${output_dir}/2.trimming/*_R1.fastq \
	-o ${output_dir}/3.mergepairs/assembled_16S.fastq \
	-l 32 -d 8 -t 7

# -l : minimum overlap between reads
# -d : maximum mismatch in overlap region
# -t: n. of threads

# Counting reads in assembled file
echo " - counting reads after joining "
grep -c '^@M' $output_dir/3.mergepairs/assembled_16S.fastq

# zipping back trimmed files

echo " - recompressing trimmed fastq files "
cd $project_home
gzip ${output_dir}/2.trimming/*.fastq

echo "JOINING: DONE"

####################################
### STEP 4: DATA FILTERING       ###
####################################
## Filter 
cd $project_home
echo " - filtering assembled reads data"
echo "we are in $(pwd). Let's do the filtering"

# Remove N from assembly (uncalled bases)
singularity run $sing_container micca filter \
	-i $output_dir/3.mergepairs/assembled_16S.fastq \
	-o $output_dir/4.filter/assembled_16S_filtered.fasta \
	--maxns 0

# count
echo " - counting reads after filtering"
grep -c '>' $output_dir/4.filter/assembled_16S_filtered.fasta

echo "FILTERING: DONE"

#####################################
### STEP 5: OTU PICKING (BINNING) ###
#####################################
## OTU picking
echo " - OTU picking"
echo "we are in $currpath. Let's do the OTU picking"
 
# pick otu
singularity run $sing_container micca otu -m denovo_unoise \
	-i $output_dir/4.filter/assembled_16S_filtered.fasta \
	-o ${output_dir}/5.otu_picking \
	-t 8 --rmchim

echo "OTU-PICKING: DONE"

#####################################
### STEP 6: CLASSIFY OTU          ###
#####################################
## classify RDP
cd $project_home

echo " - classifying OTUs"
echo "we are in $(pwd). Let's align to the database. Today it will be $dbpath"

## CLASSIFY WITH VSEARCH AND SILVA!
# QIIME compatible SILVa DB should be downloaded

singularity run $sing_container micca classify \
	-m cons \
	-i ${output_dir}/5.otu_picking/otus.fasta  \
	-o $output_dir/6.classify/taxa_SILVA.txt \
	--ref $dbpath/rep_set/rep_set_16S_only/97/silva_132_97_16S.fna \
	--ref-tax $dbpath/taxonomy/16S_only/97/taxonomy_7_levels.txt

#ALTERNATIVE
#singularity run $currpath/micca.sif micca classify -m rdp -i $currpath/Analysis/micca_its/otus.fasta --rdp-gene 16srrna -o $currpath/Analysis/micca_its/taxa.txt

echo "Classify: DONE"

echo "Excellent, your analysis with MICCA are completed. Have a nice day!"
