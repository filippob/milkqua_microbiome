#!/bin/sh                                                                       

# This is an example of a pipeline using vsearch to process data in the Mothur 16S rRNA MiSeq SOP tutorial dataset to perform initial paired-end read merging, quality filtering, chimera removal and OTU clustering.          

THREADS=1
REF="gold.fasta"
PERL=$(which perl)
VSEARCH="singularity_containers/Vsearch.sif"
MISEQSOPDATA="https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
GOLD="https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip"
currpath=$(pwd)
output_dir="temp"
config="Config"
date

cd $output_dir

echo Obtaining Mothur MiSeq SOP tutorial dataset

if [ ! -e miseqsopdata.zip ]; then
    wget $MISEQSOPDATA
fi

echo Decompressing...

unzip -u -o miseqsopdata.zip
rm -rf miseqsopdata.zip

#mv MiSeq_SOP $output_dir

echo Obtaining Gold reference database for chimera detection

if [ ! -e gold.fasta ]; then
    if [ ! -e silva.gold.bacteria.zip ]; then
        wget $GOLD
    fi

    echo Decompressing and reformatting...
    unzip -p silva.gold.bacteria.zip silva.gold.align | \
        sed -e "s/[.-]//g" > gold.fasta

fi

#mv gold.fasta $output_dir
rm -rf silva.gold.bacteria.zip

# Enter subdirectory                                                            

#cd MiSeq_SOP

cd $currpath

echo Checking FASTQ format version for one file

singularity run $VSEARCH vsearch --fastq_chars $(ls -1 $output_dir/MiSeq_SOP/*.fastq | head -1)

# Process samples                                                               

for f in $output_dir/MiSeq_SOP/*_R1_*.fastq; do

    r=$(sed -e "s/_R1_/_R2_/" <<< "$f")
    s=$(cut -d_ -f1 <<< "$f")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================

   singularity run $VSEARCH vsearch --fastq_mergepairs $f \
        --threads $THREADS \
        --reverse $r \
        --fastq_minovlen 200 \
        --fastq_maxdiffs 15 \
        --fastqout $s.merged.fastq \
        --fastq_eeout

    # Commands to demultiplex and remove tags and primers                       
    # using e.g. cutadapt may be added here.                                    

    echo
    echo Calculate quality statistics

   singularity run $VSEARCH vsearch --fastq_eestats $s.merged.fastq \
        --output $s.stats
    echo
    echo Quality filtering

    singularity run $VSEARCH vsearch --fastq_filter $s.merged.fastq \
        --fastq_maxee 0.5 \
        --fastq_minlen 225 \
        --fastq_maxlen 275 \
        --fastq_maxns 0 \
        --fastaout $s.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

   singularity run $VSEARCH vsearch --derep_fulllength $s.filtered.fasta \
        --strand plus \
        --output $s.derep.fasta \
        --sizeout \
        --uc $s.derep.uc \
        --relabel $s. \
        --fasta_width 0

done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Merge all samples

rm -f $output_dir/all.derep.fasta $output_dir/all.nonchimeras.derep.fasta
cat $output_dir/*.derep.fasta > $output_dir/all.fasta

echo
echo Dereplicate across samples and remove singletons

singularity run $VSEARCH vsearch --derep_fulllength $output_dir/all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $output_dir/all.derep.uc \
    --output $output_dir/all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" $output_dir/all.derep.fasta)

echo
echo Precluster at 98% before chimera detection

singularity run $VSEARCH vsearch --cluster_size $output_dir/all.derep.fasta \
    --threads $THREADS \
    --id 0.98 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $output_dir/all.preclustered.uc \
    --centroids $output_dir/all.preclustered.fasta

echo Unique sequences after preclustering: $(grep -c "^>" $output_dir/all.preclustered.fasta)

echo
echo De novo chimera detection

singularity run $VSEARCH vsearch --uchime_denovo $output_dir/all.preclustered.fasta \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras $output_dir/all.denovo.nonchimeras.fasta \

echo Unique sequences after de novo chimera detection: $(grep -c "^>" $output_dir/all.denovo.nonchimeras.fasta)

echo
echo Reference chimera detection

singularity run $VSEARCH vsearch --uchime_ref $output_dir/all.denovo.nonchimeras.fasta \
    --threads $THREADS \
    --db $output_dir/$REF \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --nonchimeras $output_dir/all.ref.nonchimeras.fasta

echo Unique sequences after reference-based chimera detection: $(grep -c "^>" $output_dir/all.ref.nonchimeras.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences, dereplicated

$PERL $config/map.pl $output_dir/all.derep.fasta $output_dir/all.preclustered.uc $output_dir/all.ref.nonchimeras.fasta > $output_dir/all.nonchimeras.derep.fasta

echo Unique non-chimeric, non-singleton sequences: $(grep -c "^>" $output_dir/all.nonchimeras.derep.fasta)

echo
echo Extract all non-chimeric, non-singleton sequences in each sample

$PERL $config/map.pl $output_dir/all.fasta $output_dir/all.derep.uc $output_dir/all.nonchimeras.derep.fasta > $output_dir/all.nonchimeras.fasta

echo Sum of unique non-chimeric, non-singleton sequences in each sample: $(grep -c "^>" all.nonchimeras.fasta)

echo
echo Cluster at 97% and relabel with OTU_n, generate OTU table

singularity run $VSEARCH vsearch --cluster_size $output_dir/all.nonchimeras.fasta \
    --threads $THREADS \
    --id 0.97 \
    --strand plus \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc $output_dir/all.clustered.uc \
    --relabel OTU_ \
    --centroids $output_dir/all.otus.fasta \
    --otutabout $output_dir/all.otutab.txt

echo
echo Number of OTUs: $(grep -c "^>" $output_dir/all.otus.fasta)

echo
echo Done

date
