#!/bin/sh  

### This is Mothur pipeline. Guidelines were retrieved from mothur official website, under MiSeq SOP. Althought 1st version of mothur was released in 2013, developers from michigan institute kept it updeated. Latest release is v.1.48.0 released in May 2022.

##Cite: Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013): Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. Applied and Environmental Microbiology. 79(17):5112-20.
        #REPORT THE DATE YOU ACCESSED Miseq SOP 

mothur="singularity_containers/Mothurv1.48.sif"
dataset="https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip"
db="Databases"
silva="https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.bacteria.zip"
rdp="https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip"
currpath=$(pwd)
output_dir="Analysis/mock_communities/mothur"
config="Config"

date
cd $currpath
echo "$currpath"
echo "Retrieve dataset and databases"

if [ ! -e $output_dir ]; then
	mkdir $output_dir
	chmod g+rwx $output_dir
fi

echo "Obtaining Mothur MiSeq SOP tutorial dataset"

cd $output_dir

if [ ! -e miseqsopdata.zip ]; then
    wget $dataset
fi

echo "Decompressing..."

unzip -u -o miseqsopdata.zip
rm -rf miseqsopdata.zip

cd $currpath

mv $output_dir/MiSeq_SOP/* $output_dir
rm -rf $output_dir/MiSeq_SOP

echo "Obtaining silva reference database" 
cd $db

#if [ ! -e silva.bacteria.zip ]; then
if [ ! -e silva.bacteria ]; then
        wget $silva
fi

echo "Decompressing and reformatting..."

unzip -u -o silva.bacteria.zip 
rm -rf silva.bacteria.zip
rm -rf __MACOSX

echo "Obtaining RDP reference database"

#if [ ! -e trainset9_032012.pds.tax ]; then
if [ ! -e trainset9_032012.pds.fasta ]; then
        wget $rdp
fi

echo "Decompressing and reformatting..."

unzip -u -o trainset9_032012.pds.zip 
rm -rf trainset9_032012.pds.zip

echo "Setting directories"

echo "Merging FASTQ files, creating stability.file"

#tell mothur which fastq files go together. We can do this with the make.file command. This command will use the text before the first _ of the fastq file names as the name of the sample. For this reason, it is best not to include - characters in your sample names

cd $currpath

singularity run $mothur mothur "#make.file(inputdir=$output_dir, outputdir=$output_dir, type=fastq, prefix=stability)"

#type parameter is used to indicate what file type you would like mothur to look for. Options are gz or fastq. Default=fastq.
#numcols parameter allows you to set number of columns you mothur to make in the file. Options 2 or 3. Default=3, meaning groupName forwardFastq reverseFastq
#prefix parameter allows you to enter your own prefix for the output filename. Default=stability

echo "Reducing sequencing and PCR errors"

#combine our two sets of reads for each sample and then to combine the data from all of the samples. This command will extract the sequence and quality score data from your fastq files, create the reverse complement of the reverse read and then join the reads into contigs.

singularity run $mothur mothur "#make.contigs(file=$output_dir/stability.files, outputdir=$output_dir)"

#This command will also produce several files that you will need down the road: stability.trim.contigs.fasta and stability.contigs.count_table. These contain the sequence data and group identity for each sequence. The stability.contigs.report file will tell you something about the contig assembly for each read.

#echo "Implement ambiguous bases removal" - not mandatory
#The following command implements the removal of any sequnece with ambiguous bases and any longer that 275 bp. 
#singularity run $mothur mothur "#screen.seqs(fasta=$output_dir/stability.trim.contigs.fasta, count=$output_dir/stability.contigs.count_table, maxambig=0, maxlength=275, maxhomop=8)"

#We anticipate that many of our sequences are duplicates of each other. Because it’s computationally wasteful to align the same thing a bazillion times, we’ll unique our sequences using the unique.seqs command:

echo "Processing improved sequences"

singularity run $mothur mothur  "#unique.seqs(fasta=$output_dir/stability.trim.contigs.fasta, count=$output_dir/stability.contigs.count_table, outputdir=$output_dir)"

echo "Align to reference database"

#Align our sequences to the reference alignment. Again we can make our lives a bit easier by making a database customized to our region of interest using the pcr.seqs command. To run this command you need to have the reference database (silva.bacteria.fasta) and know where in that alignment your sequences start and end. To remove the leading and trailing dots we will set keepdots to false. You could also run this command using your primers of interest. - not mandatory
#singularity run $mothur mothur "#pcr.seqs(fasta=$db/silva.bacteria/silva.bacteria.fasta, start=11895, end=25318, keepdots=F, outputdir=$db/silva.bacteria)" 

#Let’s rename it to something more useful using the rename.file command:

singularity run $mothur mothur "#rename.file(input=$db/silva.bacteria/silva.bacteria.fasta, new=silva.v4.fasta, outputdir=$db)"

#Now we have a customized reference alignment to align our sequences to. The nice thing about this reference is that instead of being 50,000 columns wide, it is now 13,425 columns wide which will save our hard drive some space and should improve the overall alignment quality. We’ll do the alignment with align.seqs:

singularity run $mothur mothur  "#align.seqs(fasta=$output_dir/stability.trim.contigs.unique.fasta, reference=$db/silva.v4.fasta, outputdir=$output_dir)"

# To make sure that everything overlaps the same region we’ll re-run screen.seqs to get sequences that start at or before position 1968 and end at or after position 11550. 

echo "Filtering overhangs"

# Now we know our sequences overlap the same alignment coordinates, we want to make sure they only overlap that region. So we’ll filter the sequences to remove the overhangs at both ends.

singularity run $mothur mothur "#filter.seqs(fasta=$output_dir/stability.trim.contigs.unique.align, vertical=T, trump=., outputdir=$output_dir)"
#singularity run $mothur mothur "#filter.seqs(fasta=$output_dir/stability.trim.contigs.unique.good.align, vertical=T, trump=., outputdir=$output_dir)"

# Because we’ve perhaps created some redundancy across our sequences by trimming the ends, we can re-run unique.seqs:

singularity run $mothur mothur "#unique.seqs(fasta=$output_dir/stability.trim.contigs.unique.filter.fasta, count=$output_dir/stability.trim.contigs.count_table, outputdir=$output_dir)"

echo "Denoising and pre-clustering"

# further de-noise our sequences is to pre-cluster the sequences using the pre.cluster command allowing for up to 2 differences between sequences. This command will split the sequences by group and then sort them by abundance and go from most abundant to least and identify sequences that are within 2 nt of each other. If they are then they get merged. We generally favor allowing 1 difference for every 100 bp of sequence

singularity run $mothur mothur "#pre.cluster(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.fasta, count=$output_dir/stability.trim.contigs.unique.filter.count_table, diffs=0, outputdir=$output_dir)"

#echo "Chimera check"
# using the VSEARCH algorithm that is called within mothur using the chimera.vsearch command. Again, this command will split the data by sample and check for chimeras. Our preferred way of doing this is to use the abundant sequences as our reference. In addition, if a sequence is flagged as chimeric in one sample, the default (dereplicate=F) is to remove it from all samples.
#singularity run $mothur mothur "#chimera.vsearch(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.fasta, count=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.count_table, dereplicate=t, outputdir=$output_dir)"

echo "Classifing using Bayesian classifier"

#classify those sequences using the Bayesian classifier with the classify.seqs command:

singularity run $mothur mothur "#classify.seqs(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.fasta, count=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.count_table, reference=$db/trainset9_032012.pds.fasta, taxonomy=$db/trainset9_032012.pds.tax, outputdir=$output_dir)"

#Now that everything is classified we want to remove our undesirables. We do this with the remove.lineage command:
#singularity run $mothur mothur "#remove.lineage(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.denovo.vsearch.fasta, count=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.denovo.vsearch.count_table, taxonomy=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.denovo.vsearch.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)"

#Now we have a couple of options for clustering sequences into OTUs. For a small dataset like this, we can do the traditional approach using dist.seqs and cluster:
#singularity run $mothur mothur "#dist.seqs(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.denovo.vsearch.pick.fasta, cutoff=0.03, outputdir=$output_dir)"

echo "Renaming files to make it smoother"

singularity run $mothur mothur "#rename.file(fasta=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.fasta, count=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.count_table, taxonomy=$output_dir/stability.trim.contigs.unique.filter.unique.precluster.pds.wang.taxonomy , prefix=final)"

echo "OTUs"

#use our cluster.split command. In this approach, we use the taxonomic information to split the sequences into bins and then cluster within each bin. The use of the cluster splitting heuristic was probably not worth the loss in clustering quality. However, as datasets become larger, it may be necessary to use the heuristic to clustering the data into OTUs. The advantage of the cluster.split approach is that it should be faster, use less memory, and can be run on multiple processors.

singularity run $mothur mothur "#cluster.split(fasta=$output_dir/final.fasta, count=$output_dir/final.count_table, taxonomy=$output_dir/final.taxonomy, taxlevel=5, cutoff=0.03, outputdir=$output_dir)"

#ow many sequences are in each OTU from each group and we can do this using the make.shared command. Here we tell mothur that we’re really only interested in the 0.03 cutoff level:

singularity run $mothur mothur "#make.shared(list=$output_dir/final.opti_mcc.list, count=$output_dir/final.count_table, label=0.03, outputdir=$output_dir)"

#We probably also want to know the taxonomy for each of our OTUs. We can get the consensus taxonomy for each OTU using the classify.otu command:

singularity run $mothur mothur "#classify.otu(list=$output_dir/final.opti_mcc.list, count=$output_dir/final.count_table, taxonomy=$output_dir/final.taxonomy, label=0.03, outputdir=$output_dir)"

echo "ASVs"

#OTUs generally represent sequences that are not more than 3% different from each other. In contrast, ASVs (aka ESVs) strive to differentiate sequences into separate OTUs if they are different from each other. The method built into mothur for identifying ASVs is pre.cluster. We did this above and then removed chimeras and contaminant sequences. We can convert the fasta and count_table files we used to form OTUs to a shared file using the make.shared command. This results in a shared and list file. The shared file we can use like the shared file from forming OTUs or phylotypes. The list file we can use to generate a consensus taxonomy for each ASV.

singularity run $mothur mothur "#classify.otu(list=$output_dir/final.asv.list, count=$output_dir/final.count_table, taxonomy=$output_dir/final.taxonomy, label=ASV, outputdir=$output_di)"

echo "DONE!"
echo "Excellent, your analysis with Mothur v.1.48  are completed. Have a nice day!"
