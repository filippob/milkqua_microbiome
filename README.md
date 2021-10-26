# milkqua_microbiome
Repo to store code for the analysis of 16S rRNA-gene sequencing data in the MILKQUA project.

Initially, the following tools will be explored and used for training:
- Qiime2
- Conda environments
- Singularity
- Nextflow

## how do Qiime 2 steps relate to Qiime 1.9 steps?

From Qiime 1.9:

1. extract barcodes
2. join paried-end reads
3. assemble a single Fasta file with all clean/filtered sequences (mutliple_split_library.py)
4. OTU picking (usually closed OTU picking)

Qiime 2 steps include:

1. import data
2. 
