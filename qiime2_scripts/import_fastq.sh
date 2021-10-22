
#qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path $HOME/MILKQUA/data/subset \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --output-path $HOME/MILKQUA/Analysis/demux-paired-end.qza

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $HOME/MILKQUA/Config/subset_manifest.tsv \
  --output-path $HOME/MILKQUA/Analysis/paired-end-demux.qza \
  --input-format PairedEndFastqManifestPhred33V2
