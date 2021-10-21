
qiime tools import  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $HOME/MILKQUA/data/subset \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path $HOME/MILKQUA/Analysis/demux-paired-end.qza

