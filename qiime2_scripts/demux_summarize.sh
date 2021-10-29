#!/bin/sh

dir=$HOME/MILKQUA/Analysis

qiime demux summarize \
  --i-data $dir/paired-end-demux.qza \
  --o-visualization $dir/demux_seqs.qzv


