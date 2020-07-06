#!/bin/bash
export LANG=en_US.utf-8
export LC_ALL=en_US.utf-8
export MPLBACKEND='Agg'

module load contrib/anaconda/anaconda4.4.0
source activate /sw/contrib/qiime2/2019.1

## import demultiplexed fastq files
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]'  --input-path fastq_manifest_novaseq2.csv --input-format PairedEndFastqManifestPhred33 --output-path gelada-paired-end-demux-novaseq2.qza

## summarize and visualize
qiime demux summarize --i-data gelada-paired-end-demux-novaseq2.qza --o-visualization gelada-demux-novaseq2.qzv

