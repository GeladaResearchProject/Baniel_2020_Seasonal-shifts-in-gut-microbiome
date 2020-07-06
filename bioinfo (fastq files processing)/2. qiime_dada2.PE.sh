#!/bin/bash
export LANG=en_US.utf-8
export LC_ALL=en_US.utf-8
export MPLBACKEND='Agg'

module load contrib/anaconda/anaconda4.4.0
source activate /sw/contrib/qiime2/2019.1

## Trim reads and denoise based on the visualization from above (DADA2 denoise)
qiime dada2 denoise-paired --i-demultiplexed-seqs gelada-paired-end-demux-novaseq2.qza --p-n-threads 24 --p-trim-left-f 20 --p-trim-left-r 20 --p-trunc-len-f 220 --p-trunc-len-r 180 --o-table gelada-table.qza --o-representative-sequences gelada-rep-seqs.qza --p-max-ee 2 --verbose --p-trunc-q 2 --o-denoising-stats gelada-stats-dada2.qza

qiime metadata tabulate --m-input-file gelada-stats-dada2.qza --o-visualization gelada-stats-dada2.qzv

qiime feature-table summarize --i-table gelada-table.qza --o-visualization gelada-table.qzv --m-sample-metadata-file gelada-sample-metadata.txt

qiime feature-table tabulate-seqs --i-data gelada-rep-seqs.qza --o-visualization gelada-rep-seqs.qzv

qiime alignment mafft --i-sequences gelada-rep-seqs.qza --o-alignment gelada-aligned-rep-seqs.qza --p-n-threads 24

qiime alignment mask --i-alignment gelada-aligned-rep-seqs.qza --o-masked-alignment gelada-masked-aligned-rep-seqs.qza

qiime phylogeny fasttree --i-alignment gelada-masked-aligned-rep-seqs.qza --o-tree gelada-unrooted-tree.qza --p-n-threads 24

qiime phylogeny midpoint-root --i-tree gelada-unrooted-tree.qza --o-rooted-tree gelada-rooted-tree.qza

#alpha diversity (check sampling depth)
qiime diversity alpha-rarefaction \
  --i-table gelada-table.qza \
  --i-phylogeny gelada-rooted-tree.qza \
  --p-max-depth 50000 \
  --m-metadata-file gelada-sample-metadata.txt \
  --o-visualization gelada-alpha-rarefaction.qzv
  
qiime diversity core-metrics-phylogenetic --i-phylogeny gelada-rooted-tree.qza --i-table gelada-table.qza --p-sampling-depth 5000 --m-metadata-file gelada-sample-metadata.txt --output-dir gelada-core-metrics-results

qiime feature-classifier classify-sklearn \
  --i-classifier silva-132-99-515-806-nb-classifier.qza \
  --i-reads gelada-rep-seqs.qza \
  --o-classification gelada-taxonomy-silva.qza
  
qiime metadata tabulate \
  --m-input-file gelada-taxonomy-silva.qza \
  --o-visualization gelada-taxonomy-silva.qzv
  
qiime taxa barplot \
  --i-table gelada-table.qza \
  --i-taxonomy gelada-taxonomy-silva.qza \
  --m-metadata-file gelada-sample-metadata.txt \
  --o-visualization gelada-taxa-bar-plots-silva.qzv

qiime feature-table filter-samples \
--i-table gelada-table.qza \
--p-min-features 1 \
--o-filtered-table gelada-table-filter.qza
