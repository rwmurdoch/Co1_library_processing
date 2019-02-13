#!/bin/bash

#########################################
## Adapting the Yongchao script to ###
## willcox's data
## Jan 18, 2019
#########################################

## conda activate qiime2-2018.11 ##

## data import ##

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path reads \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path willcoxreads.qza

qiime demux summarize \
--i-data willcoxreads.qza \
--o-visualization willcox-demux.qzv \
--verbose


## trim the primers (Zeale)
## stdout is directed to a logfile

qiime cutadapt trim-paired \
--i-demultiplexed-sequences willcoxreads.qza \
--p-cores 8 \
--p-front-f AGATATTGGAACWTTATATTTTATTTTTGG \
--p-front-r WACTAATCAATTWCCAAATCCTCC \
--o-trimmed-sequences willcox-data-cutadapt \
--verbose > cutadapt_log.txt

###############################################
# training the feature classifier
###############################################

## ONLY RUN THIS ONE TIME, IT TAKES HOURS##
## These scripts are not relevant to this dataset!
## Contact me (RWMurdoch) for info on how to deal with the BOLD classifiers

## I'll work with the  full-length silva132 database (seqs and taxonomy)
## downloaded from: https://forum.qiime2.org/t/silva-132-classifiers/3698/4 (distributed by Caporaso)
## Database is first trimmed using the same primers used for the study (341F/785R)

#qiime feature-classifier extract-reads \
#  --i-sequences 99-otus.qza \
#  --p-f-primer CCTACGGGNGGCWGCAG \
#  --p-r-primer GACTACHVGGGTATCTAATCC \
#  --o-reads ref-seqs.qza \
#  --verbose

#qiime feature-classifier fit-classifier-naive-bayes \
#--i-reference-reads ref-seqs.qza \
#--i-reference-taxonomy 7_level_taxonomy.qza \
#--o-classifier classifier.qza


#####################################
## dada2 and classification
####################################

# The PCR product was ~157 bp.  Manual truncation to 140bp
# Leaves large overlap without large risk of readthrough contamination

qiime dada2 denoise-paired \
--i-demultiplexed-seqs willcox-data-cutadapt.qza \
--o-table willcox-dada2-table \
--o-representative-sequences willcox-dada2_seqs \
--o-denoising-stats willcox-dada2-stats \
--p-n-threads 8 \
--p-trunc-len-f 140 \
--p-trunc-len-r 140

###########################################
# Size filtering of reads/ASVs #
#########################################

# this is a custom step which is not built into qiime2
# basic strategy is to export the dada2_seqs, apply a size filter in R
# and re-import under the same name

qiime tools export \
  --input-path willcox-dada2_seqs.qza \
  --output-path unfiltered.ASV.reps.fa

#upper and lower size limits must be altered in the script directly
Rscript size.filter.seqs.q2.R

qiime tools import \
  --type FeatureData[Sequence] \
  --input-path new.seqs.fa \
  --output-path willcox-dada2_seqs

qiime feature-table summarize \
--i-table willcox-dada2-table.qza \
--o-visualization willcox-dada2-seq-stats \
--verbose



#not enough ram (64Gb system) to use default setting for reads-per-batch

qiime feature-classifier classify-sklearn \
  --i-classifier BOLD.classifier.6.2017.qza \
  --i-reads willcox-dada2_seqs.qza \
  --p-n-jobs -1 \
  --p-reads-per-batch 10 \
  --o-classification willcox-dada2-taxonomy

qiime metadata tabulate \
  --m-input-file willcox-dada2-taxonomy.qza \
  --o-visualization willcox-dada2-taxonomy

#remove singletons
qiime feature-table filter-features \
  --i-table willcox-dada2-table.qza \
  --p-min-frequency 2 \
  --o-filtered-table willcox-dada2-table-filtered



####################################
# taxonomy visualization
###############################

qiime taxa barplot \
  --i-table willcox-dada2-table-filtered.qza \
  --i-taxonomy willcox-dada2-taxonomy.qza \
  --m-metadata-file willcox_metadata.csv \
  --o-visualization willcox-dada2-taxa-bar-plots

qiime feature-table summarize \
  --i-table willcox-dada2-table-filtered.qza \
  --o-visualization willcox-dada2-table.qzv \
  --m-sample-metadata-file willcox_metadata.csv

qiime feature-table tabulate-seqs \
  --i-data willcox-dada2_seqs.qza \
  --o-visualization willcox-dada2_seqs.qzv

#######################################
# preparation for alpha div metrics
#######################################

qiime alignment mafft \
  --i-sequences willcox-dada2_seqs.qza \
  --o-alignment willcox-dada2_seqs-aligned.qza \
  --p-n-threads 8

qiime alignment mask \
 --i-alignment willcox-dada2_seqs-aligned.qza \
--o-masked-alignment willcox-dada2_seqs-aligned-masked.qza

qiime phylogeny fasttree \
 --i-alignment willcox-dada2_seqs-aligned-masked.qza \
 --o-tree willcox-unrooted-tree-dada2.qza

qiime phylogeny midpoint-root \
--i-tree willcox-unrooted-tree-dada2.qza \
--o-rooted-tree willcox-rooted-tree-dada2.qza

################################################
# generate core diversity metrics #
# pay close attention to the sampling depth command #
###################################################

#qiime diversity core-metrics-phylogenetic \
#  --i-phylogeny willcox-rooted-tree-dada2.qza \
#  --i-table willcox-dada2-table-filtered.qza \
#  --p-sampling-depth 30000 \
#  --m-metadata-file willcox_metadata.csv \
#  --output-dir willcox-core-metrics-results-dada2

## visualization of diversity metrics ##

#qiime diversity alpha-rarefaction \
#  --i-table willcox-dada2-table-filtered.qza \
#  --i-phylogeny willcox-rooted-tree-dada2.qza \
#  --p-min-depth 1 \
#  --p-max-depth 10000 \
#  --m-metadata-file willcox_metadata.csv \
#  --o-visualization willcox-open-alpha-rarefaction-dada2.qzv

###############################################
# exporting data
# you can easily export seqs and taxonomy directly out of
# their corresponding .qza files using the qiime tools export feature.
# to get the feature table, you export the table in biom format and then convert
# into a simple text file
#################################

qiime tools export --input-path willcox-dada2-table-filtered.qza --output-path exports
qiime tools export --input-path willcox-dada2_seqs.qza --output-path exports
qiime tools export --input-path willcox-dada2-taxonomy.qza --output-path exports

biom convert -i exports/feature-table.biom \
-o exports/otu_table.txt --to-tsv

## combine everything into a single table
## this uses the R script create.OTU.table.R
## script is available here: https://github.com/rwmurdoch/project.scripts/blob/master/create.OTU.table.R
## download/copy the script into your project directory; it will read and write in the "export" directory

Rscript create.OTU.table.R

## combine various outputs into a results package

tar -cf results.tar.gz \
exports \
willcox-dada2-taxa-bar-plots.qzv \
willcox_metadata.csv \
willcox-dada2-seq-stats.qzv \
willcox-demux.qzv

#willcox-core-metrics-results-dada2 \
#willcox-open-alpha-rarefaction-dada2.qzv \
