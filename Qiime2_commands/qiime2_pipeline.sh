#!/bin/sh

#PBS -l nodes=1
#PBS -l walltime=24:0:0
## #PBS -A  pobrien

### - job basename ---------------------------------------------------------
jobname=qiime2_16S
### ------------------------------------------------------------------------

## #PBS -S /bin/bash
cd $PBS_O_WORKDIR
export TMPDIR=/export/scratch/pobrien
TD=$(mktemp -d --tmpdir=/export/scratch/pobrien)


### ---------------------------------------- ###


########################
### QIIME 2 ANALYSIS ###
########################


# Commands used to create ASV table and taxonomic table used for the analysis in the manuscript on phylosymbiosis in coral reef invertebrates
# Paul A. O'Brien
# paul.obrien@my.jcu.edu.au


### ---------------------------- ###


module load qiime2/2018.4

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.txt \
  --output-path paired-end-demux.qza \
  --source-format PairedEndFastqManifestPhred33


qiime demux summarize \
  --i-data paired-end-demux.qza \
  --o-visualization demux.qzv

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs paired-end-demux.qza \
  --p-trunc-len-f 248 \
  --p-trunc-len-r 231 \
  --p-trim-left-f 12 \
  --p-trim-left-r 12 \
  --p-n-threads 4 \
  --o-table table-dada2.qza \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-denoising-stats stats-dada2.qza 
   
qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

qiime feature-table summarize \
  --i-table table-dada2.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs-dada2.qza \
  --o-visualization rep-seqs.qzv

qiime alignment mafft \
  --i-sequences rep-seqs-dada2.qza \
  --o-alignment aligned-rep-seqs.qza 
  
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza

qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza 

qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza  

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table Output_files/table-dada2.qza \
  --p-sampling-depth 3500 \
  --m-metadata-file metadata.tsv \
  --output-dir Core-metrics-results-dada2

### repeat with filtered table (see filter_data script) ---------

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny ../Output_files/Tree_files/Rooted-tree/rooted-tree.qza \
  --i-table ../Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Table-dada2-filtered/table-dada2-filtered.qza \
  --p-sampling-depth 3300 \
  --m-metadata-file ../Metadata_files/metadata.tsv \
  --output-dir ../Core-metrics-results-dada2-filtered

### -----------	

qiime diversity alpha-group-significance \
  --i-alpha-diversity Core-metrics-results-dada2/faith_pd_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization faith-pd-group-significance_dada2.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity Core-metrics-results-dada2/shannon_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization shannon_vector.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity Core-metrics-results-dada2/observed_otus_vector.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization observed_otus_vector.qzv

qiime diversity alpha-rarefaction \
  --i-table table-dada2.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 52000 \
  --p-min-depth 2000 \
  --m-metadata-file metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix Core-metrics-results-dada2/bray_curtis_distance_matrix.qza \
  --m-metadata-file metadata.tsv \
  --m-metadata-column Taxonomy \
  --o-visualization Core-metrics-results-dada2/bray_curtis_taxonomy_significance.qzv \
  --p-pairwise

## download the feature classifier for taxonomic assignment. Same primers as used for our data

wget \
  -O "gg-13-8-99-515-806-nb-classifier.qza" \
  "https://data.qiime2.org/2018.8/common/gg-13-8-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

## Repeat taxonomic assignment with silva DB

 wget -O "silva_132_99_515F_806R_nb_classifier.qza" "https://data.qiime2.org/2018.8/common/silva-132-99-515-806-nb-classifier.qza"

qiime feature-classifier classify-sklearn \
  --i-classifier silva_132_99_515F_806R_nb_classifier.qza \
  --i-reads rep-seqs-dada2.qza \
  --o-classification taxonomy_silva.qza

qiime metadata tabulate \
  --m-input-file taxonomy_silva.qza \
  --o-visualization taxonomy_silva.qzv



