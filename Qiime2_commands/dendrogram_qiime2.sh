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

### ------------------------------------------ ###


# Commands used to create dendrograms used for the analysis in the manuscript on phylosymbiosis in coral reef invertebrates
# Paul A. O'Brien
# paul.obrien@my.jcu.edu.au


module load qiime2/2018.4


###################
### Filter data ###
###################


# Remove Chloroplast, Mitochondria and Eukarytoic sequences from data

qiime taxa filter-table \
  --i-table Output_files/Feature_tables/table-dada2.qza \
  --i-taxonomy Output_files/Taxonomy_silva/taxonomy_silva.qza \
  --p-exclude Mitochondria,Chloroplast,Eukaryota \
  --o-filtered-table table-dada2-f.qza

# Remove unknown gorgonian samples
qiime feature-table filter-samples \
  --i-table table-dada2-f.qza \
  --m-metadata-file Metadata_files/metadata.tsv \
  --p-where "\"sample-id\" NOT IN ('Gor1','Gor2','Gor3')" \
  --o-filtered-table table-NoG-f.qza

# Remove blanks and seawater samples as irrelevant to host phylogeny
qiime feature-table filter-samples \
  --i-table table-NoG-f.qza \
  --m-metadata-file Metadata_files/metadata.tsv \
  --p-where "\"Taxonomy\" NOT IN ('Blank','Seawater')" \
  --o-filtered-table table-NoGBS-f.qza


### --------------------------------------------------- ###


# Subset ASV table to host group of interest, including outgroup species
# Then...
# Filter ASVs that are only present 2 times or less
# Filter ASVs that are only in one sample


### Coral hosts plus outgroup ###


qiime feature-table filter-samples \
  --i-table table-NoGBS-f.qza \
  --m-metadata-file coral_metadata.tsv \  # Create metadata file with samples to keep
  --o-filtered-table coral-only-table-f.qza

qiime feature-table filter-features \
  --i-table coral-only-table-f.qza \
  --p-min-frequency 3 \
  --p-min-samples 2 \
  --o-filtered-table table-coral-no-doubles-min-2-f.qza


### Repeat for other host taxa ###


#########################
### Create dendrogram ###
#########################


### Coral ###

# Group samples by host species
qiime feature-table group \
  --i-table Output_files/Feature_tables/coral-table-no-doubles-min-2-f.qza \
  --p-axis sample \
  --m-metadata-file Metadata_files/coral_metadata.tsv \
  --m-metadata-column Species \
  --p-mode mean-ceiling \
  --o-grouped-table table-grouped-coral-f.qza

# Cluster species into dendrogram
qiime diversity beta-rarefaction \
  --i-table table-grouped-coral-f.qza \
  --p-metric braycurtis \
  --p-clustering-method upgma \
  --p-iterations 1000 \
  --m-metadata-file Metadata_files/coral_grouped_metadata.tsv \  # Create grouped metadata file
  --p-sampling-depth 14000 \  # Lowest number of sequences in grouped sample
  --o-visualization coral_den_grouped_bc_upgma.qzv

# Dendrogram file can be exported from the .qzv file at https://view.qiime2.org/

### Repeat for other hosts ###


