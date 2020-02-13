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


# Commands used to create dendrograms for phylosymbiosis analysis in coral reef invertebrates
# Paul A. O'Brien
# paul.obrien@my.jcu.edu.au


module load qiime2/2018.4


###################
### Filter data ###
###################


# Subset ASV table to host group of interest, including outgroup species
# Filter ASVs that are only present 2 times or less
# Filter ASVs that are only in one sample


### Coral hosts plus outgroup ###


qiime feature-table filter-samples \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/table-NoGBS-f.qza \
  --m-metadata-file /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Metadata_files/coral_metadata.tsv \  # Create metadata file with samples to keep
  --o-filtered-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-only-table-f.qza

qiime feature-table filter-features \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-only-table-f.qza \
  --p-min-frequency 3 \
  --p-min-samples 2 \
  --o-filtered-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/table-coral-no-doubles-min-2-f.qza


### Repeat for other host taxa ###


######################
# Fragment Insertion #
######################

module load qiime2/2019.1

## Perform fragment insertion of 16S sequences using the SEPP algorithm against the Greengenes 13_8 99% tree.
#  Use tree for UniFrac analysis

wget \
  -O "sepp-refs-silva-128.qza" \
  "https://data.qiime2.org/2019.10/common/sepp-refs-silva-128.qza"

qiime fragment-insertion sepp \
  --i-representative-sequences /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Rep_seqs/Rep-seqs-dada2/rep-seqs-dada2.qza \
  --o-tree /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Tree_files/insertion-tree.qza \
  --o-placements /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Tree_files/insertion-placements.qza

qiime fragment-insertion filter-features \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Table-dada2-filtered/table-dada2-filtered.qza \
  --i-tree /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Tree_files/Fragment_insertion/insertion-tree.qza \
  --o-filtered-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/tree-filtered_table.qza \
  --o-removed-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/tree-removed_table.qza

qiime tools export \
  --input-path /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/tree-filtered_table.qza \
  --output-path /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Tree-filtered

qiime tools export \
  --input-path /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/tree-removed_table.qza \
  --output-path /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Tree-removed



#########################
### Create dendrogram ###
#########################

module load qiime2/2018.4

### Coral ###

# Group samples by host species
qiime feature-table group \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/table-coral-no-doubles-min-2-f.qza \
  --p-axis sample \
  --m-metadata-file /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Metadata_files/coral_metadata_wCar.tsv \
  --m-metadata-column Species \
  --p-mode mean-ceiling \
  --o-grouped-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-wCar-no-doubles-min-2-grouped-f.qza 


# Cluster species into dendrogram

# Bray-Curtis 
qiime diversity beta-rarefaction \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-wCar-no-doubles-min-2-grouped-f.qza \
  --p-metric braycurtis \
  --p-clustering-method upgma \
  --p-iterations 1000 \
  --m-metadata-file /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Metadata_files/coral_grouped_metadata_wCar.tsv \ # Create grouped metadata file
  --p-sampling-depth 10000 \  # Lowest number of sequences in grouped sample
  --o-visualization coral_den_grouped_bc_upgma.qzv

# Weighed UniFrac
qiime diversity beta-rarefaction \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-wCar-no-doubles-min-2-grouped-f.qza \
  --i-phylogeny /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Tree_files/Fragment_insertion/insertion_tree_trimmed.qza \
  --p-metric weighted_unifrac \
  --p-clustering-method upgma \
  --p-iterations 1000 \
  --m-metadata-file /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Metadata_files/coral_grouped_metadata_wCar.tsv \
  --p-sampling-depth 10000 \
  --o-visualization coral_den_FI-tree_trim_wuf_upgma.qzv

qiime diversity beta-rarefaction \
  --i-table /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Filtered_no_mitochondira_chloroplast_eukaryota/Feature_tables/Coral/coral-wCar-no-doubles-min-2-grouped-f.qza \
  --i-phylogeny /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Output_files/Tree_files/Fragment_insertion/insertion_tree_trimmed.qza \
  --p-metric unweighted_unifrac \
  --p-clustering-method upgma \
  --p-iterations 1000 \
  --m-metadata-file /export/home/l-p/pobrien/PhD/Ramaciotti_OBR5432-95032942/Metadata_files/coral_grouped_metadata_wCar.tsv \
  --p-sampling-depth 10000 \
  --o-visualization coral_den_FI-tree_trim_uwuf_upgma.qzv



### Repeat for other hosts ###

