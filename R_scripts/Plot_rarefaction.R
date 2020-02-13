## Plot rarefaction curve for phylosymbiosis MS

# Paul A. O'Brien
# paul.obrien@my.jcu.edu.au

# Note: this script uses the function ggrare, which can be found at:
# https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/R/

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis")

## Load packages
library(ggplot2)
library(phyloseq)
library(tidyr)
library(reshape2)


#Import data
asv_table <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/table-f-NoGor.txt", 
                        sep = '\t', row.names = 1, header = T, strip.white = T)
asv.mat <- as.matrix(asv_table)
taxonomy <- read.table("~/Documents/R/Ramaciotti_16S_dataset/taxonomy_silva.tsv", 
                       sep = '\t', row.names = 1, header = T, strip.white = T)
taxonomy$Confidence <- NULL
metadata <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Metadata/metadata_noGor.tsv", 
                       sep = '\t', row.names = 1, header = T, strip.white = T)


# Separate taxonomy into different columns
tax_sep <- separate(taxonomy, Taxon, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
                    sep = ";", remove = TRUE, convert = FALSE, extra = "warn", fill = "warn")

# Replace all NA's
tax_sep$Phylum[is.na(tax_sep$Phylum)] <- "p__"
tax_sep$Class[is.na(tax_sep$Class)] <- "c__"
tax_sep$Order[is.na(tax_sep$Order)] <- "o__"
tax_sep$Family[is.na(tax_sep$Family)] <- "f__"
tax_sep$Genus[is.na(tax_sep$Genus)] <- "g__"
tax_sep$Species[is.na(tax_sep$Species)] <- "s__"

# 'Phyloseq-ize' the data
asv.tab <- otu_table(as.matrix(asv_table[-grep("control", colnames(asv_table), ignore.case = T)]), # also remove zymo control samples
                     taxa_are_rows = T)

tax.tab <- tax_table(as.matrix(tax_sep))

sample.data <- sample_data(metadata[-grep("control", row.names(metadata),ignore.case = T),]) # also remove zymo control samples

physeq <- phyloseq(asv.tab, tax.tab, sample.data)
physeq

# Remove unwanted samples
physeq.sub <- subset_samples(physeq, Taxonomy!="Mollusc")

#Change order
sample_data(physeq.sub)$Taxonomy = factor(sample_data(physeq.sub)$Taxonomy,
                                          levels = c("Seawater","Coral","Soft Coral","Gorgonian",
                                                     "Ascidian", "Sponge", "Blank"))

# Plot
cbPalette <- c("Seawater" = "#56B4E9", "Coral" = "#009E73", "Soft Coral" = "#F0E442", "Gorgonian" = "#0072B2",  
               "Ascidian" = "#D55E00", "Sponge" = "#CC79A7", "Blank" = "#999999") # colourblind friendly colour palette

rarefaction_plot <- ggrare(physeq.sub, step = 500, color = "Taxonomy", label = NULL, se = FALSE) #without

rarefaction_plot_f <- rarefaction_plot + facet_wrap(~Taxonomy) + 
  scale_color_manual(values = cbPalette) + guides(colour = guide_legend(title = "Host")) + 
  scale_x_continuous(breaks = c(3500, 30000, 60000, 90000)) +
  theme_bw()
rarefaction_plot_f

ggsave("Rarefaction_final.pdf", width = 10, height = 8, dpi = "print")




