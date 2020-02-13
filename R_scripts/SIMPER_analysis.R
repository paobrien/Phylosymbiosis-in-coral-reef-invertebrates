## Use SIMPER analysis to identify ASVs that are contributing to incongruency between sites/species

# Paul A. O'Brien
# paul.obrien@my.jcu.edu.au

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

library(vegan)

taxonomy <- read.table("~/Documents/R/Ramaciotti_16S_dataset/taxonomy_silva.tsv", sep = "\t", header = T, strip.white = T) # For ASV IDs

## Coral ----
asv_tablec <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/coral-table-f.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
meta_datac <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Metadata/coral_metadata.tsv", sep = '\t', header = T, strip.white = T)

asv_tablec <- as.data.frame(t(asv_tablec))

#SIMPER
simp.testc <- simper(asv_tablec, meta_datac$Species)
sum.simpc <- summary(simp.testc)

# Pull out the comparison you want to look at
pcyl <- sum.simpc$`P. cylindrica_R_P. cylindrica_P` 
shyst <- sum.simpc$`S. hysterix_R_S. hysterix_P`

write.csv(pcyl, file = "p.cylindra_simp.csv")
write.csv(shyst, file = "s.hyst_simp.csv")

## Octocoral ----
asv_tableo <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/octocoral-table-f.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
meta_datao <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Metadata/octocoral_metadata.tsv", sep = '\t', header = T, strip.white = T)
asv_tableo <- as.data.frame(t(asv_tableo))

#SIMPER
simp.testo <- simper(asv_tableo, meta_datao$Species)
sum.simpo <- summary(simp.testo)

#pull out the comparison you want to look at
sinul <- sum.simpo$`Sinularia spR_Sinularia spP`
sarc <- sum.simpo$`Sarcophyton spR_Sarcophyton spP`
gorg <- sum.simpo$`I. hippurus_Pinnigorgia sp`

write.csv(sinul, file = "sin_simp.csv")
write.csv(sarc, file = "sarc.csv")
write.csv(gorg, file = "gorg.csv")

