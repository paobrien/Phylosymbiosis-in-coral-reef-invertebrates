## Calculate and plot Relative Abundance for microbiome data

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

library(ggplot2)
library(forcats)
library(reshape2)
library(dplyr) 

## Plot phylum level ----
asv_table_p <- read.csv("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/phylyum_tax_f.csv", 
                      sep = ',', header = T, row.names = 1, strip.white = T) # phylum table created using qiime2


# Remove unwanted samples
asv_table_p <- asv_table_p[-grep("control", row.names(asv_table_p), ignore.case = T),] # sequencing controls
asv_table_p <- asv_table_p[!(asv_table_p$Taxonomy=="Mollusc"),] # Mollusc removed from analysis as only one species


#subset metadata 
meta_data <- subset(asv_table_p, select=c("Taxonomy", "Species", "Reproduction", # first create metadata table
                                        "Reef", "Zone", "Date")) 

asv_table_p2 <- as.matrix(asv_table_p[, -c(63:68)]) # Remove metadata columns


# Turn absolute values to a percentage. Divide value by sum of row, multiply by 100
asv_percent <- asv_table_p2 / rowSums(asv_table_p2) * 100

# Join with metadata
asv_rel_abun <- cbind(asv_percent, meta_data)

# reshape to long format
asv_table_long <- melt(asv_rel_abun, variable.name = "Phyla")
head(asv_table_long, 20)

# Remove long names
unique(asv_table_long$Phyla)

asv_table_long$Phyla <- gsub("D_0__Archaea.D_1__", "\\1", asv_table_long$Phyla) # the \\1 sets the replacement string to that capture group
asv_table_long$Phyla <- gsub("Unassigned.__", "Unassigned", asv_table_long$Phyla)
asv_table_long$Phyla <- gsub("D_0__Archaea.__", "Unknown Archaea", asv_table_long$Phyla)
asv_table_long$Phyla <- gsub("D_0__Bacteria.D_1__", "\\1", asv_table_long$Phyla)
asv_table_long$Phyla <- gsub("D_0__Bacteria.__", "Unknown Bacteria", asv_table_long$Phyla)
asv_table_long$Phyla <- gsub("D_0__Bacteria.__", "Unknown Bacteria", asv_table_long$Phyla)
asv_table_long$Phyla <- gsub("Marinimicrobia..SAR406.clade.", "SAR406", asv_table_long$Phyla)

# Group by taxonomy and summarise
group_taxon <- asv_table_long %>% 
  mutate(Taxonomy = fct_relevel(Taxonomy, "Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge", "Blank")) %>%
  group_by(Taxonomy, Phyla)
group_taxon

taxon_summary <- summarise(group_taxon, Percentage = mean(value, na.rm = T))
taxon_summary

# Plot relative abundance
bubble <- ggplot(taxon_summary, aes(x = Taxonomy, y = Phyla)) +
  geom_point(aes(size = Percentage), colour = "#000000", fill = "#56B4E9", shape = 21) + 
  scale_size(range = range(taxon_summary$Perentage)) +
  scale_size_area(max_size = 10) + 
  labs(x = "\nHost", y = "Relative Abundance (Phyla)\n") +
  guides(colour = guide_legend(title = "Host")) +
  theme_bw() +
  theme(axis.text=element_text(size=10),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
bubble

ggsave(filename = "Rel_Abun_Phyla_final.pdf", device = "pdf", width = 27, height = 20, units = "cm", dpi = "print")

## Plot family level (Top 25) ----
asv_table_f <- read.csv("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/family_tax_f.csv", sep = ',', header = T, row.names = 1, strip.white = T)

# Remove unwanted samples 
asv_table_f <- asv_table_f[-grep("control", row.names(asv_table_f), ignore.case = T),] # Remove sequencing controls
asv_table_f <- asv_table_f[!(asv_table_f$Taxonomy=="Mollusc"),] 

# Subset metadata
meta_data <- subset(asv_table_f, select=c("Taxonomy", "Species", "Reproduction",
                                        "Reef", "Zone", "Date")) 

asv_table_f2 <- as.matrix(asv_table_f[, -c(1228:1233)])


# Calculate Top 25 abundant family #

# Calculcate column totals
total <- colSums(asv_table_f2)
asv_w_total <- rbind(asv_table_f2, total)

# re-order columns by total
asv_sort <- asv_w_total[,order(-asv_w_total[which(rownames(asv_w_total) == 'total'), ])]

## Convert counts to percentages on data matrix & remove row total
asv_percent <- asv_sort[-which(rownames(asv_sort) == 'total'),] / rowSums(asv_sort[-which(rownames(asv_sort) == 'total'),]) * 100

# create new column "other" which is the sum of all the taxa not in top 25, remove other taxa
Other <- rowSums(asv_percent[,-c(1:25)])

# Combine data
asv_top25 <- cbind(asv_percent[,1:25], Other)
asv_rel_abunF <- cbind(asv_top25, meta_data)
asv_table_longF <- melt(asv_rel_abunF, variable.name = "Family")

# Remove long names
unique(asv_table_longF$Family)

asv_table_longF$Family <- gsub("D_0__Ba.*D_4__", "\\1", asv_table_longF$Family) 
asv_table_longF$Family <- gsub("D_0__Ba.*D_2__", "\\1", asv_table_longF$Family)
asv_table_longF$Family <- gsub("Alphaproteobacteria.__.__", "Alphaproteobacteria unknown ", asv_table_longF$Family)
asv_table_longF$Family <- gsub("D_0__Ar.*D_4__", "\\1", asv_table_longF$Family)
asv_table_longF$Family <- gsub("Gammaproteobacteria.__.__", "Gammaproteobacteria unknown", asv_table_longF$Family)
asv_table_longF$Family <- gsub("uncultured.alpha.proteobacterium", "Alphaproteobacteria uncultured 2", asv_table_longF$Family)
asv_table_longF$Family <- gsub("D_0__Bacteria.__.__.__.__", "Bacteria unknown", asv_table_longF$Family)
asv_table_longF$Family <- gsub("Unassigned.__.__.__.__", "Unassigned", asv_table_longF$Family)
asv_table_longF$Family <- gsub("Alphaproteobacteria.D_3__uncultured.__", "Alphaproteobacteria uncultured", asv_table_longF$Family)
asv_table_longF$Family <- gsub("uncultured.bacterium", "Bacteria uncultured", asv_table_longF$Family)

# Group by Taxonomy
group_taxonF <- asv_table_longF %>% 
  mutate(Taxonomy = fct_relevel(Taxonomy, "Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge", "Blank")) %>%
  group_by(Taxonomy, Family)
group_taxonF

taxon_summaryF <- summarise(group_taxonF, Percentage = mean(value, na.rm = T))
taxon_summaryF <- taxon_summaryF %>% mutate(Family = fct_relevel(Family, "Other", after = Inf))

# Plot family relative abundance
bubblef <- ggplot(taxon_summaryF, aes(x = Taxonomy, y = Family)) +
  geom_point(aes(size = Percentage), colour = "#000000", fill = "#56B4E9", shape = 21) + 
  scale_size(range = range(taxon_summaryF$Percentage)) +
  scale_size_area(max_size = 18) +
  labs(x = "\nHost", y = "Relative Abundance (Family)\n") +
  guides(colour = guide_legend(title = "Host")) +
  theme_bw() +
  theme(axis.text=element_text(size=11, colour = "black"),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
bubblef


## Plot subset by host taxonomy ----
asv_table_a <- asv_table_f[which(asv_table_f$Taxonomy=="Ascidian" | asv_table_f$Species=="C. foliascens"),]
asv_table_c <- asv_table_f[which(asv_table_f$Taxonomy=="Coral" | asv_table_f$Species=="C. foliascens"),]
asv_table_s <- asv_table_f[which(asv_table_f$Taxonomy=="Sponge" | asv_table_f$Species=="Cladiella sp"),]
asv_table_o <- asv_table_f[which(asv_table_f$Taxonomy=="Soft Coral" | asv_table_f$Taxonomy=="Gorgonian" | asv_table_f$Species=="C. foliascens"),]

asv_tables <- list(asv_table_a, asv_table_c, asv_table_o, asv_table_s)

metadata_all <- list()
for (i in seq_along(asv_tables)){
  metadata_all[[i]] <- subset(asv_tables[[i]], select=c("Species", "Reproduction", "Reef", "Zone", "Date")) 
}

matrix_all <- list()
for (i in seq_along(asv_tables)) {
  matrix_all[[i]] <- as.matrix(asv_tables[[i]][,-c(1228:1233)])
}

# Get column totals
total_a <- colSums(matrix_all[[1]])
total_c <- colSums(matrix_all[[2]])
total_o <- colSums(matrix_all[[3]])
total_s <- colSums(matrix_all[[4]])

asv_a <- rbind(matrix_all[[1]], total_a)
asv_c <- rbind(matrix_all[[2]], total_c)
asv_o <- rbind(matrix_all[[3]], total_o)
asv_s <- rbind(matrix_all[[4]], total_s)

asv_tables_total <- list(asv_a, asv_c, asv_o, asv_s)

# re-order columns by total
for (i in seq_along(asv_tables_total)){
  rownames(asv_tables_total[[i]]) <- gsub("total_.*", "total", rownames(asv_tables_total[[i]]))
}

asv_sort_all <- list()
for (i in seq_along(asv_tables_total)) {
  asv_sort_all[[i]] <- asv_tables_total[[i]][,order(-asv_tables_total[[i]][which(rownames(asv_tables_total[[i]]) == 'total'), ])]
}

## Convert counts to percentages on data matrix & remove row total
asv_pc_all <- list()
for (i in seq_along(asv_sort_all)){
  asv_pc_all[[i]] <- asv_sort_all[[i]][-which(rownames(asv_sort_all[[i]]) == 'total'),] / 
    rowSums(asv_sort_all[[i]][-which(rownames(asv_sort_all[[i]]) == 'total'),]) * 100
}


# create new column "other" which is the sum of all the taxa not in top 25, remove other taxa
other_all <- list()
for (i in seq_along(asv_pc_all)){
  other_all[[i]] <- rowSums(asv_pc_all[[i]][,-c(1:25)])
  }

asv_25_all <- list()
for (i in  seq_along(asv_pc_all)){
  asv_25_all[[i]] <- asv_pc_all[[i]][,1:25]
  asv_25_all[[i]] <- cbind(asv_25_all[[i]], other_all[[i]])
}
  
# Combine with metadata
rel_abun_all <- list()
for (i in seq_along(asv_25_all)) {
  rel_abun_all[[i]] <- cbind(asv_25_all[[i]], metadata_all[[i]])
}

rel_long_all <- lapply(rel_abun_all, melt, variable.name = "Family")



# Format names

unique(rel_long_all[[4]]$Family)

for (i in seq_along(rel_long_all)) {
  rel_long_all[[i]]$Family <- gsub("D_0__Ar.*D_4__", "\\1", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("D_0__Ba.*D_4__", "\\1", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("D_0__Ba.*D_3__", "\\1", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("D_0__Ba.*D_2__", "\\1", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("D_0__", "\\1",  rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("Gammaproteobacteria.__.__", "Gammaproteobacteria unknown", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("Alphaproteobacteria.__.__", "Alphaproteobacteria unknown", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("uncultured.gamma.proteobacterium", "Gammaproteobacteria uncultured", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("Puniceispirillales.Incertae.Sedis", "Puniceispirillales unknown", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("KI89A.clade.__", "KI89A", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("UBA10353.marine.group.__", "UBA10353", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("Entomoplasmatales.Incertae.Sedis", "Entomoplasmatales", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("uncultured.Desulfovibrionales.bacterium", "Desulfovibrionales unknown", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub(".__", "", rel_long_all[[i]]$Family)
  rel_long_all[[i]]$Family <- gsub("V26", "Other", rel_long_all[[i]]$Family)
}


# Split into taxa
rel_a <- rel_long_all[[1]]
rel_c <- rel_long_all[[2]]
rel_o <- rel_long_all[[3]]
rel_s <- rel_long_all[[4]]


# Group and summarise by species

# Ascidian
group_taxon_a <- rel_a %>%  
  mutate(Species = fct_relevel(Species, "L. patella", "D. molle", "P. aurata", "C. foliascens")) %>%
  group_by(Species, Family)
group_taxon_a

taxon_summary_a <- summarise(group_taxon_a, Proportion = mean(value, na.rm = T))
taxon_summary_a <- taxon_summary_a %>% mutate(Family = fct_relevel(Family, "Other", after = Inf))

# Coral
group_taxon_c <- rel_c %>% 
  mutate(Taxonomy = fct_relevel(Species, "P. damicornis", "P. verrucosa", "S. hysterix (RR)", "S. hysterix (PI)",
                                "S. pistillata", "D. heliopora", "E. mammiformis", "P. speciosa", "P. cactus", "A. formosa", 
                                "A. hyacinthus","P. massive", "P. cylindrica (PI)", "P. cylindrica (RR)", "C. foliascens")) %>%
  group_by(Taxonomy, Family)
group_taxon_c
taxon_summary_c <- summarise(group_taxon_c, Proportion = mean(value, na.rm = T))
taxon_summary_c <-taxon_summary_c %>% mutate(Family = fct_relevel(Family, "Other", after = Inf))

# Octocoral 
group_taxon_o <- rel_o %>% 
  mutate(Taxonomy = fct_relevel(Species, "Sarcophyton sp (PI)", "Sarcophyton sp (RR)", "Sinularia sp (PI)", 
                                "Sinularia sp2 (PI)", "Sinularia sp (RR)", "Cladiella sp", "Pinnigorgia sp", "I. hippurus", 
                                "Heteroxenia sp", "Briareum sp", "Briareum sp2", "Clavularia sp", "C. foliascens")) %>%
  group_by(Taxonomy, Family)
group_taxon_o
taxon_summary_o <- summarise(group_taxon_o, Proportion = mean(value, na.rm = T))
taxon_summary_o <- taxon_summary_o %>% mutate(Family = fct_relevel(Family, "Other", after = Inf))


# Sponge
group_taxon_s <- rel_s %>% 
  mutate(Taxonomy = fct_relevel(Species, "I. ramosa (DR)", "I. ramosa (RR)", "Ircinia sp", "Coscinoderma sp",
                                "C. foliascens", "C. singaporense", "Cladiella sp")) %>%
  group_by(Taxonomy, Family)
group_taxon_s
taxon_summary_s <- summarise(group_taxon_s, Proportion = mean(value, na.rm = T))
taxon_summary_s <- taxon_summary_s %>% mutate(Family = fct_relevel(Family, "Other", after = Inf))


# Plot 

# Ascidian
bubblea <- ggplot(taxon_summary_a, aes(x = Species, y = Family, color = Species)) +
  geom_point(aes(size = Proportion)) + 
  scale_size(range = range(taxon_summary_a$Proportion)) +
  scale_size_area(max_size = 10) +
  scale_colour_manual(values = c("#56B4E9", "#56B4E9", "#009E73", "#D55E00")) +
  labs(x = "\nSpecies", y = "Relative Abundance (Family)\n") +
  theme_bw() +
  theme(axis.text=element_text(size=11), axis.text.x=element_text(angle = 0), axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))  + guides(color=FALSE)
bubblea

ggsave(filename = "Acid_Bubble_final.pdf", device = "pdf", width = 27, height = 20, units = "cm", dpi = "print")


# Coral
bubblec <- ggplot(taxon_summary_c, aes(x = Taxonomy, y = Family, color = Taxonomy)) +
  geom_point(aes(size = Proportion)) + 
  scale_size(range = range(taxon_summary_c$Proportion)) +
  scale_size_area(max_size = 10) +
  scale_colour_manual(values = c("#009E73", "#009E73", "#0072B2", "#0072B2", "#0072B2", "#56B4E9",
                                 "#56B4E9", "#F0E442", "#CC79A7", "#D55E00", "#D55E00", "#E69F00", "#E69F00", "#E69F00", "#999999")) +
  labs(x = "\nTaxa", y = "Relative Abundance (Family)\n") +
  theme_bw() +
  theme(axis.text=element_text(size=11), axis.text.x=element_text(angle = 90), axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))  + guides(color=FALSE)
bubblec

ggsave(filename = "coral_bubble_final.pdf", device = "pdf", width = 27, height = 20, units = "cm", dpi = "print")


# Octocoral
bubbleo <- ggplot(taxon_summary_o, aes(x = Taxonomy, y = Family, color = Taxonomy)) +
  geom_point(aes(size = Proportion)) + 
  scale_size(range = range(taxon_summary_o$Proportion)) +
  scale_size_area(max_size = 10) +
  scale_colour_manual(values = c("#0072B2", "#0072B2", "#56B4E9", "#56B4E9", "#56B4E9", "#CC79A7", 
                                 "#009E73", "#009E73", "#F0E442", "#E69F00", "#E69F00", "#D55E00", "#999999")) +
  labs(x = "\nTaxa", y = "Relative Abundance (Family)\n") +
  theme_bw() +
  theme(axis.text=element_text(size=11), axis.text.x=element_text(angle = 90), axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13)) + guides(color=FALSE)
bubbleo

ggsave(filename = "octocoral_bubble_final.pdf", device = "pdf", width = 27, height = 20, units = "cm", dpi = "print")


# Sponge
bubbles <- ggplot(taxon_summary_s, aes(x = Taxonomy, y = Family, color = Taxonomy)) +
  geom_point(aes(size = Proportion)) + 
  scale_size(range = range(taxon_summary_s$Proportion)) +
  scale_size_area(max_size = 10) +
  scale_colour_manual(values = c("#0072B2", "#0072B2", "#CC79A7", "#56B4E9", "#009E73", "#E69F00", "#999999")) +
  labs(x = "\nSpecies", y = "Relative Abundance (Family)\n") +
  theme_bw() +
  theme(axis.text=element_text(size=11), axis.text.x=element_text(angle = 90), axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))  + guides(color=FALSE)
bubbles

ggsave(filename = "Sponge_bubble_final.pdf", device = "pdf", width = 27, height = 20, units = "cm", dpi = "print")


## Calculate summary statistics ----

# get mean and SD
group_list <- list(group_taxon, group_taxonF, group_taxon_a, group_taxon_c, group_taxon_o, group_taxon_s)
sum_list <- list()
for (i in seq_along(group_list)){
  sum_list[[i]] <- summarise(group_list[[i]], average = mean(value, na.rm = T), SD = sd(value, na.rm = T), n = n())
}

sum_df <- lapply(sum_list, as.data.frame)

# get standard error
se <- list()
for (i in seq_along(sum_df)) {
  se[[i]] <- sum_df[[i]]$SD/sqrt(sum_df[[i]]$n) 
  sum_df[[i]] <- cbind(sum_df[[i]], se[[i]])
}


View(sum_df[[6]])



