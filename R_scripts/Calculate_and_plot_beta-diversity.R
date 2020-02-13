## Calculate and Bray-Curtis beta-diversity statistics and plot using NMDS

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

library(vegan)
library(ggplot2)
library(gg3D)

cbPalette <- c("Seawater" = "#56B4E9", "Coral" = "#009E73", "Soft Coral" = "#F0E442", 
               "Gorgonian" = "#0072B2", "Ascidian" = "#D55E00", "Sponge" = "#CC79A7") # colourblind friendly colour palette


## Plot ----

# Import data
asv_table <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/table-f-NoBl_Gor.txt", sep = '\t', row.names = 1, header = T, strip.white = T)
metadata <- read.table("~/Documents/R/Ramaciotti_16S_dataset/Metadata/metadata_NoBlankOrGor.tsv", sep = '\t', row.names = 1, header = T, strip.white = T)

asv_table_wide <- t(asv_table) # tanspose to wide format

# remove unwanted samples
asv_table_wide <- asv_table_wide[-grep("Clam", rownames(asv_table_wide)), ] 
metadata <- metadata[!(metadata$Taxonomy == "Mollusc"),]


# standardise data and create distance matrix. Note: wisconsin standardisation with bray-curtis optimal for this data set
asv_wis <- wisconsin(asv_table_wide)

# calculate NMDS scores
# k=3 for third dimension; decreases stress score and better represents the pattern in the dissimilarity matrix
bray.nmds <- metaMDS(asv_wis, distance = "bray", trymax = 200, autotransform = FALSE, k = 3, plot = TRUE) 
bray.nmds

goodness(bray.nmds) # large values = poor fit
stressplot(bray.nmds)


# Extract scores from metaMDS and return as data matrix
nmds_scores <- as.data.frame(scores(bray.nmds))

#add metadata
nmds_scores_meta <- cbind(nmds_scores, metadata)
nmds_scores_meta$Sample <- rownames(nmds_scores_meta)

# reorder levels
nmds_scores_meta$Taxonomy <- factor(nmds_scores_meta$Taxonomy,
                                    levels = c("Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge"),ordered = TRUE)

theta=0
phi=90
plot2.2 <- ggplot(nmds_scores_meta, aes(x=NMDS1, y=NMDS2, z=NMDS3)) +
  scale_shape_manual(values = c(21,24,22)) +
  scale_fill_manual(values = cbPalette, name="Host") +
  guides(fill = guide_legend(override.aes=list(shape=21), title = "Host")) +
  labs(color = "Host") +
  axes_3D(theta = theta, phi=phi) +
  stat_3D(aes(fill = Taxonomy, shape=Zone), size=3.5, alpha = 0.99) +
  #labs_3D(labs=c("NMDS1", "NMDS2", "NMDS3"),  hjust=c(-2.2,1,1), vjust=c(0.5, 0.5, 0.5), angle=c(0, 90, 0)) +
  theme_void() 
plot2.2 + annotate("text", x= c(-0.17, -0.24, -0.32, -0.17), y = c(-0.23, 0, -0.28, -0.28), 
                   label = c("NMDS1", "NMDS2", "NMDS3", "Stress = 0.17"))

ggsave(filename = "NMDS_3D_Final_corrected2_autoT.pdf", width = 20, height = 15, units = "cm", device = "pdf", dpi = "print")


## Statistical analysis on dist matrix ----
#PERMANOVA
data.perm <- cbind(asv_wis, metadata)
bray_dist <- vegdist(asv_wis, method = "bray")

#By taxonomy
perm_taxa <- adonis(bray_dist~Taxonomy, data.perm, permutations = 999)
perm_taxa

# By species
perm_sp <- adonis(bray_dist~Species, data.perm, permutations = 999)
perm_sp

# By site
perm_site <- adonis(bray_dist~Site, data.perm, permutations = 999)
perm_site

# By reef zone
perm_zone <- adonis(bray_dist~Zone, data.perm, permutations = 999)
perm_zone

# By collection date
perm_date <- adonis(bray_dist~Date, data.perm, permutations = 999)
perm_date




