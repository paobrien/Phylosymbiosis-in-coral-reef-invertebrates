## Calculate and plot species richness and shannon diversity for alpha diversity metrics in microbiome data

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

library(ggplot2)
library(forcats)
library(reshape2)
library(dplyr) 
library(car)
library(multcomp)
library(vegan)

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

cbPalette <- c("#56B4E9", "#009E73", "#F0E442", 
               "#0072B2", "#D55E00", "#CC79A7") # colourblind friendly colour palette

# Import ASV table
asv_table <- read.delim("~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/rarefied_table-f-NoBl_Gor.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table <- t(asv_table)

# Import metadata
metadata <- read.delim("~/Documents/R/Ramaciotti_16S_dataset/Metadata/metadata_NoGor_NoB.tsv", sep = '\t', header = T, row.names = 1, strip.white = T)

## Plot and calculate statistics for species richness (in this case, ASV richness) ----
richness <- apply(asv_table >0, 1, sum)
rich_meta <- cbind(richness, metadata)

# Remove unwanted samples and order factors for plot
rich_meta <-rich_meta[!(rich_meta$Taxonomy=="Mollusc"),]
rich_meta$Taxonomy <- factor(rich_meta$Taxonomy,
                                 levels = c("Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge"),ordered = TRUE)

# Plot as boxplot
box <- ggplot(rich_meta, aes(x=Taxonomy, y=richness, fill=Taxonomy)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  scale_x_discrete(name = "\nTaxa") +
  scale_y_continuous(name = "ASV count\n") +
  scale_fill_manual(values=cbPalette) +
  guides(fill = guide_legend(title = "Taxa")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13)) 
box + geom_text(data=rich.groups.df, aes(x = Taxonomy, y = Ymax+25, label = letters),vjust=0) # for significance

# Calculate statistics
rich_meta$fourth <-(rich_meta$richness ^ (1/4)) # Transform data
residualPlots(lm(rich_meta$fourth~rich_meta$Taxonomy)) # Check assumptions 

asv_fourth.lm <- lm(fourth ~ Taxonomy, rich_meta)
plot(asv_fourth.lm)

summary(asv_fourth.lm) # summarise linear model. Note: taxonomic names are retained when factor levels are not re-ordered
anova(asv_fourth.lm) # for anova

#post hoc test
glht.rich.sum <- summary(glht(asv_fourth.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))
rich.groups <- cld(glht.rich.sum)
rich.groups.df <- fortify(rich.groups)
colnames(rich.groups.df) <- c("Taxonomy", "letters")

ymax <- tapply(rich_meta$richness, rich_meta$Taxonomy, max)
rich.groups.df$Ymax <- ymax # add to plot above



## Plot and calculate statistics for shannon diversity ----
shannon <- diversity(asv_table, index = "shannon")
shan_meta <- cbind(shannon, metadata)

# Remove unwanted samples and reorder factors for plot 
shan_meta <-shan_meta[!(shan_meta$Taxonomy=="Mollusc"),]
shan_meta$Taxonomy <- factor(shan_meta$Taxonomy,
                             levels = c("Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge"),ordered = TRUE)

box2 <- ggplot(shan_meta, aes(x=Taxonomy, y=shannon, fill=Taxonomy)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  scale_x_discrete(name = "\nTaxa") +
  scale_y_continuous(name = "Shannon-Wiener Index (H')\n") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(title = "Taxa")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
box2 + geom_text(data=shan.groups.df, aes(x = Taxonomy, y = Ymax+0.1, label = letters),vjust=0) # for significance letters


# Calculate statistics
residualPlots(lm(shan_meta$shannon~shan_meta$Taxonomy)) # check assumptions 
shannon.lm <- lm(shannon ~ Taxonomy, shan_meta)
plot(shannon.lm)

summary(shannon.lm) 
anova(shannon.lm)

#post hoc test
glht.shan.sum <- summary(glht(shannon.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))
shan.groups <- cld(glht.shan.sum)
shan.groups.df <- fortify(shan.groups)
colnames(shan.groups.df) <- c("Taxonomy", "letters")

ymax <- tapply(shan_meta$shannon, shan_meta$Taxonomy, max)
shan.groups.df$Ymax <- ymax # add to above plot for significance



## Combine plots into one figure ----
rich_meta2 <- rich_meta[,!(colnames(rich_meta) %in% "fourth")]
colnames(rich_meta2)[colnames(rich_meta2)=="richness"] <- "alpha_div"
colnames(shan_meta)[colnames(shan_meta)=="shannon"] <- "alpha_div"
rich_shan <- rbind(rich_meta2, shan_meta)

Alpha <- c(rep("ASV Richness", (nrow(rich_shan)/2)), rep("Shannon-Wiener Index", (nrow(rich_shan)/2))) 
rich_shan <- cbind(rich_shan, Alpha)

rich_shan$Taxonomy <- factor(shan_meta$Taxonomy,
                             levels = c("Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge"),ordered = TRUE)

box3 <- ggplot(rich_shan, aes(x=Taxonomy, y=alpha_div)) +
  geom_boxplot(colour = "#000000", fill = "#56B4E9", outlier.shape = 16, outlier.size = 2) +
  scale_x_discrete(name = "\nHost") +
  scale_y_continuous(name = "") +
  scale_fill_manual(values=cbPalette) +
  guides(fill = guide_legend(title = "Host")) +
  theme_bw() +
  theme(axis.text=element_text(size=11, colour = "black"),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13)) 
box3 + facet_wrap(~Alpha, nrow = 2, ncol = 1, scales = "free_y") 

  
ggsave("Alpha_diversity_final.pdf", width = 10, height = 8, device = "pdf", dpi = "print")

