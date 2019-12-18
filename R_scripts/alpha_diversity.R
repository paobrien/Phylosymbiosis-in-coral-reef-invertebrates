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


# Import ASV table
asv_table <- read.delim("../Asv_tables/Filtered/rarefied_table-f-NoBl_Gor.txt", sep = '\t', header = T, row.names = 1, strip.white = T)
asv_table <- t(asv_table)

# Import metadata
metadata <- read.delim("../Metadata/metadata_NoGor_NoB.tsv", sep = '\t', header = T, row.names = 1, strip.white = T)

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
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(title = "Taxa")) +
  theme_bw() +
  theme(axis.text=element_text(size=11),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
box

# Calculate statistics
rich_meta$fourth <-(rich_meta$richness ^ (1/4)) # Transform data
residualPlots(lm(rich_meta$fourth~rich_meta$Taxonomy)) # Check assumptions 
asv_fourth.lm <- lm(fourth ~ Taxonomy, rich_meta)
plot(asv_fourth.lm)

summary(asv_fourth.lm) # summarise linear model. Note: taxonomic names a retained when factor levels are not re-ordered
anova(asv_fourth.lm) # for anova

#post hoc test
summary(glht(asv_fourth.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))


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
box2


# Calculate statistics
residualPlots(lm(shan_meta$shannon~shan_meta$Taxonomy)) # check assumptions 
shannon.lm <- lm(shannon ~ Taxonomy, shan_meta)
plot(shannon.lm)

summary(shannon.lm) 
anova(shannon.lm)

#post hoc test
summary(glht(shannon.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))

## Combine plots into one figure ----
rich_meta2 <- rich_meta[,!(colnames(rich_meta) %in% "fourth")]
colnames(rich_meta2)[colnames(rich_meta2)=="richness"] <- "alpha_div"
colnames(shan_meta)[colnames(shan_meta)=="shannon"] <- "alpha_div"
rich_shan <- rbind(rich_meta2, shan_meta)

Alpha <- c(rep("ASV Richness", (nrow(rich_shan)/2)), rep("Shannon-Wiener Index", (nrow(rich_shan)/2))) 
rich_shan <- cbind(rich_shan, Alpha)

rich_shan$Taxonomy <- factor(shan_meta$Taxonomy,
                             levels = c("Seawater", "Coral", "Soft Coral", "Gorgonian", "Ascidian", "Sponge"),ordered = TRUE)

box3 <- ggplot(rich_shan, aes(x=Taxonomy, y=alpha_div, fill=Taxonomy)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
  scale_x_discrete(name = "\nHost") +
  scale_y_continuous(name = "") +
  scale_fill_brewer(palette = "Dark2") +
  guides(fill = guide_legend(title = "Host")) +
  theme_bw() +
  theme(axis.text=element_text(size=11, colour = "black"),axis.title=element_text(size = 13),
        legend.text=element_text(size=11), legend.title=element_text(size=13))
box3 + facet_wrap(~Alpha, nrow = 2, ncol = 1, scales = "free_y")














