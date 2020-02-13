## Compare intraspecic with interspecific variation in the microbiome ##
# Paul A O'Brien
# paul.obrien@my.jcu.edu.au

library(vegan)
library(dplyr)
library(ggplot2)
library(multcomp)

## Set up data ----

# Import data
setwd("/Users/Paul/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered")
asv_files <- list.files(pattern = "*table-f.txt")
asv_data <- list()
for (i in seq_along(asv_files)) {
  asv_data[[i]] <- read.table(file = asv_files[i], sep = '\t', row.names = 1, header = T, strip.white = T)
}

names(asv_data) <- gsub("-table-f.txt", "", asv_files)

setwd("/Users/Paul/Documents/R/Ramaciotti_16S_dataset/Metadata")
meta_files <- list.files(pattern = "_metadata.tsv")
meta_data <- list()
for (i in seq_along(meta_files)) {
  meta_data[[i]] <- read.table(file = meta_files[i], sep = '\t', header = T, strip.white = T)
}

names(meta_data) <- gsub(".tsv", "", meta_files)

# Transpose ASV data
asv_t <- lapply(asv_data, t)

# Standardise data and create long format distance matrix 
asv_wis <- lapply(asv_t, wisconsin)
bray.dist <- lapply(asv_wis, vegdist, method = "bray")
dist.mat <- lapply(bray.dist, as.matrix)  

dist.melt <- list()
for (i in seq_along(dist.mat)) {
  dist.melt[[i]] <- data.frame(X=colnames(dist.mat[[i]])[col(dist.mat[[i]])], Y=rownames(dist.mat[[i]])[row(dist.mat[[i]])], dist=c(dist.mat[[i]]))
}  

names(dist.melt) <- names(asv_data)

# Create new column with intraspecific vs. interspecific
# Ascidian
dist.melt[[1]]$intra.sp <- ifelse(dist.melt[[1]]$X %in% c("Acid1", "Acid2", "Acid3") & dist.melt[[1]]$Y %in% c("Acid1", "Acid2", "Acid3") |
                           dist.melt[[1]]$X %in% c("Lis1", "Lis2", "Lis3", "Lis4", "Lis5") & dist.melt[[1]]$Y %in% c("Lis1", "Lis2", "Lis3", "Lis4", "Lis5") |
                           dist.melt[[1]]$X %in% c("Pol1", "Pol2", "Pol3", "Pol4", "Pol5") & dist.melt[[1]]$Y %in% c("Pol1", "Pol2", "Pol3", "Pol4", "Pol5"), "intraspecific", "interspecific")

# Coral
dist.melt[[2]]$intra.sp <- ifelse(dist.melt[[2]]$X %in% c("Ahy2", "Ahy3", "Ahy4") & dist.melt[[2]]$Y %in% c("Ahy2", "Ahy3", "Ahy4") |
                     dist.melt[[2]]$X %in% c("CS18", "CS44", "CS51") & dist.melt[[2]]$Y %in% c("CS18", "CS44", "CS51") |
                     dist.melt[[2]]$X %in% c("CS59", "CS60", "CS61", "CS62", "CS63", "Por1", "Por2", "Por3") & 
                     dist.melt[[2]]$Y %in% c("CS59", "CS60", "CS61", "CS62", "CS63", "Por1", "Por2", "Por3") |
                     dist.melt[[2]]$X %in% c("CS64", "CS65", "CS66", "CS67", "CS68", "Hist1", "Hist2", "Hist3", "Shy1", "Shy2", "Shy3") & 
                     dist.melt[[2]]$Y %in% c("CS64", "CS65", "CS66", "CS67", "CS68", "Hist1", "Hist2", "Hist3", "Shy1", "Shy2", "Shy3") |
                     dist.melt[[2]]$X %in% c("Dip1", "Dip2", "Dip3") & dist.melt[[2]]$Y %in% c("Dip1", "Dip2", "Dip3") |
                     dist.melt[[2]]$X %in% c("Eck1", "Eck2", "Eck3", "Eck4", "Eck5") & dist.melt[[2]]$Y %in% c("Eck1", "Eck2", "Eck3", "Eck4", "Eck5") |
                     dist.melt[[2]]$X %in% c("For1", "For2", "For3", "For4", "For5") & dist.melt[[2]]$Y %in% c("For1", "For2", "For3", "For4", "For5") |
                     dist.melt[[2]]$X %in% c("Pav1", "Pav2", "Pav3", "Pav4", "Pav5") & dist.melt[[2]]$Y %in% c("Pav1", "Pav2", "Pav3", "Pav4", "Pav5") |
                     dist.melt[[2]]$X %in% c("Pdam1", "Pdam2", "Pdam3") & dist.melt[[2]]$Y %in% c("Pdam1", "Pdam2", "Pdam3") |
                     dist.melt[[2]]$X %in% c("Pm1", "Pm2", "Pm3", "Pm4", "Pm5", "Pmas1", "Pmas2", "Pmas3") & 
                     dist.melt[[2]]$Y %in% c("Pm1", "Pm2", "Pm3", "Pm4", "Pm5", "Pmas1", "Pmas2", "Pmas3") |
                     dist.melt[[2]]$X %in% c("Sty1", "Sty2", "Sty3", "Sty4", "Sty5") & dist.melt[[2]]$Y %in% c("Sty1", "Sty2", "Sty3", "Sty4", "Sty5") |
                     dist.melt[[2]]$X %in% c("Ver1", "Ver2", "Ver3", "Ver4", "Ver5") & dist.melt[[2]]$Y %in% c("Ver1", "Ver2", "Ver3", "Ver4", "Ver5"), "intraspecific", "interspecific")

# Octocoral
dist.melt[[3]]$intra.sp <- ifelse(dist.melt[[3]]$X %in% c("Br1", "Br2", "Br3", "Br4", "Br5") & dist.melt[[3]]$Y %in% c("Br1", "Br2", "Br3", "Br4", "Br5") |
                                    dist.melt[[3]]$X %in% c("Br2.1", "Br2.2", "Br2.3", "Br2.4", "Br2.5") & dist.melt[[3]]$Y %in% c("Br2.1", "Br2.2", "Br2.3", "Br2.4", "Br2.5") |
                                    dist.melt[[3]]$X %in% c("Cav1", "Cav2", "Cav3", "Cav4", "Cav5") & dist.melt[[3]]$Y %in% c("Cav1", "Cav2", "Cav3", "Cav4", "Cav5") |
                                    dist.melt[[3]]$X %in% c("CS74", "CS75", "CS76", "CS77", "CS78", "Sin1", "Sin2", "Sin3", "Sin4", "Sin5") & 
                                    dist.melt[[3]]$Y %in% c("CS74", "CS75", "CS76", "CS77", "CS78", "Sin1", "Sin2", "Sin3", "Sin4", "Sin5") |
                                    dist.melt[[3]]$X %in% c("CS79", "CS80", "CS81", "CS82", "CS83", "Sar1", "Sar2", "Sar3") & 
                                    dist.melt[[3]]$Y %in% c("CS79", "CS80", "CS81", "CS82", "CS83", "Sar1", "Sar2", "Sar3") |
                                    dist.melt[[3]]$X %in% c("Hel1", "Hel2", "Hel3") & dist.melt[[3]]$Y %in% c("Hel1", "Hel2", "Hel3") |
                                    dist.melt[[3]]$X %in% c("Isis1", "Isis2", "Isis3", "Isis4", "Isis5") & dist.melt[[3]]$Y %in% c("Isis1", "Isis2", "Isis3", "Isis4", "Isis5") |
                                    dist.melt[[3]]$X %in% c("Lobo1", "Lobo2", "Lobo3", "Lobo4", "Lobo5") & dist.melt[[3]]$Y %in% c("Lobo1", "Lobo2", "Lobo3", "Lobo4", "Lobo5") |
                                    dist.melt[[3]]$X %in% c("Pin1", "Pin2", "Pin3", "Pin4", "Pin5") & dist.melt[[3]]$Y %in% c("Pin1", "Pin2", "Pin3", "Pin4", "Pin5") |
                                    dist.melt[[3]]$X %in% c("Usc1", "Usc2", "Usc3", "Usc4", "Usc5") & dist.melt[[3]]$Y %in% c("Usc1", "Usc2", "Usc3", "Usc4", "Usc5"), "intraspecific", "interspecific")

# Sponge 
dist.melt[[4]]$intra.sp <- ifelse(dist.melt[[4]]$X %in% c("Car1", "Car2", "Car3", "Car4", "Car5") & dist.melt[[4]]$Y %in% c("Car1", "Car2", "Car3", "Car4", "Car5") |
                      dist.melt[[4]]$X %in% c("CS54", "CS55", "CS56", "CS57", "CS58") & dist.melt[[4]]$Y %in% c("CS54", "CS55", "CS56", "CS57", "CS58") |
                      dist.melt[[4]]$X %in% c("CS69", "CS70", "CS71", "CS72", "CS73", "Irc1", "Icr2", "Irc3", "Irc4", "Irc5") & 
                      dist.melt[[4]]$Y %in% c("CS69", "CS70", "CS71", "CS72", "CS73", "Irc1", "Icr2", "Irc3", "Irc4", "Irc5") |
                      dist.melt[[4]]$X %in% c("IrB1", "IrB2", "IrB3", "IrB4", "IrB5") & dist.melt[[4]]$Y %in% c("IrB1", "IrB2", "IrB3", "IrB4", "IrB5") |
                      dist.melt[[4]]$X %in% c("Spg1", "Spg2", "Spg3") & dist.melt[[4]]$Y %in% c("Spg1", "Spg2", "Spg3"), "intraspecific", "interspecific")




# Remove 0's and arcsine tranformation to normalise data
for (i in seq_along(dist.melt)) {
  dist.melt[[i]] <- dist.melt[[i]][!(dist.melt[[i]]$dist ==0),]
}

for (i in seq_along(dist.melt)) {
  dist.melt[[i]]$arcsin.t <- asin(sqrt(dist.melt[[i]]$dist))*180/pi
}


## Statistical analysis ----

# T-test to see if there is a significant difference between interspecific and intraspecific variation
t_test <- list()
for (i in seq_along(dist.melt)){
  t_test[[i]] <- t.test(arcsin.t~intra.sp, dist.melt[[i]], var.equal = FALSE)
}


# ANOVA to test if intraspecific or interspecific variation is different among the groups
dist.df <- bind_rows(dist.melt, .id = "Taxonomy") # bind lists to a dataframe
intraspecific <- dist.df[dist.df$intra.sp == "intraspecific",]
interspecific <- dist.df[dist.df$intra.sp == "interspecific",]

intraspecific$Taxonomy <- as.factor(intraspecific$Taxonomy)
interspecific$Taxonomy <- as.factor(interspecific$Taxonomy)

intra.lm <- lm(arcsin.t ~ Taxonomy, intraspecific) # Set up linear model for anova
par(mfrow = c(2,2))
plot(intra.lm) # check assumptions

summary(intra.lm)
anova(intra.lm)

inter.lm <- lm(arcsin.t ~ Taxonomy, interspecific)
par(mfrow = c(2,2))
plot(inter.lm)

summary(inter.lm)
anova(inter.lm)

#post hoc test
summary(glht(intra.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))
summary(glht(inter.lm, linfct = mcp(Taxonomy = "Tukey")), test = adjusted("bonferroni"))


## Plot data ----

p <- ggplot(dist.df, aes(x=Taxonomy, y=arcsin.t, fill=intra.sp)) + 
  geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) + 
  ylim(0, NA) +
  ylab("Bray-curtis dissimilarity (arcsine transformed)\n") + xlab("\nHost") +
  scale_fill_brewer(palette = "Set1") +
  labs(fill = "") + 
  theme_bw() +
  theme(axis.text=element_text(size=15, colour = "black"),axis.title=element_text(size = 18),
        legend.text=element_text(size=15)) 
p + annotate("text", x = 1, y = 2, label = "t=19.1, p < 0.001") + annotate("text", x = 2, y = 2, label = "t=13.5, p < 0.001") +
  annotate("text", x = 3, y = 2, label = "t=18.8, p < 0.001") + annotate("text", x = 4, y = 2, label = "t=34.8, p < 0.001") 



