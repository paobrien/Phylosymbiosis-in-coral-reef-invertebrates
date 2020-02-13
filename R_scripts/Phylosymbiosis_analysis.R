## Test for phylosymbiosis by; 
## a) comparing topologies between host phylogenetic tree and microbial dendrogram
## b) Correlate host phylogenetic distance with microbial dissimilarity

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# Note this script requires the function RFmeasures (Mazel et al. 2018) which can be found at "https://github.com/FloMazel/Phylosymbiosis-Ecological-model"

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(ggplot2)
library(GUniFrac)

## a) comparing topologies between host phylogenetic tree and microbial dendrogram ----


## Import host phylogeny (created using Mr Bayes)

host_filepath <- "~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/"
host_files <- list.files(host_filepath)

host_trees <- list()
for (f in seq_along(host_files)) {
  host_trees[[f]] <- read.nexus(paste0(host_filepath, file = host_files[[f]]))
}

host_trees[[1]] <- host_trees[[1]]$con_50_majrule # fix import error for ascidians

# Root host phylogenies
for (i in seq(1, 3)){
  host_trees[[i]] <- root(host_trees[[i]], "Carteriospongia_foliascens" , resolve.root = T) # set root
}
host_trees[[4]] <- root(host_trees[[4]], "Cladiella_sp_1399" , resolve.root = T) # set root

# Edit tip labels
tiplabels <- list()
for (i in seq_along(host_trees)){
  tiplabels[[i]] <- host_trees[[i]]$tip.label
}

 
label1 <- c("D. molle",  "L. patella", "P. aurata", "C. foliascens")
label2 <- c("P. speciosa", "P. cylindrica (RR)","S. hystrix (RR)", "E. mammiformis", 
            "P. cactus","P. damicornis","P. massive","P. cylindrica (PI)",
            "S. hystrix (PI)","S. pistillata","P. verrucosa","C. foliascens",
            "A. formosa", "A. hyacinthus","D. heliopora")
label3 <- c("Briareum sp2", "Briareum sp", "Clavularia sp",
            "Sinularia sp (RR)", "Sinularia sp2 (PI)", "Sinularia sp (PI)",
            "Sarcophyton sp (RR)", "Sarcophyton sp (PI)", "Cladiella sp", 
            "Heteroxenia sp","Pinnigorgia sp", "I. hippuris", 
            "C. foliascens")
label4 <- c("Coscinoderma sp", "I. ramosa (RR)", "C. foliascens", "Ircinia sp",
            "I. ramosa (DR)", "C. singaporense", "Cladiella sp")

lbs <- list(label1, label2, label3, label4) 
newlbs <- list()
for (i in seq_along(tiplabels)){
  newlbs[[i]] <- data.frame(tiplabels[[i]], lbs[[i]]) 
}

newlbs 

for (i in seq_along(newlbs)){
  colnames(newlbs[[i]]) <- c("label", "label2") 
}

for (i in seq_along(host_trees)){
  host_trees[[i]]$tip.label <- as.character(newlbs[[i]][,2])
}



# Import dendrograms (created using QIIME2 - Bray Curtis, weighted and unweighted UniFrac)
den_filepath <- "~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/"

a_den_files <- list.files(den_filepath, pattern = "^a")
c_den_files <- list.files(den_filepath, pattern = "^c")
o_den_files <- list.files(den_filepath, pattern = "^o")
s_den_files <- list.files(den_filepath, pattern = "^s")

a_den <- list()
c_den <- list()
o_den <- list()
s_den <- list()

for (i in seq_along(a_den_files)) {
  a_den[[i]] <- read.tree(paste0(den_filepath, file = a_den_files[i]))
}

for (i in seq_along(c_den_files)) {
  c_den[[i]] <- read.tree(paste0(den_filepath, file = c_den_files[i]))
}

for (i in seq_along(o_den_files)) {
  o_den[[i]] <- read.tree(paste0(den_filepath, file = o_den_files[i]))
}

for (i in seq_along(s_den_files)) {
  s_den[[i]] <- read.tree(paste0(den_filepath, file = s_den_files[i]))
}

dendrograms <- list(a_den, c_den, o_den, s_den)
names(dendrograms) <- c("a_den", "c_den", "o_den", "s_den")

# Edit tip labels
for (j in seq_along(dendrograms)){
  temp_den<-dendrograms[[j]]
  for (i in seq_along(temp_den)) {
    temp_den[[i]]$tip.label <- gsub("_", " ", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("OIRS", "(PI)", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("CS", "(RR)", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("I.r", "I. r", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("Pachyseris sp", "P. speciosa", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("damnicornis", "damicornis", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("S.hysterix", "S. hystrix", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("P.cylindrica", "P. cylindrica", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("Sarcophytonsp", "Sarcophyton sp", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("Sinulariasp", "Sinularia sp", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("Lobophytum sp", "Sinularia sp2 (PI)", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("I. hippurus", "I. hippuris", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("Xenia sp", "Heteroxenia sp", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("I. ramosa Dav", "I. ramosa (DR)", temp_den[[i]]$tip.label)
    temp_den[[i]]$tip.label <- gsub("I. blue", "Ircinia sp", temp_den[[i]]$tip.label)
   }
  dendrograms[[j]] <- temp_den
}

# check if label match
a_den <- dendrograms[[1]]
c_den <- dendrograms[[2]]
o_den <- dendrograms[[3]]
s_den <- dendrograms[[4]]

den_labs <- list(a_den[[1]], c_den[[1]], o_den[[1]], s_den[[1]])

for (i in seq_along(host_trees)){
  print(is.element(host_trees[[i]]$tip.label, den_labs[[i]]$tip.label))
}


# Test if tree and dendrogram are congruent
for (i in seq_along(host_trees)) {
  print(is.binary(host_trees[[i]])) # ensure host trees are binary
}

RF_a <- list()
for (i in seq_along(a_den)) {
  RF_a[[i]] <- RFmeasures(host_trees[[1]], a_den[[i]], nRandom = 9999)
}

RF_c <- list()
for (i in seq_along(c_den)) {
  RF_c[[i]] <- RFmeasures(host_trees[[2]], c_den[[i]], nRandom = 9999)
}

RF_o <- list()
for (i in seq_along(o_den)) {
  RF_o[[i]] <- RFmeasures(host_trees[[3]], o_den[[i]], nRandom = 9999)
}

RF_s <- list()
for (i in seq_along(s_den)) {
  RF_s[[i]] <- RFmeasures(host_trees[[4]], s_den[[i]], nRandom = 9999)
}



## b) Correlate host phylogenetic distance with microbial dissimilarity ----

# create distance matrix of host phylogenetic distances
host_tree_mat <- list()
for (i in seq_along(host_trees)) {
  host_tree_mat[[i]] <- cophenetic.phylo(host_trees[[i]]) 
}

# Create BC and UniFrac distance matrices of ASV table (collapsed by host species in qiime2)

# Import and format table
asv_filepath <- "~/Documents/R/Ramaciotti_16S_dataset/Asv_tables/Filtered/Grouped_by_species/"
asv_files <- list.files(asv_filepath)

asv <- list()
for (f in asv_files){
  asv[[f]] <-read.table(paste0(asv_filepath, f), sep = '\t', row.names = 1, header = T, strip.white = T)
}

asv_t <- lapply(asv, t)

# Edit row names
for (i in seq_along(asv_t)) {
  row.names(asv_t[[i]]) <- gsub("\\.\\.", "\\. ", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("_CS", " (RR)", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("_OIRS", " (PI)", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("Pachyseris.sp", "P. speciosa", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("hysterix", "hystrix", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("damnicornis", "damicornis", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("m\\.", "m ", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("n\\.", "n ", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("a\\.", "a ", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("Xenia sp", "Heteroxenia sp", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("I. hippurus", "I. hippuris", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("Lobophytum sp", "Sinularia sp2 (PI)", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("I. blue", "Ircinia sp", row.names(asv_t[[i]]))  
  row.names(asv_t[[i]]) <- gsub("I. ramosa_Dav", "I. ramosa (DR)", row.names(asv_t[[i]]))  
  }

# Check matching row names
for (i in seq_along(asv_t)){
  print(is.element(row.names(asv_t[[i]]), row.names(host_tree_mat[[i]])))
}


# align rows between dataframes
for (i in seq_along(asv_t)) {
  asv_t[[i]] <- asv_t[[i]][match(row.names(host_tree_mat[[i]]), rownames(asv_t[[i]])),] 
}

# standardise data
asv_wis <- lapply(asv_t, wisconsin) # more approriate for Bray-Curtis, but check both

asv_r <- list()
asv_r[[1]] <- rrarefy(asv_t[[1]], min(apply(asv_t[[1]], 1, sum))) # Rarefy more appropriate for UniFrac
asv_r[[2]] <- rrarefy(asv_t[[2]], min(apply(asv_t[[2]], 1, sum)))
asv_r[[3]] <- rrarefy(asv_t[[3]], min(apply(asv_t[[3]], 1, sum)))
asv_r[[4]] <- rrarefy(asv_t[[4]], min(apply(asv_t[[4]], 1, sum)))

# Test for Bray-Curtis correlation
bc_dist <- lapply(asv_wis, vegdist, method = "bray")
bc_mat  <- lapply(bc_dist, as.matrix)

mantel_stats <- list()
for (i in seq_along(bc_mat)) {
  mantel_stats[[i]] <- mantel(bc_mat[[i]], host_tree_mat[[i]], method = "pearson", permutations = 9999)
}

# Check stats
mantel_stats

# For UniFrac
asv_tree <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Bacteria_rep_seqs/insertion_tree_trimmed.nwk") # import bacteria phylogeny

UF <- list()
for (i in seq_along(asv_r)) {
  UF[[i]] <- GUniFrac(asv_r[[i]], asv_tree, alpha = c(0, 0.5,1))
}

wUF <- list()  # weighted
uwUF <- list() # unweighted
for (i in seq_along(UF)) {
  wUF[[i]] <- UF[[i]]$unifracs[, , "d_1"]
  uwUF[[i]] <- UF[[i]]$unifracs[, , "d_UW"]
}

# Test for correlation
mantel_wUF <- list()
mantel_uwUF <- list()
for (i in 1:4) {
  mantel_wUF[[i]] <- mantel(wUF[[i]], host_tree_mat[[i]], method = "pearson", permutations = 9999)
  mantel_uwUF[[i]] <- mantel(uwUF[[i]], host_tree_mat[[i]], method = "pearson", permutations = 9999)
}

# Check stats
mantel_wUF
mantel_uwUF

# Plot using Plot Phylosymbiosis script



