## Test for phylosymbiosis by; 
## a) comparing topologies between host phylogenetic tree and microbial dendrogram
## b) Correlate host phylogenetic distance with microbial dissimilarity
## c) Plot host phylogeny against microbial dendrogram

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

# Note this script requires the function RFmeasures (Mazel et al. 2018) which can be found at "https://github.com/FloMazel/Phylosymbiosis-Ecological-model"

library(ape)
library(vegan)
library(phangorn)
library(ggtree)
library(ggplot2)
library(cowplot)

## a) comparing topologies between host phylogenetic tree and microbial dendrogram ----

# Import host phylogeny (created using Mr Bayes)
host_tree <- read.nexus("../Phylogenies/Coral/coral_concatenated_msa.nex.con.tre")
host_tree_r <- root(host_tree, "Carteriospongia_foliascens" , resolve.root = T) # set root

# Edit tip labels
tiplabels <- host_tree_r$tip.label
label2 <- c("P. speciosa", "P. cylindrica (RR)","S. hystrix (RR)", "E. mammiformis", 
            "P. cactus","P. damicornis","P. massive","P. cylindrica (PI)",
            "S. hystrix (PI)","S. pistillata","P. verrucosa","C. foliascens",
            "A. formosa", "A. hyacinthus","D. heliopora")

d <- data.frame(label = tiplabels, label2 = label2)
d

host_tree_r$tip.label <- d[,2]
host_tree_r$tip.label

# Import dendrogram (created using QIIME2)
dendrogram <- read.tree("Dendrograms/Filtered/coral_wCar_den_grouped_bc_upgma-min-2-f.tre")

# Edit tip labels
tiplabels2 <- dendrogram$tip.label
label3 <- c("C. foliascens", "P. damicornis", "P. verrucosa", "P. cylindrica (RR)", "P. massive", 
           "A. hyacinthus", "A. formosa", "P. cactus", "P. speciosa", "S. hystrix (RR)", "S. pistillata", 
           "S. hystrix (PI)", "P. cylindrica (PI)", "D. heliopora", "E. mammiformis")

d2 <- data.frame(label = tiplabels2, label2 = label3)
d2


dendrogram$tip.label <- d2[,2]

# Test if tree and dendrogram are congruent
RFmeasures(host_tree_r, dendrogram, nRandom = 9999)


## b) Correlate host phylogenetic distance with microbial dissimilarity ----

# create distance matrix of host phlogenetic distances
host_tree_mat <- cophenetic.phylo(host_tree_r)

# Create distance matrix of ASV table (collapsed by host species in qiime2)
asv_table <- read.table("../Asv_tables/Filtered/Grouped_by_species/coral-wCar-no-doubles-min-2-grouped-f.txt", 
                    sep = '\t', row.names = 1, header = T, strip.white = T)
asv_table_t <- t(asv_table)
table_wis <- wisconsin(asv_table_t) # standardise ASV data

host_names <- row.names(host_tree_mat) # set row names
row.names(table_wis)

host_names <- host_names[c(14, 1, 2, 3, 15, 4, 13, 9, 5, 6, 7, 8, 10, 11, 12)] # change order to correspond with dataframe above (would love to hear of a better way to do this!)
row.names(table_wis) <- host_names 

table_wis <- table_wis[match(row.names(host_tree_mat), rownames(table_wis)),] # Align rows between dataframes

bray_dist <- vegdist(table_wis, method = "bray") # Create distance matrix
bray_dist <- as.matrix(bray_dist)

# Check the two dataframes
head(bray_dist[1:5, 1:5])
head(host_tree_mat[1:5, 1:5])

# Test for correlation
mantel(bray_dist, host_tree_mat, method = "pearson", permutations = 9999)


## c) Plot host phylogeny against microbial dendrogram ----

# Reload phylogeny
host_tree <- read.nexus("../Phylogenies/Coral/coral_concatenated_msa.nex.con.tre")
host_tree_r <- root(host_tree, "Carteriospongia_foliascens" , resolve.root = T) # set root

tiplabels <- host_tree_r$tip.label
label2 <- c("P. speciosa", "P. cylindrica (RR)","S. hystrix (RR)", "E. mammiformis", 
            "P. cactus","P. damicornis","P. massive","P. cylindrica (PI)",
            "S. hystrix (PI)","S. pistillata","P. verrucosa","C. foliascens",
            "A. formosa", "A. hyacinthus","D. heliopora")

d <- data.frame(label = tiplabels, label2 = label2)

# Reload dendrogram
dendrogram <- read.tree("Dendrograms/Filtered/coral_wCar_den_grouped_bc_upgma-min-2-f.tre")

tiplabels2 <- dendrogram$tip.label
label3 <- c("C. foliascens", "P. damicornis", "P. verrucosa", "P. cylindrica (RR)", "P. massive", 
            "A. hyacinthus", "A. formosa", "P. cactus", "P. speciosa", "S. hystrix (RR)", "S. pistillata", 
            "S. hystrix (PI)", "P. cylindrica (PI)", "D. heliopora", "E. mammiformis")

d2 <- data.frame(label = tiplabels2, label2 = label3)

# Plot tree
host_tree_p <- ggtree(host_tree_r, branch.length = "none") %<+% d + geom_tiplab(offset = 0.15, aes(label=label2)) + geom_nodepoint() +
  geom_tippoint(aes(colour = label2), size = 4) + geom_treescale(x=1.5, y=1.4, offset = 0.05) + theme_tree2() + xlim(0,9) + theme_tree()
host_tree_p


# Plot dendrogram
dendrogram_p <- ggtree(dendrogram) %<+% d2 + geom_tiplab(offset = -.21, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  scale_x_reverse(limits=c(0.75,0)) + theme_tree()
dendrogram_p

# Beutify tree
h_nodes <- ggtree(host_tree_r, branch.length = "none") + geom_text(aes(label=node), hjust=-.3) # for nodes
h_nodes

host_tree_p <- ggtree::rotate(host_tree_p, 28)
host_tree_p <- ggtree::rotate(host_tree_p, 17)
host_tree_p <- ggtree::rotate(host_tree_p, 18)
host_tree_p <- ggtree::rotate(host_tree_p, 23)
host_tree_p <- ggtree::rotate(host_tree_p, 24)
host_tree_p <- ggtree::rotate(host_tree_p, 19)
host_tree_p <- ggtree::rotate(host_tree_p, 20)
host_tree_p <- flip(host_tree_p, 7,8)
host_tree_p <- flip(host_tree_p, 8,2)
host_tree_p

# Add clades if needed
w_clades <- host_tree_p +  geom_cladelabel(node = 23, label = "Complex", offset = 2.5, offset.text = 0.15, angle = 90) +
  geom_cladelabel(node = 17, label = "Robust", offset = 2.5, offset.text = 0.15, angle = 90)
w_clades


# Beautify dendrogram
d_nodes <- ggtree(dendrogram) + geom_text(aes(label=node), hjust=-.3) # for node labels
d_nodes

dendrogram_p <- ggtree::rotate(dendrogram_p, 17)
dendrogram_p <- ggtree::rotate(dendrogram_p, 19)
dendrogram_p <- ggtree::rotate(dendrogram_p, 25)
dendrogram_p <- ggtree::rotate(dendrogram_p, 20)
dendrogram_p

# Combine plots and add results
phylo_p <- ggdraw() + 
  draw_plot(w_clades, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(dendrogram_p, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("\nMantel r = 0.35, p = 0.039 \n\n nRF = 0.69, p < 0.001", size = 12, x = 0.5, y = 0.325)
phylo_p







