## Plot dendrogram against host phylogeny 

# Paul A. O'Brien
# E. paul.obrien@my.jcu.edu.au

setwd("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/")

library(ggtree)
library(ggplot2)
library(cowplot)
library(ape)


## Ascidian ----

## Host tree
a_tree <- read.nexus("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/ascidian_concatenate_msa.nex.con.tre")
a_tree <- a_tree$con_50_majrule
a_tree_r <- root(a_tree, "Carteriospongia_foliascens", resolve.root = T)

# Edit tips
tiplabels <- a_tree_r$tip.label  
a <- data.frame(label = tiplabels, label2 = c("D. molle (PI)",  "L. patella (DR)", "P. aurata (BR)", "C. foliascens (DR)"))

# Set colours
a_cols <- c("D. molle (PI)" = "#56B4E9", "L. patella (DR)" = "#56B4E9", "P. aurata (BR)" = "#009E73", "C. foliascens (DR)" = "#D55E00")

host_tree_a <- ggtree(a_tree_r, branch.length = "none") %<+%  a + 
  geom_tiplab(aes(label=label2), offset = 0.1) + 
  geom_nodepoint() + geom_tippoint(aes(color = label2), size=4) + 
  geom_nodelab(aes(label = label), nudge_x = 0.15) +
  geom_treescale(x=1.5, y=1.3, offset = 0.05) + theme_tree2() + xlim(0, 4) + theme_tree() +
  scale_color_manual(values = a_cols) +
  #geom_text(data=a_dat, aes(label=label), nudge_x = 0.15) 
host_tree_a

# To format node labels - replace geom_nodelab with geom_text above
a_dat <- host_tree_a$data # ggtree object
a_dat <- a_dat[!a_dat$isTip,]
a_dat$label <- as.numeric(a_dat$label)
a_dat$label <- round(a_dat$label, 2)


## Bray Curtis
a_bc <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/ascidian_wCar_den_grouped_bc_upgma-min-2-f.tre")

tiplabels <- a_bc$tip.label  
a_b <- data.frame(label = tiplabels, label2 = c("C. foliascens (DR)", "P. aurata (BR)", "D. molle (PI)",  "L. patella (DR)"))

den_a_bc <- ggtree(a_bc) %<+% a_b + geom_tiplab(offset = -.21, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(nudge_x = -0.02) + 
  scale_x_reverse(limits=c(0.75,0)) + theme_tree() +
  scale_color_manual(values = a_cols)
den_a_bc

# Combine dendrogram and host phylogeny 
asc1 <- ggdraw() + 
  draw_plot(host_tree_a, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_a_bc, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = -0.03, p = 0.63 \n\n nRF = 0, p < 0.001", size = 12, x = 0.5, y = 0.325)
asc1  

ggsave(filename = "Ascidian_bc_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

# Weighted UF
a_wuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/ascidian_den_FI-tree_trim_wuf_upgma.tre")

tiplabels <- a_wuf$tip.label
a_w <- data.frame(label = tiplabels, label2 = c("C. foliascens (DR)", "P. aurata (BR)", "D. molle (PI)",  "L. patella (DR)"))

den_a_wuf <- ggtree(a_wuf) %<+% a_w + geom_tiplab(offset = -.125, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.425,0)) + theme_tree() +
  scale_color_manual(values = a_cols)
den_a_wuf

# Combine dendrogram and host phylogeny 
asc2 <- ggdraw() + 
  draw_plot(host_tree_a, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_a_wuf, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.46, p = 0.17 \n\n nRF = 0, p < 0.001", size = 12, x = 0.5, y = 0.325)
asc2  

ggsave(filename = "Ascidian_wUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")


# Unweighted UF
a_uw <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/ascidian_den_FI-tree_trim_uwuf_upgma.tre")

tiplabels <- a_uw$tip.label
a_u <- data.frame(label = tiplabels, label2 = c("C. foliascens (DR)", "D. molle (PI)", "L. patella (DR)", "P. aurata (BR)"))

anodes <- ggtree(a_uw) + geom_text(aes(label=node), hjust=-.3) # node numbers

den_a_uw <- ggtree(a_uw) %<+% a_u + geom_tiplab(offset = -.175, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab() +
  scale_x_reverse(limits=c(0.62,0)) + theme_tree() +
  scale_color_manual(values = a_cols) 
  #geom_text(data=den_a_dat, aes(label=label), nudge_x = -0.03) 
den_a_uw 

# Format node labels and change as above
den_a_dat <- den_a_uw$data
den_a_dat <- den_a_dat[!den_a_dat$isTip,]
den_a_dat$label <- as.numeric(den_a_dat$label)
den_a_dat$label <- round(den_a_dat$label, 3)

# Beautify tree
den_a_uw <- ggtree::rotate(den_a_uw, 7)

# Combine dendrogram and host phylogeny 
asc3 <- ggdraw() + 
  draw_plot(host_tree_a, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_a_uw, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.17, p = 0.46 \n\n nRF = 0.5, p = 0.34", size = 12, x = 0.5, y = 0.325)
asc3  

ggsave(filename = "Ascidian_uwUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")



## Coral ----

## HOST TREE
c_tree <- read.nexus("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/coral_concatenated_msa.nex.con.tre")
c_tree_r <- root(c_tree, "Carteriospongia_foliascens", resolve.root = T)

# Edit tips
tiplabels <- c_tree_r$tip.label  
tiplabels2 <- c("P. speciosa (OR)", "P. cylindrica (RR)","S. hystrix (RR)", "E. mammiformis (PI)", 
                "P. cactus (PI)","P. damicornis (PI)","P. massive (PI)","P. cylindrica (PI)",
                "S. hystrix (PI)","S. pistillata (BR)","P. verrucosa (RB)","C. foliascens (DR)",
                "A. formosa (PI)", "A. hyacinthus (PI)","D. heliopora (PI)")
c <- data.frame(label = tiplabels, label2 = tiplabels2)

# Set colours
c_cols <- c("P. speciosa (OR)" = "#F0E442", "P. cylindrica (RR)" = "#E69F00", "S. hystrix (RR)" = "#0072B2", "E. mammiformis (PI)" = "#56B4E9", 
            "P. cactus (PI)" = "#CC79A7", "P. damicornis (PI)" = "#009E73", "P. massive (PI)"  = "#E69F00", "P. cylindrica (PI)"  = "#E69F00",
            "S. hystrix (PI)"  = "#0072B2", "S. pistillata (BR)" = "#0072B2", "P. verrucosa (RB)" = "#009E73", "C. foliascens (DR)" = "#999999",
            "A. formosa (PI)" = "#D55E00", "A. hyacinthus (PI)" = "#D55E00","D. heliopora (PI)" = "#56B4E9")


host_tree_c <- ggtree(c_tree_r, branch.length = "none") %<+%  c + 
  geom_tiplab(aes(label=label2), offset = 0.15) + 
  geom_nodepoint() + geom_tippoint(aes(color = label2), size=4) + 
  geom_nodelab(aes(label = label), nudge_x = 0.05) +
  geom_treescale(x=1.5, y=1.4, offset = 0.05) + theme_tree2() + xlim(0, 9.5) + theme_tree() +
  scale_color_manual(values = c_cols) 
host_tree_c

cnodes <- ggtree(c_tree_r, branch.length = "none") + geom_text(aes(label=node), hjust=-.3) # for node numbers

# beautify tree
host_tree_c <- ggtree::rotate(host_tree_c, 29)
host_tree_c <- ggtree::rotate(host_tree_c, 17)
host_tree_c <- ggtree::rotate(host_tree_c, 18)
host_tree_c <- ggtree::rotate(host_tree_c, 19)
host_tree_c <- ggtree::rotate(host_tree_c, 23)
host_tree_c <- ggtree::rotate(host_tree_c, 24)

w_clades <- host_tree_c +  geom_cladelabel(node = 23, label = "Complexa", offset = 3, offset.text = 0.15, angle = 90) +
  geom_cladelabel(node = 17, label = "Robusta", offset = 3, offset.text = 0.15, angle = 90)
w_clades



## BRAY-CURTIS
c_bc <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/coral_wCar_den_grouped_bc_upgma-min-2-f.tre")

tiplabels <- c_bc$tip.label 
tiplabels2 <- c("C. foliascens (DR)", "P. damicornis (PI)", "P. verrucosa (RB)", "P. cylindrica (RR)", 
                "P. massive (PI)", "A. hyacinthus (PI)", "A. formosa (PI)", "P. cactus (PI)",
                "P. speciosa (OR)", "S. hystrix (RR)", "S. pistillata (BR)", "S. hystrix (PI)",
                "P. cylindrica (PI)", "D. heliopora (PI)", "E. mammiformis (PI)")

c_b <- data.frame(label = tiplabels, label2 = tiplabels2)


den_c_bc <- ggtree(c_bc) %<+% c_b + geom_tiplab(offset = -.24, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  #geom_nodelab(aes(label = label), nudge_x = -0.02) + # to see if values need formatting
  scale_x_reverse(limits=c(0.8,0)) + theme_tree() +
  scale_color_manual(values = c_cols) 
den_c_bc 

cbnodes <- ggtree(c_bc) + geom_text(aes(label=node), hjust=-.3) # for node numbers

# Beautify tree
den_c_bc <- ggtree::rotate(den_c_bc, 17)
den_c_bc <- ggtree::rotate(den_c_bc, 19)
den_c_bc <- ggtree::rotate(den_c_bc, 25)
den_c_bc <- ggtree::rotate(den_c_bc, 20)

# Format node labels
cb_dat <- den_c_bc$data
cb_dat <- cb_dat[!cb_dat$isTip,]
cb_dat$label <- as.numeric(cb_dat$label)
cb_dat$label <- round(cb_dat$label, 2)

# Add node labels
c_bc_nodes <- den_c_bc + geom_text(data=cb_dat, aes(label=label), nudge_x = -0.02) 


# Combine
cor1 <- ggdraw() + 
  draw_plot(w_clades, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(c_bc_nodes, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.37, p = 0.02 \n\n nRF = 0.69, p < 0.001", size = 12, x = 0.5, y = 0.325)
cor1  

ggsave(filename = "Coral_bc_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## WEIGHTED UF
c_wuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/coral_den_FI-tree_trim_wuf_upgma.tre")

tiplabels <- c_wuf$tip.label 
tiplabels2 <- c("C. foliascens (DR)", "P. cactus (PI)", "P. damicornis (PI)", "P. verrucosa (RB)", 
                "P. cylindrica (RR)", "P. massive (PI)", "A. hyacinthus (PI)", "A. formosa (PI)",  
                "P. speciosa (OR)", "S. pistillata (BR)", "D. heliopora (PI)", "E. mammiformis (PI)", 
                "S. hystrix (RR)", "S. hystrix (PI)", "P. cylindrica (PI)")

c_w <- data.frame(label = tiplabels, label2 = tiplabels2)

den_c_w<- ggtree(c_wuf) %<+% c_w + geom_tiplab(offset = -.18, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  #geom_nodelab(aes(label = label), nudge_x = -0.02) +   
  scale_x_reverse(limits=c(0.475,0)) + theme_tree() +
  scale_color_manual(values = c_cols) 
den_c_w

cwnodes <- ggtree(c_wuf) + geom_text(aes(label=node), hjust=-.3) # for node numbers

# Beautify tree
den_c_w <- ggtree::rotate(den_c_w, 17)
den_c_w <- ggtree::rotate(den_c_w, 18)
den_c_w <- ggtree::rotate(den_c_w, 19)
den_c_w <- ggtree::rotate(den_c_w, 21)
den_c_w <- ggtree::rotate(den_c_w, 22)
den_c_w <- ggtree::rotate(den_c_w, 26)

# Format node labels
cw_dat <- den_c_w$data
cw_dat <- cw_dat[!cw_dat$isTip,]
cw_dat$label <- as.numeric(cw_dat$label)
cw_dat$label <- round(cw_dat$label, 2)

# Add node labels
c_wuf_nodes <- den_c_w + geom_text(data=cw_dat, aes(label=label), nudge_x = -0.02)

cor2 <- ggdraw() + 
  draw_plot(w_clades, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(c_wuf_nodes, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.38, p = 0.01 \n\n nRF = 0.69, p < 0.001", size = 12, x = 0.5, y = 0.325)
cor2  

ggsave(filename = "Coral_wUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## UNWEIGHTED UF
c_uw <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/coral_den_FI-tree_trim_uwuf_upgma.tre")

tiplabels <- c_uw$tip.label 
tiplabels2 <- c("C. foliascens (DR)", "P. speciosa (OR)", "P. cactus (PI)", "E. mammiformis (PI)", 
                "P. massive (PI)", "S. pistillata (BR)", "S. hystrix (PI)", "D. heliopora (PI)", 
                "P. cylindrica (PI)", "S. hystrix (RR)", "P. verrucosa (RB)", "P. cylindrica (RR)", 
                "A. formosa (PI)", "A. hyacinthus (PI)", "P. damicornis (PI)")

c_u <- data.frame(label = tiplabels, label2 = tiplabels2)

den_c_u <- ggtree(c_uw) %<+% c_u + geom_tiplab(offset = -.22, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  #geom_nodelab(aes(label = label), nudge_x = -0.02) +   
  scale_x_reverse(limits=c(0.7,0)) + theme_tree() +
  scale_color_manual(values = c_cols) 
den_c_u

cunodes <- ggtree(c_uw) %<+% c_u + geom_tiplab(offset = -.21, aes(label=label2))+ 
  geom_text(aes(label=node), hjust=-.3) # for node numbers

# Beautify tree
den_c_u <- ggtree::rotate(den_c_u, 25)
den_c_u <- ggtree::rotate(den_c_u, 27)
den_c_u <- ggtree::rotate(den_c_u, 21)
den_c_u <- ggtree::rotate(den_c_u, 24)

# Format node labels
cu_dat <- den_c_u$data
cu_dat <- cu_dat[!cu_dat$isTip,]
cu_dat$label <- as.numeric(cu_dat$label)
cu_dat$label <- round(cu_dat$label, 2)

# Add node labels
c_uwuf_nodes <- den_c_u + geom_text(data=cu_dat, aes(label=label), nudge_x = -0.02)

cor3 <- ggdraw() + 
  draw_plot(w_clades, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(c_uwuf_nodes, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.42, p = 0.03 \n\n nRF = 0.92, p = 0.15", size = 12, x = 0.5, y = 0.325)
cor3  

ggsave(filename = "Coral_uwUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## Octocoral ----

## HOST TREE
o_tree <- read.nexus("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/octocoral_concatenated_msa_curated2.nex.con.tre")
o_tree_r <- root(o_tree, "Carteriospongia_foliascens", resolve.root = T)

# Edit tips
tiplabels <- o_tree_r$tip.label  
tiplabels2 <- c("Briareum sp2 (PR)", "Briareum sp (PR)", "Clavularia sp (PR)",
                "Sinularia sp (RR)", "Sinularia sp2 (PI)", "Sinularia sp (PI)",
                "Sarcophyton sp (RR)", "Sarcophyton sp (PI)", "Cladiella sp (RB)", 
                "Heteroxenia sp (PI)","Pinnigorgia sp (RB)", "I. hippuris (RB)", 
                "C. foliascens (DR)")

o <- data.frame(label = tiplabels, label2 = tiplabels2)

# Set colours
o_cols <- c("Briareum sp2 (PR)" = "#E69F00", "Briareum sp (PR)" = "#E69F00", "Clavularia sp (PR)" = "#D55E00",
            "Sinularia sp (RR)" = "#56B4E9", "Sinularia sp2 (PI)" = "#56B4E9", "Sinularia sp (PI)" = "#56B4E9",
            "Sarcophyton sp (RR)" = "#0072B2", "Sarcophyton sp (PI)" = "#0072B2", "Cladiella sp (RB)" = "#CC79A7", 
            "Heteroxenia sp (PI)" = "#F0E442","Pinnigorgia sp (RB)" = "#009E73", "I. hippuris (RB)" = "#009E73", 
            "C. foliascens (DR)" = "#999999")


host_tree_o <- ggtree(o_tree_r, branch.length = "none") %<+%  o + 
  geom_tiplab(aes(label=label2), offset = 0.15) + 
  geom_nodepoint() + geom_tippoint(aes(color = label2), size=4) + 
  geom_nodelab(aes(label = label), nudge_x = 0.05) +
  geom_treescale(x=1.5, y=1.4, offset = 0.05) + theme_tree2() + xlim(0, 12.5) + theme_tree() +
  scale_color_manual(values = o_cols) 
host_tree_o

onodes <- ggtree(o_tree_r, branch.length = "none") %<+% o + geom_tiplab(offset = -.21, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 

# Beautify tree
host_tree_o <- ggtree::rotate(host_tree_o, 20)
host_tree_o <- ggtree::rotate(host_tree_o, 23)
host_tree_o <- ggtree::rotate(host_tree_o, 21)
host_tree_o <- ggtree::rotate(host_tree_o, 24)



## BRAY CURTIS
o_bc <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/octocoral_wCar_den_grouped_bc_upgma-min-2-f.tre")

tiplabels <- o_bc$tip.label 
tiplabels2 <- c("C. foliascens (DR)", "I. hippuris (RB)", "Clavularia sp (PR)", "Briareum sp (PR)", 
                "Briareum sp2 (PR)", "Sinularia sp (RR)", "Sinularia sp2 (PI)", "Sinularia sp (PI)",
                "Sarcophyton sp (RR)", "Sarcophyton sp (PI)", "Pinnigorgia sp (RB)", "Heteroxenia sp (PI)",  
                "Cladiella sp (RB)")

o_b <- data.frame(label = tiplabels, label2 = tiplabels2)


den_o_bc <- ggtree(o_bc) %<+% o_b + geom_tiplab(offset = -.26, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.8,0)) + theme_tree() +
  scale_color_manual(values = o_cols) 
den_o_bc 

o_b_nodes <- ggtree(o_bc, branch.length = "none") %<+% o + geom_tiplab(offset = -.21, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 

# Beautify tree
den_o_bc <- ggtree::rotate(den_o_bc, 16)
den_o_bc <- ggtree::rotate(den_o_bc, 17)
den_o_bc <- ggtree::rotate(den_o_bc, 18)
den_o_bc <- ggtree::rotate(den_o_bc, 25)


# Combine
oct1 <- ggdraw() + 
  draw_plot(host_tree_o, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_o_bc, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.23, p = 0.02 \n\n nRF = 0.64, p < 0.001", size = 12, x = 0.5, y = 0.325)
oct1  

ggsave(filename = "Octocoral_BC_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## WEIGHTED UF
o_wuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/octocoral_den_FI-tree_trim_wuf_upgma.tre")

tiplabels <- o_wuf$tip.label 
tiplabels2 <- c("C. foliascens (DR)", "I. hippuris (RB)", "Briareum sp2 (PR)", "Briareum sp (PR)", 
                "Clavularia sp (PR)", "Cladiella sp (RB)", "Heteroxenia sp (PI)", "Pinnigorgia sp (RB)",
                "Sinularia sp2 (PI)", "Sinularia sp (PI)", "Sinularia sp (RR)", 
                "Sarcophyton sp (RR)", "Sarcophyton sp (PI)")

o_w <- data.frame(label = tiplabels, label2 = tiplabels2)


den_o_w <- ggtree(o_wuf) %<+% o_w + geom_tiplab(offset = -.16, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.5,0)) + theme_tree() +
  scale_color_manual(values = o_cols) 
den_o_w 

o_w_nodes <- ggtree(o_wuf, branch.length = "none") %<+% o + geom_tiplab(offset = -.21, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 

# Beautify tree
den_o_w <- ggtree::rotate(den_o_w, 19)
den_o_w <- ggtree::rotate(den_o_w, 24)
den_o_w <- ggtree::rotate(den_o_w, 20)
den_o_w <- ggtree::rotate(den_o_w, 16)

# Combine
oct2 <- ggdraw() + 
  draw_plot(host_tree_o, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_o_w, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.36, p < 0.001 \n\n nRF = 0.82, p = 0.02", size = 12, x = 0.5, y = 0.325)
oct2

ggsave(filename = "Octocoral_wUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

# UNWEIGHTED UF
o_uwuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/octocoral_den_FI-tree_trim_uwuf_upgma.tre")

tiplabels <- o_uwuf$tip.label 
tiplabels2 <- c("I. hippuris (RB)", "Clavularia sp (PR)", "Briareum sp (PR)", "Briareum sp2 (PR)", 
                "C. foliascens (DR)", "Sinularia sp (PI)", "Cladiella sp (RB)", "Heteroxenia sp (PI)",
                "Sarcophyton sp (PI)", "Sinularia sp (RR)", "Sarcophyton sp (RR)","Sinularia sp2 (PI)", 
                "Pinnigorgia sp (RB)")

o_u <- data.frame(label = tiplabels, label2 = tiplabels2)


den_o_u <- ggtree(o_uwuf) %<+% o_u + geom_tiplab(offset = -.24, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  #geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.7,0)) + theme_tree() +
  scale_color_manual(values = o_cols) 
den_o_u 

o_u_nodes <- ggtree(o_uwuf, branch.length = "none") %<+% o + geom_tiplab(offset = -.21, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 

# Beautify tree
den_o_u <- ggtree::rotate(den_o_u, 22)
den_o_u <- ggtree::rotate(den_o_u, 24)
den_o_u <- ggtree::rotate(den_o_u, 17)

# Format node labels
o_u_dat <- den_o_u$data
o_u_dat <- o_u_dat[!o_u_dat$isTip,]
o_u_dat$label <- as.numeric(o_u_dat$label)
o_u_dat$label <- round(o_u_dat$label, 2)

# Add formatted node labels
o_uwuf_nodes <- den_o_u + geom_text(data=o_u_dat, aes(label=label), nudge_x = -0.02)


# Combine
oct3 <- ggdraw() + 
  draw_plot(host_tree_o, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(o_uwuf_nodes, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.25, p < 0.001 \n\n nRF = 0.91, p = 0.24", size = 12, x = 0.5, y = 0.325)
oct3  

ggsave(filename = "Octocoral_uwUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## Sponge ----

## HOST TREE
s_tree <- read.nexus("~/Documents/R/Ramaciotti_16S_dataset/Phylogenies/Final_host/sponge_concatenated_msa_curated.nex.con.tre")
s_tree_r <- root(s_tree, "Cladiella_sp_1399", resolve.root = T)

# Edit tips
tiplabels <- s_tree_r$tip.label  
tiplabels2 <- c("Coscinoderma sp (RR)", "I. ramosa (RR)", "C. foliascens (DR)", "Ircinia sp (BR)",
                "I. ramosa (DR)", "C. singaporense (PI)", "Cladiella sp (RB)")
s <- data.frame(label = tiplabels, label2 = tiplabels2)


# Set colours
s_cols <- c("Coscinoderma sp (RR)" = "#56B4E9", "I. ramosa (RR)" = "#0072B2", "C. foliascens (DR)" = "#009E73", 
            "Ircinia sp (BR)" = "#CC79A7", "I. ramosa (DR)" = "#0072B2", "C. singaporense (PI)" = "#E69F00", 
            "Cladiella sp (RB)" = "#999999")


host_tree_s <- ggtree(s_tree_r, branch.length = "none") %<+%  s + 
  geom_tiplab(aes(label=label2), offset = 0.15) + 
  geom_nodepoint() + geom_tippoint(aes(color = label2), size=4, size = 1) + 
  geom_nodelab(aes(label = label), nudge_x = 0.05) +
  geom_treescale(x=1.5, y=1.3, offset = 0.05) + theme_tree2() + xlim(0, 8.75) + theme_tree() +
  scale_color_manual(values = s_cols) 
host_tree_s

snodes <- ggtree(s_tree_r, branch.length = "none") %<+% s + geom_tiplab(offset = -.21, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 


## BRAY CURTIS
s_bc <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/sponge_wCla_den_grouped_bc_upgma-min-2-f.tre")

tiplabels <- s_bc$tip.label 
tiplabels2 <- c("Cladiella sp (RB)", "C. singaporense (PI)", "C. foliascens (DR)", 
                "Coscinoderma sp (RR)", "Ircinia sp (BR)", "I. ramosa (RR)", "I. ramosa (DR)")

s_b <- data.frame(label = tiplabels, label2 = tiplabels2)


den_s_bc <- ggtree(s_bc) %<+% s_b + geom_tiplab(offset = -.30, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.85,0)) + theme_tree() +
  scale_color_manual(values = s_cols) 
den_s_bc 


# Combine
spg1 <- ggdraw() + 
  draw_plot(host_tree_s, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_s_bc, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.71, p < 0.001 \n\n nRF = 0.2, p < 0.001", size = 12, x = 0.5, y = 0.325)
spg1

ggsave(filename = "Sponge_BC_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## WEIGHTED UF
s_wuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/sponge_den_FI-tree_trim_wuf_upgma.tre")

tiplabels <- s_wuf$tip.label 
tiplabels2 <- c("Cladiella sp (RB)", "I. ramosa (RR)", "I. ramosa (DR)", 
                "Coscinoderma sp (RR)", "Ircinia sp (BR)",
                "C. foliascens (DR)", "C. singaporense (PI)")

s_w <- data.frame(label = tiplabels, label2 = tiplabels2)


den_s_w <- ggtree(s_wuf) %<+% s_w + geom_tiplab(offset = -.25, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.65,0)) + theme_tree() +
  scale_color_manual(values = s_cols) 
den_s_w 

s_w_nodes <- ggtree(s_wuf, branch.length = "none") %<+% s + geom_tiplab(offset = .05, aes(label=label2))+ # for node numbers
  geom_text(aes(label=node), hjust=-.3) 

# Beautify tree
den_s_w <- ggtree::rotate(den_s_w, 10)
den_s_w <- ggtree::rotate(den_s_w, 11)
den_s_w <- ggtree::rotate(den_s_w, 13)


# Combine
spg2 <- ggdraw() + 
  draw_plot(host_tree_s, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_s_w, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.78, p = 0.006 \n\n nRF = 0.4, p = 0.01", size = 12, x = 0.5, y = 0.325)
spg2

ggsave(filename = "Sponge_wUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")

## UNWEIGHTED UF
s_uwuf <- read.tree("~/Documents/R/Ramaciotti_16S_dataset/Phylosymbiosis/Dendrograms/Filtered/sponge_den_FI-tree_trim_uwuf_upgma.tre")

tiplabels <- s_uwuf$tip.label 
tiplabels2 <- c("Cladiella sp (RB)", "C. foliascens (DR)", "C. singaporense (PI)",
                "Coscinoderma sp (RR)", "Ircinia sp (BR)", 
                "I. ramosa (RR)", "I. ramosa (DR)")

s_u <- data.frame(label = tiplabels, label2 = tiplabels2)

den_s_u <- ggtree(s_uwuf) %<+% s_u + geom_tiplab(offset = -.26, aes(label=label2)) +
  geom_nodepoint() + geom_tippoint(aes(colour = label2), size = 4) + 
  geom_nodelab(aes(label = label), nudge_x = -0.02) +
  scale_x_reverse(limits=c(0.75,0)) + theme_tree() +
  scale_color_manual(values = s_cols) 
den_s_u


# Combine
spg3 <- ggdraw() + 
  draw_plot(host_tree_s, x = 0, y = .35, width = .5 , height = .6) + draw_label("Host Phylogeny", size=13, x=0.1, y=0.97) +
  draw_plot(den_s_u, x = .5, y = .35, width = .5, height = .6) + draw_label("Microbial dendrogram", size=13, x=0.9, y=0.97) +
  draw_label("Mantel r = 0.75, p = 0.03 \n\n nRF = 0.4, p = 0.01", size = 12, x = 0.5, y = 0.325)
spg3

ggsave(filename = "Sponge_uwUF_final.pdf", device = "pdf", width = 27, height = 25, units = "cm", dpi = "print")


