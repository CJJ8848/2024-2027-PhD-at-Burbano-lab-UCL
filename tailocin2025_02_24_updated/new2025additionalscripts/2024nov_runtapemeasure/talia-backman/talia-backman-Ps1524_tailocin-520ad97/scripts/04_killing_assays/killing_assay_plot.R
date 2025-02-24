library(ggplot2)
library(ape)
library(phytools)
library(ggtree)
library(dplyr)
library(viridis)
setwd('/Users/cuijiajun/Desktop/others/tmphernan/2024nov_runtapemeasure/talia-backman/talia-backman-Ps1524_tailocin-520ad97/')
# read in data
# pseudomonas phylogeny
tree <- read.tree("./data/phylogeny_data/ps_1524_uncollapsed_5_2018.nwk")
# killing assay dat
dat2 <- read.csv("./data/killing_assay_data/3replicates_killing_assay_dat.csv")

# subset tree to only the tester strains used in dat
uneeded <- which(!(tree$tip.label) %in% dat2$Tester_Strain)
tree <- drop.tip(tree, uneeded)

# clean dat2
# remove all rows with NAs 
dat2 <- na.omit(dat2)
# subset dat2 to only the data we need
dat2 <- dat2[,c(1,2,7)]
# make binary column
dat2$Final_qualitative[dat2$Final_qualitative == "Resistant"] <- 0
dat2$Final_qualitative[dat2$Final_qualitative == "Sensitive"] <- 1
dat2$Final_qualitative[dat2$Final_qualitative == "Partially Sensitive"] <- 0
dat2$Final_qualitative[dat2$Final_qualitative == "Mostly Sensitive"] <- 1
# convert values in results to numeric
dat2$Final_qualitative = as.numeric(dat2$Final_qualitative)
# go from long to wide format for plotting
wide_dat2 <- reshape(dat2, idvar = "Tester_Strain", timevar = "Tailocin_donor_strain", direction = "wide")
# make first column row names, then remove that column
rownames(wide_dat2) <- wide_dat2[,1]
wide_dat2 <- wide_dat2[ ,c(2:16)]
# clean column names 
colnames(wide_dat2) <- sub("Final_qualitative.","", colnames(wide_dat2))
# convert to matrix for plotting
dat_matrix2 <- as.matrix(wide_dat2)

# plot with a heatmap on the phylogeny
p <- ggtree(tree) + geom_tiplab(size = 8)
p3 <- gheatmap(p, dat_matrix2, colnames = TRUE, legend_title = "Results", offset = .05, color = "black", 
               colnames_angle = 30, colnames_position = "top", font.size = 8,hjust = 0, width = 7) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 + scale_fill_viridis(discrete = FALSE)

# remove no MMC controls from plot
wide_dat2 <- wide_dat2[ ,c(1,3,5,7,9,11)]
dat_matrix2 <- as.matrix(wide_dat2)

# plot without MMC controls
p <- ggtree(tree) + geom_tiplab(size = 8)
p3 <- gheatmap(p, dat_matrix2, colnames = TRUE, legend_title = "Results", offset = .05, color = "black", 
               colnames_position = "top", font.size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 + scale_fill_viridis(discrete = FALSE)

# remove bad data (killed itself: p23.B8 and p26.D6) from plot
wide_dat2 <- wide_dat2[,c(1,2,3,6)]
dat_matrix2 <- as.matrix(wide_dat2)
# plot
p <- ggtree(tree) + geom_tiplab(size = 8)
p3 <- gheatmap(p, dat_matrix2, colnames = TRUE, legend_title = "Results", offset = .05, color = "black", 
               colnames_position = "top", font.size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 + scale_fill_viridis(discrete = FALSE)

#only p25.A12
dat_matrix2 <- dat_matrix2[, "p25.A12", drop = FALSE]

# without any labels
p <- ggtree(tree) 
p3 <- gheatmap(p, dat_matrix2, colnames = FALSE, legend_title = "Results", offset = .05, color = "black", 
               colnames_position = "top", font.size = 8) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p3 + scale_fill_viridis(discrete = FALSE)

# export as pdf
pdf("killing_assaysbycui.pdf")
p3 + scale_fill_viridis(discrete = FALSE)
dev.off()
