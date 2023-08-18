library(ape)

args <- commandArgs(trailingOnly = TRUE)
BestTree <- args[1]
species <- args[2]

tree <- read.tree(BestTree)

png(file = paste(species,"phylogeny.png",sep = "_"), width = 1900, height = 900, units = "px")

plot(tree, cex = 0.8, show.tip.label = TRUE, type = "phylogram", main = paste("Phylogeny of", species), no.margin = FALSE)
ape::add.scale.bar(cex = 0.8)
dev.off()  

pdf(file =paste(species,"phylogeny.pdf",sep = "_"), width = 100, height = 50)

plot(tree, cex = 4,cex.main = 4, show.tip.label = TRUE, type = "phylogram", main = paste("Phylogeny of", species), no.margin = FALSE)
ape::add.scale.bar(cex=4)
dev.off()

png(file = paste(species,"phylogeny_nolabel.png",sep = "_"), width = 1900, height = 900, units = "px")

plot(tree, cex = 0.8, show.tip.label = FALSE, type = "phylogram", main = paste("Phylogeny of", species), no.margin = FALSE)
ape::add.scale.bar(cex = 0.8)
dev.off()  

pdf(file =paste(species,"phylogeny_nolabel.pdf",sep = "_"), width = 100, height = 50)

plot(tree, cex = 4,cex.main = 4, show.tip.label = FALSE, type = "phylogram", main = paste("Phylogeny of", species), no.margin = FALSE)
ape::add.scale.bar(cex=4)
dev.off()