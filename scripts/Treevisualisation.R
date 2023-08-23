##Import the necessary libraries
library(ape)
library(phytools)

##This script create 6 different picture files. 3 pdf and 3 png, which represent the phylogenetic trees of one sequence alignment
##It takes two arguments:
##  A tree file in the Newick format
##  A string "title" which is the name of the run, applied to the title of the picture files

##Argument collection
args <- commandArgs(trailingOnly = TRUE)
BestTree <- args[1]
title <- args[2]

##Reading and rooting the tree
tree <- read.tree(BestTree)
tree<-midpoint.root(tree)


############################## Creation of the different files #######################################


##Phylogramm with labels
#png file
png(file = paste(title,"phylogram.png",sep = "_"), width = 1900, height = 900, units = "px")
plot(tree, cex = 0.8, show.tip.label = TRUE, type = "phylogram", no.margin = FALSE)
ape::add.scale.bar(cex = 0.8)
dev.off() 


#pdf file
pdf(file =paste(title,"phylogram.pdf",sep = "_"), width = 100, height = 50)
plot(tree, cex = 4,cex.main = 4, show.tip.label = TRUE, type = "phylogram", no.margin = FALSE)
ape::add.scale.bar(cex=4)
dev.off()


##Fan with labels
#png file
png(file = paste(title,"fan.png",sep = "_"), width = 1900, height = 900, units = "px")
plot(tree, cex = 0.8, show.tip.label = TRUE, type = "fan", no.margin = FALSE)
ape::add.scale.bar(cex = 0.8)
dev.off()

#pdf file
pdf(file =paste(title,"fan.pdf",sep = "_"), width = 100, height = 50)
plot(tree, cex = 4,cex.main = 4, show.tip.label = TRUE, type = "fan", no.margin = FALSE)
ape::add.scale.bar(cex=4)
dev.off()


##Phylogramm without labels
#png file
png(file = paste(title,"phylogram_nolabel.png",sep = "_"), width = 1900, height = 900, units = "px")
plot(tree, cex = 0.8, show.tip.label = FALSE, type = "phylogram", no.margin = FALSE)
ape::add.scale.bar(cex = 0.8)
dev.off()  

#pdf file
pdf(file =paste(title,"phylogram_nolabel.pdf",sep = "_"), width = 100, height = 50)
plot(tree, cex = 4,cex.main = 4, show.tip.label = FALSE, type = "phylogram", no.margin = FALSE)
ape::add.scale.bar(cex=4)
dev.off()


######################################################################################################
