##Packages to install if necessary

# update.packages(ask=F)
knitr::opts_chunk$set(echo = T)
# from https://github.com/gtonkinhill/fastbaps
#BiocManager::install("ggtree", ask=F)
#BiocManager::install("phytools", ask=F)
#BiocManager::install("ggplot2", ask=F)
#BiocManager::install("treeio", ask=F)
#install.packages("ape", ask=F)
#install.packages("devtools", ask=F)
#install.packages("irlba", ask=F)

##Import the necessary libraries
library(irlba)
library(ggtree)   # ggtree v2.4.1
library(phytools) # phytools v0.7-70
library(ggplot2)  # ggplot v2_3.3.3
library(ape)      # ape v5.5
library(treeio)   # treeio v1.14.3
library(grid)
library(Rcpp)
library(RcppArmadillo)
#devtools::install_github("gtonkinhill/fastbaps", force=T, type="source")
library(fastbaps) # Load libraries


##This script uses the FastBAPS algorithm from Tonkin-Hill et al. For more information see github.com/gtonkinhill/fastbaps
##It takes in argument:
##  An alignment file in fasta format
##  A tree file in Newick format
##And it performs a Hierarchical Bayesian  clustering analysis

##Arguments collection
args <- commandArgs(trailingOnly = TRUE)
fastaname = args[[1]]
treename = args[[2]]

##Initialisation
samples = c()   
fastbapsplot = paste0("FastBAPS_", fastaname, ".png", sep="")
heatmapplot = paste0("Heatmap_", fastaname, ".png", sep="")
sparse.data <-import_fasta_sparse_nt(fastaname)


# do clustering, get dk values, then Bayesian clustering
sparse.data <- optimise_prior(sparse.data, type="optimise.symmetric")

best.partition <- best_baps_partition(sparse.data, fast_baps(sparse.data))
# partition model

##Plot the fastbaps clustering analysis
plot.df <-data.frame(id=colnames(sparse.data$snp.matrix),
                        fastbaps=best.partition, stringsAsFactors=F)
gg <- ggtree(read.tree(treename),linewidth=0.2) + geom_tiplab(hjust=0.06, size=0.8) 
   
png(fastbapsplot, height=980, width=980, units='px', res=300) # tree & groups
print(facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile,
                    aes(x=fastbaps), color="blue")) 
dev.off()

##Plot the heatmap of resampling analysis   
dendro <- as.dendrogram(fast_baps(sparse.data)) # Heatmap
   
png(heatmapplot, height=2000, width=2000, units='px', res=300) #  
print(gplots::heatmap.2(boot_fast_baps(sparse.data), dendro, dendro,
                   tracecol=NA, margins=c(8,8), key=F, cex.lab=0.1))
dev.off() 
