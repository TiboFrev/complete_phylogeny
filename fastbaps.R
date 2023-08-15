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

samples = c()
args <- commandArgs(trailingOnly = TRUE)
fastaname = args[[1]]
treename = args[[2]]
   
fastbapsplot = paste0("FastBAPS_", fastaname, ".png", sep="")
heatmapplot = paste0("Heatmap_", fastaname, ".png", sep="")

sparse.data <-import_fasta_sparse_nt(fastaname)

# do clustering, get dk values, then Bayesian clustering
sparse.data <- optimise_prior(sparse.data, type="optimise.symmetric")
# hyperparameter: 0.044 for core
best.partition <- best_baps_partition(sparse.data, fast_baps(sparse.data))
# partition model
plot.df <-data.frame(id=colnames(sparse.data$snp.matrix),
                        fastbaps=best.partition, stringsAsFactors=F)
gg <- ggtree(read.tree(treename),linewidth=0.2) + geom_tiplab(hjust=0.06, size=0.8) 
   
png(fastbapsplot, height=980, width=980, units='px', res=300) # tree & groups
print(facet_plot(gg, panel="fastbaps", data=plot.df, geom=geom_tile,
                    aes(x=fastbaps), color="blue")) 
dev.off()
   
dendro <- as.dendrogram(fast_baps(sparse.data)) # Heatmap
   
png(heatmapplot, height=2000, width=2000, units='px', res=300) #  
print(gplots::heatmap.2(boot_fast_baps(sparse.data), dendro, dendro,
                   tracecol=NA, margins=c(8,8), key=F, cex.lab=0.1))
dev.off() # end function plotstuff