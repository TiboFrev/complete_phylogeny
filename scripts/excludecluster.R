library(ape)
library(stringr)

arg <- commandArgs(trailingOnly = TRUE)
tree <- read.tree(arg[[1]])
moyenne <- mean(tree$edge.length)
SD <- sd(tree$edge.length)
index <- which(tree$edge.length >= (moyenne + 10 * SD))

if (length(index)==0){
  write("", file = "IDstoremove.txt")
} else {

for (i in 1:length(index)) {
  NODE <- c(tree$edge[index[i], 2])
  subtree <- extract.clade(tree,tree$edge[index[i], 2])
  png(file = paste("subtree_",i,".png"), width = 1900, height = 1000, units = "px")
  plot(subtree)
  dev.off() 
  for (k in 1:length(subtree$tip.label)) {
    NODE[k]<- subtree$tip.label[[k]]
    NODE[k]<- str_extract(NODE[k],"[A-Za-z]{2}\\d+\\.\\d+|[A-Za-z]{2}_\\d+\\.\\d+|[A-Za-z]\\d{5}\\.\\d+")
  }
  write(paste(NODE, collapse = "\n"), file = "IDstoremove.txt", append = TRUE)
}

}
