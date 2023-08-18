args <- commandArgs(trailingOnly = TRUE)

dups <- readLines(args)


IDs <- c()
for (i in 1:length(dups)){
  for (j in 2:length(strsplit(dups[i], ",")[[1]])){
    IDs<-c(IDs,(strsplit(dups[i], ",")[[1]][j]))
  }
}
print(IDs)
newIDs <- sapply(IDs, function(x) gsub(" ", "", x))

writeLines(newIDs, "IDstoremove.txt")
