args <- commandArgs(trailingOnly = TRUE)

dups <- readLines(args)
pattern <- "^[A-Za-z]{2}\\d+\\.\\d+$|^[A-Za-z]{2}_\\d+\\.\\d+$|^[A-Za-z]\\d{5}\\.\\d+$"


IDs <- c()
for (i in 1:length(dups)){
  for (j in 2:length(strsplit(dups[i], ",")[[1]])){
    IDs<-c(IDs,(strsplit(dups[i], ",")[[1]][j]))
  }
}
print(IDs)
newIDs <- sapply(IDs, function(x) gsub(" ", "", x))

writeLines(newIDs, "IDstoremove.txt")