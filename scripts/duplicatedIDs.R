
##This script takes in argument the output text file of seqkit::rmdup, containing accession numbers
##of duplicated sequences
##It returns a text file containing the accession numbers of the removed sequences
##with one accession number per lines


##Arguments collection
args <- commandArgs(trailingOnly = TRUE)
dups <- readLines(args) 

##Extracting the IDs in a new variable
IDs <- c()
for (i in 1:length(dups)){ #each line is one set of n duplicated sequences
  for (j in 2:length(strsplit(dups[i], ",")[[1]])){ #The first ID is the kept sequence, all the others need to be stored
    IDs<-c(IDs,(strsplit(dups[i], ",")[[1]][j]))
  }
}
print(IDs)
newIDs <- sapply(IDs, function(x) gsub(" ", "", x))

##Creating the output file
writeLines(newIDs, "IDstoremove.txt")
