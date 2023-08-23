##Import the necessary libraries
library(seqinr)

##This R script removes sequences from a fasta file, based on the input accession numbers given
##It takes in argument:
##  A text file containing one accession number per line
##  An input fasta file
##  An output name for the resulting fasta file.

##Arguments collection
args <- commandArgs(trailingOnly = TRUE)
AC_of_interest_file <- read.table(args[1])
fasta <- args[2]
output <- args[3]

##Extracting the accession numbers from the file to put them in a R vector
AC_of_interest <- c()
for (i in (1:length(AC_of_interest_file[[1]]))){
  AC_of_interest <- c(AC_of_interest,AC_of_interest_file[[1]][i])
}

##initialisation
data<-read.fasta(fasta)
pattern <- "^[A-Za-z]{2}_\\d+\\.\\d+|^[A-Za-z]{2}_\\d+" #Regular expression for a NCBI refseq accession number (with a _ between the letters and the numbers)
vec.tokeep <- c() #Vector which will contain the accession number we need to keep in the fasta

##Iterate through the fasta file
for (i in 1:length(data)){
  
  ##Each time you get to a header with a refseq accession number
  if (grepl(pattern, attr(data[[i]],"name"))){
    AC <- paste(strsplit(attr(data[[i]], "name"), "_")[[1]][1], strsplit(attr(data[[i]], "name"), "_")[[1]][2],sep = "_") ##Extract the accession number (by extracting the two first element, separated by a "_")
  }
  else{##Else, when the header contains a classic accession number without a "_"
    AC <- strsplit(attr(data[[i]], "name"), "_")[[1]][1]##extract the accession number (the first element)
  }
  if (!AC %in% AC_of_interest) {
    vec.tokeep <- c(vec.tokeep, i) ##Keep only the accession numbers that are NOT present in the AC_of_interest file
  }
}
write.fasta(sequences = data[vec.tokeep], names = names(data)[vec.tokeep], file.out = output) ##Write a new fasta without the sequences to remove
