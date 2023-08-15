library(seqinr)

args <- commandArgs(trailingOnly = TRUE)
AC_of_interest_file <- read.table(args[1])
fasta <- args[2]
output <- args[3]

AC_of_interest <- c()
for (i in (1:length(AC_of_interest_file[[1]]))){
  AC_of_interest <- c(AC_of_interest,AC_of_interest_file[[1]][i])
}

data<-read.fasta(fasta)
pattern <- "^[A-Za-z]{2}_\\d+\\.\\d+|^[A-Za-z]{2}_\\d+"
vec.tokeep <- c()
for (i in 1:length(data)){
  if (grepl(pattern, attr(data[[i]],"name"))){
    AC <- paste(strsplit(attr(data[[i]], "name"), "_")[[1]][1], strsplit(attr(data[[i]], "name"), "_")[[1]][2],sep = "_")
  }
  else{
    AC <- strsplit(attr(data[[i]], "name"), "_")[[1]][1]
  }
  if (!AC %in% AC_of_interest) {
    vec.tokeep <- c(vec.tokeep, i)
  }
}
write.fasta(sequences = data[vec.tokeep], names = names(data)[vec.tokeep], file.out = output)