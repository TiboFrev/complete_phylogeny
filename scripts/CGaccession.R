library(rentrez)

args <- commandArgs(trailingOnly = TRUE)
species <- args
sp <- paste (species, "[organism] AND complete genome [title]")

query <- entrez_search(db="nuccore", term=sp ,retmax=99999) #extracting the IDs of every sequence that matches the query
IDs <- query$ids #storing the IDs


sequences <- rentrez::entrez_fetch(db="nuccore", id=IDs, rettype="fasta") #With the IDs, extracting the fasta sequences
titre = paste("full",species,"genomes.fasta",sep = "_")
write(sequences, titre, sep="\n") #one fasta file is created, with every single sequence


cat(IDs)
flush.console() #The variable IDs is stored in the next steps, in the shell script (in a text file called IDs.txt)