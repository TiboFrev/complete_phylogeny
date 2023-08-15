library(rentrez)

args <- commandArgs(trailingOnly = TRUE)
species <- args
sp <- paste (species, "[organism] AND complete genome [title]")

query <- entrez_search(db="nuccore", term=sp ,retmax=99999)

IDs <- query$ids


sequences <- rentrez::entrez_fetch(db="nuccore", id=IDs, rettype="fasta")
titre = paste("full",species,"genomes.fasta",sep = "_")
write(sequences, titre, sep="\n")


cat(IDs)
flush.console()