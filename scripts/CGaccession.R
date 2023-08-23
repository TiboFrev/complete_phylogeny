##Import the necessary libraries

library(rentrez)


##This script takes in argument a string "species" which is the name of an organism
##and a string "title" which is the name of the run (names of the files)
##It returns a fasta file with every complete genomic sequences of the input organism available in Genbank

##Arguments collection
args <- commandArgs(trailingOnly = TRUE)
species <- args[1] #ex : LSDV, Sheeppox_virus
title <- args[2] #ex : 22_08_23

##Collecting the IDs of the sequences
sp <- paste (species, "[organism] AND complete genome [title]") ##Query used to ask the database
query <- entrez_search(db="nuccore", term=sp ,retmax=99999, use_history = TRUE) #extracting the IDs of every sequence that matches the query
IDs <- query$ids  #storing the IDs

##Collecting the sequences
sequences <- rentrez::entrez_fetch(db="nuccore", web_history=query$web_history, rettype="fasta") #With the IDs, extracting the fasta sequences
titre = paste(title,"full_genomes.fasta",sep = "_") #Name of the fasta file
write(sequences, titre, sep="\n",append=TRUE) #one fasta file is created, with every single sequence


cat(IDs)
flush.console() #The variable IDs is stored in the next steps, in the shell script (in a text file called IDs.txt)
