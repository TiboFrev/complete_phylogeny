##Import the necessary libraries
library(restez)
library(rentrez)


##This script takes 3 arguments:
##  A string "species" which is the name of the organism of interest
##  A string "title" which is the name of the run, used to name files
##  A list of NCBI's GenBank sequence IDs
##It returns a fasta file with headers containing at least the accession number, the country and the date of collection


##Arguments collection
args <- commandArgs(trailingOnly = TRUE)
species <- args[1] 
title <- args[2] 
IDs <- args[3:length(args)]


##Filtering the IDs to avoid collecting useless information
Toremove <- readLines("IDstoremove.txt") #This file contains the accession numbers of duplicated sequences
for (i in 1:length(Toremove)){
  a<-entrez_summary(db = "nuccore", id = Toremove[i])$uid #For each accession number, get the corresponding ID
  if (a %in% IDs){ 
    IDs <- setdiff(IDs, a) #Remove the duplicates' IDs from the variable IDs
  }
}


##Accessing with restez the country and the collection date for each ID.
metadata <- sapply(IDs, function(ID){
  #First : get the information of the genbank file in text.
  suppressMessages(genbank_data <- entrez_fetch(db="nucleotide", id=ID,rettype = "text", retmode="text"))
  
  #Extracting the "features" part
  features <- gb_extract(genbank_data, what = c("features"))
  
  #Extracting the collection date (4 possible ways):
    # if there is a "collection_date" feature :
    collection_date <- features[[1]]$collection_date
  
    #if there is a collected_by feature
    if (is.null(collection_date) && !is.null(features[[1]]$collected_by) ){
      collection_date <- features[[1]]$collected_by
    }
    
    #if the date is in the "note" feature
    else if(is.null(collection_date)&& !is.null(features[[1]]$note)){
      if(length(features[[1]]$note[grepl("(19|20)[0-9]{2}", features[[1]]$note)])!=0){
        collection_date <- regmatches(features[[1]]$note[grepl("(19|20)[0-9]{2}", features[[1]]$note)], regexpr("(19|20)[0-9]{2}", features[[1]]$note))[1]}
      }
    
    #if the date is in the "isolate" feature
    else if(is.null(collection_date)&& !is.null(features[[1]]$isolate)){
      if(length(features[[1]]$isolate[grepl("(19|20)[0-9]{2}", features[[1]]$isolate)])!=0){
        collection_date <- regmatches(features[[1]]$isolate[grepl("(19|20)[0-9]{2}", features[[1]]$isolate)], regexpr("(19|20)[0-9]{2}", features[[1]]$isolate))[1]
      }
    
    }
  #Extraction of the country
  country <- features[[1]]$country
  
  print(paste ("Collecting the information for ", ID,  ", please wait"))
  print(country)
  print(collection_date)
  return(c(country,collection_date))
})
print("Information successfully collected")


##Opening the file, created previously in the pipeline by CGaccession.R, which contains the sequences
##Opening a new file to create new headers
new_fasta <- "transientfile.fasta"
old_fasta = paste(title,"full_genomes.fasta",sep = "_")
file_conn <- file(old_fasta, "r") 
output_conn <- file(new_fasta, "w")

# Iterate the fasta line by line
i <- 1
while (length(line <- readLines(file_conn, n = 1)) > 0) {
  if (startsWith(line, ">")) {#At each sequence header
    modified_line <- paste(line, gsub("c\\(|\\)", "", metadata[i]), sep = "_")#Add the metadata information
    i <- i+1
  } else {
    modified_line <- line
  }
  writeLines(modified_line, output_conn)
}

close(file_conn)
close(output_conn)

##Now suppress the unnecessary information in the headlines
fasta <- readLines("transientfile.fasta")
fasta_modified <- gsub("complete genome|isolate|genomic sequence|_NULL| strain|\"|\'|:|,|\\)|\\(", "", fasta) #Removing non relevant information
fasta_modified2 <- gsub(" |  ", "_", fasta_modified) #Getting rid of the spaces, for the tree labels later
writeLines(fasta_modified2, old_fasta)
file.remove("transientfile.fasta")
