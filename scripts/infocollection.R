library(restez)
library(rentrez)

args <- commandArgs(trailingOnly = TRUE)
species <- args[1]
IDs <- args[2:length(args)]

Toremove <- readLines("IDstoremove.txt")

for (i in 1:length(Toremove)){
  a<-entrez_summary(db = "nuccore", id = Toremove[i])$uid
  if (a %in% IDs){
    IDs <- setdiff(IDs, a)
  }
}

titre = paste("full",species,"genomes.fasta",sep = "_")

suppressMessages(org <- entrez_fetch(db="nucleotide", id=IDs[[1]],rettype = "text", retmode="text"))
organism <- gb_extract(org, what = c("organism"))

metadata <- sapply(IDs, function(ID){
  suppressMessages(genbank_data <- entrez_fetch(db="nucleotide", id=ID,rettype = "text", retmode="text"))
  
  features <- gb_extract(genbank_data, what = c("features"))
  
  collection_date <- features[[1]]$collection_date
  if (is.null(collection_date) && !is.null(features[[1]]$collected_by) ){
    collection_date <- features[[1]]$collected_by
  }
  if(is.null(collection_date)&& !is.null(features[[1]]$note)){
    if(length(features[[1]]$note[grepl("(19|20)[0-9]{2}", features[[1]]$note)])!=0){
      collection_date <- regmatches(features[[1]]$note[grepl("(19|20)[0-9]{2}", features[[1]]$note)], regexpr("(19|20)[0-9]{2}", features[[1]]$note))[1]}
  }
  if(is.null(collection_date)&& !is.null(features[[1]]$isolate)){
    if(length(features[[1]]$isolate[grepl("(19|20)[0-9]{2}", features[[1]]$isolate)])!=0){
      collection_date <- regmatches(features[[1]]$isolate[grepl("(19|20)[0-9]{2}", features[[1]]$isolate)], regexpr("(19|20)[0-9]{2}", features[[1]]$isolate))[1]
    }
    
  }
  country <- features[[1]]$country
  
  print(paste ("Collecting the information for ", ID,  ", please wait"))
  print(country)
  print(collection_date)
  return(c(country,collection_date))
})

print("Information successfully collected")
new_fasta <- "transientfile.fasta"

file_conn <- file(titre, "r")
output_conn <- file(new_fasta, "w")

# Iterate the fasta line by line
i <- 1
while (length(line <- readLines(file_conn, n = 1)) > 0) {
  if (startsWith(line, ">")) {
    # if line starts with ">", add the metadata
    modified_line <- paste(line, gsub("c\\(|\\)", "", metadata[i]), sep = "_")
    i <- i+1
  } else {
    modified_line <- line
  }
  
  writeLines(modified_line, output_conn)
}

close(file_conn)
close(output_conn)

##Now suppress the unecessary information in the headlines
fasta <- readLines("transientfile.fasta")
fasta_modified <- gsub("complete genome|isolate|genomic sequence|_NULL| strain|\"|:|,|\\)|\\(", "", fasta)
fasta_modified2 <- gsub(species,"",fasta_modified)
fasta_modified3 <- gsub(organism,"",fasta_modified2)
fasta_modified4 <- gsub(" |  ", "_", fasta_modified3)
writeLines(fasta_modified4, titre)
file.remove("transientfile.fasta")

print(paste("The multi-fasta", titre, "has been created in your working directory"))