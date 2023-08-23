### complete_phylogeny

complete_phylogeny is a shell script which uses mafft (v7.520) and RAxML-ng (v8.2.X), and which accesses complete genomes of one or several organisms available in Genbank to create phylogenetic trees and clustering analysis with FastBAPS
This shell script uses every other scripts available in TiboFrev/complete_phylogeny/scripts

How to use : ./complete_phylogeny.sh -s ORG1,ORG2,ORGn -t NAME_OF_RUN"
Options : -r TRUE removes the annex files created by RAxML-ng
          -h display help message
          -d Having the -d option lets only one step of distance checking in the tree. Useful if you have several species. Not recommanded if you only have one species

### CGaccession.R

Simple R script which, from the name of the organism, return a fasta file every sequence corresponding to the query  "species [organism] AND complete genome [title]"
and a list of the sequences IDs in the variable "IDs"

###duplicatedIDs.R

Just an annex function which collects the accession number of duplicated sequences (found by seqkit rmdup), and stores them in a file called "IDstoremove.txt

###infocollection.R

R script which uses the restez API to access the country and date of collection of the genomes corresponding to the "IDs", and then modify the headers of the fasta file
to fit these information.

###excludecluster.R

Excludecluster.R is a Rscript which uses the tree file (created by RAxML-ng) to identify the nodes in the tree which are very divergent (tree$edge.length >= (mean(edge.length) + 10 * SD(edge.length)))
and identify the accession numbers of the trees that are in these subtrees. They are then sent to excludeID.sh

###excludeID.sh

This shell basically is here to collect the arguments for the Rscript excludeID.R. It takes in argument a text file with the accession numbers that are to remove, an input fasta file, and an output name.

###exludeID.R

This R script is the functional part of the previous script. It uses regular expression to identify the accession numbers in each header of the fasta file, and proceed to the removal of 
the corresponding sequences.

###Treevisualisation.R

This just creates 6 png/pdf files of phylogenetic trees, using the R package ape

###fastbaps.R

This creates a clustering analysis of the tree files using the Fastbaps algorithm (see Tonkin-Hill et al. 2016)
