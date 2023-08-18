### complete_phylogeny

complete_phylogeny is a shell script which uses mafft (v7.520) and RAxML-ng (v8.2.X), and which accesses complete genomes of an organism available in Genbank to create phylogenetic trees and clustering analysis
with FastBAPS
This shell script uses every other scripts available in TiboFrev/complete_phylogeny

How to use : ./complete_phylogeny -s [name of the organism]
Other option : -r TRUE removes the annex files created by RAxML-ng

### CGaccession.R

Simple R script which, from the name of the organism, return a fasta file every sequence corresponding to the query  "species [organism] AND complete genome [title]"
and a list of the sequences IDs in the variable "IDs"

###duplicatedIDs.R

Just an annex function which collects the ID of duplicated sequences (found by seqkit rmdup), in order to remove them from the variable "IDs"

###infocollection.R

R script which uses the restez API to access the country and date of collection of the genomes corresponding to the "IDs".

###
