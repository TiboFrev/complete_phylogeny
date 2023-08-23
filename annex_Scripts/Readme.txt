###

These scripts are miscellaneous scripts which were used to analyse my different results on LSDV phylogeny. Some are really specific but can be adapted for some other works.


###Countinggaps.py

Create a csv file with the percentage of gaps in each sequence of a multi-alignment file

###csv3seq.py

Extract different informations from the output csv file created by the recombination detection tool 3seq.

###specificregions.py

Extract from a multifasta the subsequences corresponding to specific regions

###segmentedgenomes.py

From an alignment, and with one aligned genome as reference, create n alignment corresponding to each region of length x.
For example, from an alignment, create n alignments of each 1kb portion of the genomes.

###snp_positions.py

From a vcf file, extract the number of SNP per 1kb in the genomes, and summarise the results in a text file

###snpplot.R

Creates a density plot of the SNPs across the genome. Also add regions of recombination 
