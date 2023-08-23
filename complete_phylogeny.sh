
#!/usr/bin/env Rscript

help(){
echo "this function collect the complete genome sequences of an organism (or several organisms) in the NCBI nuccore database, align them and create a tree file (using mafft and RAxML)"
echo "Typical usage : ./complete_phylogeny.sh -s ORG1,ORG2,ORGn -t NAME_OF_RUN -r -d"
echo ""
echo "Options :"
echo "-s is the species option. write each species you want the full genomes from, separated by a comma. This argument is mandatory. Example : -s LSDV,'Sheeppox virus',ASFV will get the sequences of LSDV SPV and AFSV"
echo "-t is the name of your run, to name the files. This argument is mandatory"
echo "-r Having the -r option removes the annex files created by RAxML"
echo "-d Having the -d option lets only one step of distance checking in the tree. Useful if you have several species. Not recommanded if you only have one species"
echo "-h display this help message"
}

RAXML="/mnt/lustre/RDS-live/freville/LSDV/RAxML/raxml-ng" ##Change the path to yours. Using raxml-ng

##Initialisation
ORG=() 
REMOVE_FILES="FALSE"
dist_check="TRUE"


##options management
while getopts 's:drh:t:' opt ; do
case $opt in
  s) ORG_STRING="$OPTARG" ;;
  r) REMOVE_FILES="TRUE" ;; 
  h) help
     exit 0  
     ;;
  d) dist_check="FALSE" ;;
  t) TITLE="$OPTARG" ;;
  \?)
     echo""
     echo "Invalid option: -$OPTARG" >&2
     echo""
     help
     exit 1
     ;;

esac
done

# Checking the recquired arguments
if [ -z "$ORG_STRING" ] || [ -z "$TITLE" ] ; then
  echo "Missing required arguments. Please refer to the help message (-h) for usage."
  exit 1
fi


##Divide the ORG_STRING in a table
IFS=',' read -ra ORG <<< "$ORG_STRING"  

##Accession via GenBank to the sequences and the metadata
module load R/latest
for value in "${ORG[@]}"; do #iterate thorugh each species
   value=${value// /_}
   echo "Processing : $value"
   IDs=$(Rscript CGaccession.R "$value" $TITLE) #This accesses the sequences in a fasta file (see next line)
   fastafile="${TITLE}_full_genomes.fasta"
   
   seqkit rmdup -s <"$fastafile"> "${fastafile}2" -D duplicatedseq.txt #Removing the duplicates from the fasta file
   Rscript duplicatedIDs.R "duplicatedseq.txt" #Getting the accession numbers of the duplicates in IDstoremove.txt
   
   rm "$fastafile"
   rm duplicatedseq.txt
   mv "${fastafile}2" "$fastafile"

   Rscript infocollection.R $value $TITLE $IDs #collection of the country and date of collection
done
echo "A multi fasta file name ${fastafile} has been created in your working directory"


##first multi-alignement
mafft --auto "$fastafile" > "${TITLE}_MSA.fasta" 
MSAfile="${TITLE}_MSA.fasta"

##first tree construction
$RAXML --model GTR+G --msa "$MSAfile" 
TREEFILE="${TITLE}_MSA.fasta.raxml.bestTree"

##Creating the result folder
result_folder="${TITLE}"
mkdir -p "$result_folder"
  

#################################### BEGINNING OF THE DISTANCE CHECKING LOOP #################################################  
 
breakpoint="NO"
passage=0

while [ "$breakpoint" = "NO" ];
do
  if [ -s IDstoremove.txt ]; then ##
    rm "IDstoremove.txt" 
  fi
  let "passage=passage+1"
  echo $passage
  
  Rscript excludecluster.R "$TREEFILE" "$passage" ##Checking the distance to extract the divergent genomes
  echo "exclude cluster done"
  if [ -s IDstoremove.txt ]; then ##If there are genomes to extract
    ./excludeID.sh -t "IDstoremove.txt" -i $fastafile -o "${passage}_${TITLE}_reduced.fasta" #Remove them from the fasta file
    echo "exclude ID done"
    
    find . -name "${MSAfile}.raxml.*" -not -name "${MSAfile}.raxml.bestTree" -exec rm {} \;
    
    fastafile="${passage}_${TITLE}_reduced.fasta"
    
    ##Align the new fasta file (without the divergent ones)
    echo "alignment starting"
    mafft --auto "$fastafile" > "${passage}_${TITLE}_reduced_MSA.fasta"
    MSAfile="${passage}_${TITLE}_reduced_MSA.fasta"
    
    
    ##Create the new tree
    echo "raxml starting"
    $RAXML --model GTR+G --msa "$MSAfile" 
    TREEFILE="${MSAfile}.raxml.bestTree"
    
    ##If -d option, exit the loop (only one distance check to avoid excluding too much, or separating two species)
    if [ "$dist_check" = "FALSE" ]; then
      breakpoint="YES"
    fi
  
  else
    breakpoint="YES"
  fi
done  

  
###################################### ENDING OF THE DISTANCE CHECKING LOOP ##################################################


##Creating 6 tree image files 
Rscript Treevisualisation.R "$TREEFILE" $TITLE


##FASTBAPS ANALYSIS##

#Fastbaps requires small names for the genomes. This line replace the entire header by the only accession number in the alignment file
awk '/^>/{match($0, /[A-Za-z]{2}[0-9]+\.[0-9]+|[A-Za-z]{2}_[0-9]+\.[0-9]+|[A-Za-z][0-9]{5}\.[0-9]+/, arr); print ">" arr[0]; next} 1' "$MSAfile" > "${TITLE}_simplified_MSA.fasta"
$RAXML --model GTR+G --msa "${TITLE}_simplified_MSA.fasta" #create the tree with simplified headers
Rscript fastbaps.R "${TITLE}_simplified_MSA.fasta" "${TITLE}_simplified_MSA.fasta.raxml.bestTree" #Run the fastbaps algorithm




###################################### MOVING THE RESULT FILES IN THE RESULT FOLDER ##################################################

mv "${TITLE}_simplified_MSA.fasta" "$result_folder"
mv "${TITLE}_full_genomes.fasta" "$result_folder"
mv "${TITLE}_MSA.fasta" "$result_folder"
mv "${TITLE}_MSA.fasta.raxml.bestTree" "$result_folder"
for file in *_"$TITLE"_reduced.fasta; do
  mv "$file" "$result_folder"
done
for file in *_"$TITLE"_reduced_MSA.fasta; do
  mv "$file" "$result_folder"
done
mv "${TITLE}_phylogram.png" "$result_folder"
mv "${TITLE}_phylogram.pdf" "$result_folder"
mv "${TITLE}_fan.pdf" "$result_folder"
mv "${TITLE}_fan.png" "$result_folder"
mv "${TITLE}_phylogram_nolabel.pdf" "$result_folder"
mv "${TITLE}_phylogram_nolabel.png" "$result_folder"
mv "FastBAPS_${TITLE}_simplified_MSA.fasta.png" "$result_folder"
mv "Heatmap_${TITLE}_simplified_MSA.fasta.png" "$result_folder"
mv subtree_* "$result_folder"
  
find . -name "*_${TITLE}_*_.raxml.*" -exec mv {} "$result_folder" \;
find . -name "${TITLE}_simplified_MSA.fasta.raxml.*" -exec mv {} "$result_folder" \;
  

##If option -r  
if [ "$REMOVE_FILES" = "TRUE" ]; then ##remove the RAxML annex files 
  cd $result_folder
  find . -name "*_${TITLE}_*_.raxml.*" -not -name "*_${TITLE}_*_.raxml.bestTree" -exec rm {} \;
  find . -name "${TITLE}_simplified_MSA.fasta.raxml.*" -not -name "${TITLE}_simplified_MSA.fasta.raxml.bestTree" -exec rm {} \;

fi

######################################################################################################################################






