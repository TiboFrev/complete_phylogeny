
#!/usr/bin/env Rscript

help(){
echo "this function collect the complete genome sequences of an organism in the NCBI nuccore database, align them and create a tree file (using mafft and RAxML)"
echo ""
echo "Options :"
echo "-s is the species option. -s LSDV will collect the sequences of LSDV. This option is required"
echo "-r = TRUE or FALSE (default). It removes the annex files created by RAxML when equals to TRUE"
echo "-h display this help message"
}

RAXML="/mnt/lustre/RDS-live/freville/LSDV/RAxML/raxml-ng" ##Change the path to yours 



while getopts ':s:r:h' opt ; do
case $opt in
  s) ORG="$OPTARG" ;;
  r) REMOVE_FILES="$OPTARG" ;;
  h) help
     exit 0  
     ;;
  \?)
     echo""
     echo "Invalid option: -$OPTARG" >&2
     echo""
     help
     exit 1
     ;;


esac
done

ORG=${ORG// /_}

module load R/latest
IDs=$(Rscript CGaccession.R $ORG)

echo $IDs


fastafile="full_${ORG}_genomes.fasta"

seqkit rmdup -s <"$fastafile"> "${fastafile}2" -D duplicatedseq.txt

Rscript duplicatedIDs.R "duplicatedseq.txt"

rm "$fastafile"
mv "${fastafile}2" "$fastafile"

Rscript infocollection.R $ORG $IDs 

mafft --auto "$fastafile" > "full_${ORG}_genomesMSA.fasta"

MSAfile="full_${ORG}_genomesMSA.fasta"
$RAXML --model GTR+G --msa "$MSAfile" 
TREEFILE="full_${ORG}_genomesMSA.fasta.raxml.bestTree"

breakpoint="NO"
passage=0
while [ "$breakpoint" = "NO" ];
do
  if [ -s IDstoremove.txt ]; then
    rm "IDstoremove.txt" 
    echo "removing done"
  fi
  let "passage=passage+1"
  echo $passage
  
  Rscript excludecluster.R "$TREEFILE" "$passage"
  echo "exclude cluster done"
  if [ -s IDstoremove.txt ]; then
    ./excludeID.sh -t "IDstoremove.txt" -i $MSAfile -o "reduced_${ORG}_genomesMSA.fasta"
    echo "exclude ID done"
  
  echo "raxml starting"
  $RAXML --model GTR+G --msa "reduced_${ORG}_genomesMSA.fasta" --prefix "${ORG}_${passage}_tree"
  MSAfile="reduced_${ORG}_genomesMSA.fasta"
  TREEFILE="${ORG}_${passage}_tree.raxml.bestTree"
  
  else
    breakpoint="YES"
  fi
done  

  

  Rscript Treevisualisation.R "$TREEFILE" $ORG


  awk '/^>/{match($0, /[A-Za-z]{2}[0-9]+\.[0-9]+|[A-Za-z]{2}_[0-9]+\.[0-9]+|[A-Za-z][0-9]{5}\.[0-9]+/, arr); print ">" arr[0]; next} 1' "$MSAfile" > "simplified_reduced_${ORG}_genomesMSA.fasta"

  $RAXML --model GTR+G --msa "simplified_reduced_${ORG}_genomesMSA.fasta"

  Rscript fastbaps.R "simplified_reduced_${ORG}_genomesMSA.fasta" "simplified_reduced_${ORG}_genomesMSA.fasta.raxml.bestTree"


  result_folder="${ORG}_folder"
  mkdir -p "$result_folder"

  mv "$fastafile" "$result_folder"
  mv "full_${ORG}_genomesMSA.fasta" "$result_folder"
  mv ""$MSAfile"" "$result_folder"
  mv "simplified_reduced_${ORG}_genomesMSA.fasta" "$result_folder"
  mv "${ORG}_phylogram.png" "$result_folder"
  mv "${ORG}_phylogram.pdf" "$result_folder"
  mv "${ORG}_fan.pdf" "$result_folder"
  mv "${ORG}_fan.png" "$result_folder"
  mv "${ORG}_phylogram_nolabel.pdf" "$result_folder"
  mv "${ORG}_phylogram_nolabel.png" "$result_folder"
  mv "FastBAPS_simplified_reduced_${ORG}_genomesMSA.fasta.png" "$result_folder"
  mv "Heatmap_simplified_reduced_${ORG}_genomesMSA.fasta.png" "$result_folder"
  mv subtree_* "$result_folder"
  
  find . -name "ASFV_*_tree.raxml.*" -exec mv {} "$result_folder" \;
  find . -name "simplified_reduced_${ORG}_genomesMSA.fasta.raxml.*" -exec mv {} "$result_folder" \;
  find . -name "full_${ORG}_genomesMSA.fasta.raxml.*" -exec mv {} "$result_folder" \;
  
  if [ "$REMOVE_FILES" = "TRUE" ]; then ##remove the RAxML annex files 
      cd $result_folder
      find . -name "ASFV_*_tree.raxml.*" -not -name "ASFV_${passage}_tree.raxml.bestTree" -exec rm {} \;
      find . -name "simplified_reduced_${ORG}_genomesMSA.fasta.raxml.*" -not -name "simplified_reduced_${ORG}_genomesMSA.fasta.raxml.bestTree" -exec rm {} \;
      find . -name "full_${ORG}_genomesMSA.fasta.raxml.*" -not -name "full_${ORG}_genomesMSA.fasta.raxml.bestTree" -exec rm {} \;
      
  else
      cd $result_folder
      find . -name "ASFV_*_tree.raxml.*" -not -name "ASFV_${passage}_tree.raxml.*" -exec rm {} \;
  fi
