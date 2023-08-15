
#!/usr/bin/env Rscript


help() {
  echo "This function removes the sequences corresponding to particular headers in a multifasta file."
  echo "Usage: ./excludeID.sh -t headers.txt -i inputfile.fasta -o outputfile.fasta"
  echo "The headers.txt is a text file with one accession number per line (WITH THE VERSION e.g. NW631933.1) and an empty line in the end."
  echo ""
  echo "Options:"
  echo "-t: Path to the text file with the headers"
  echo "-i: Path to the input multifasta file"
  echo "-o: Name of the output multifasta file"
  echo "-h: Display this help message"
}

while getopts 't:i:o:h' opt ; do
  case $opt in
    h) help
       exit 0  
       ;;
    t) IDs="$OPTARG" ;;
    i) inputfile="$OPTARG" ;;
    o) outputfile="$OPTARG" ;;
    \?)
       echo ""
       echo "Invalid option: -$OPTARG" >&2
       echo ""
       help
       exit 1
       ;;
  esac
done


# Checking the recquired arguments
if [ -z "$IDs" ] || [ -z "$inputfile" ] || [ -z "$outputfile" ]; then
  echo "Missing required arguments. Please refer to the help message (-h) for usage."
  exit 1
fi

# checking that the input files exist
if [ ! -f "$IDs" ]; then
  echo "Headers file not found: $IDs"
  exit 1
fi

if [ ! -f "$inputfile" ]; then
  echo "Input multifasta file not found: $inputfile"
  exit 1
fi


#remove the sequences

module load R/latest

Rscript excludeIDs.R $IDs $inputfile $outputfile