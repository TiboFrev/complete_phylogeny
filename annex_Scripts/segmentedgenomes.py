###This script takes in argument one alignment, one file containing the aligned version of a reference genome, and one integer
###It returns n fasta files containing the alignment for n successive regions.
###For example, if cut = 1000, creates n files of aligned 1kb regions

##Argument collection
refgenome = open("C:\\Users\\frevillet\\Desktop\\fasta_files\\66_genomes\\referencegenomev66.fasta","r")
cut = int(input("Size of the wanted segment of genomes :"))
fasta=open("C:\\Users\\frevillet\\Desktop\\fasta_files\\66_genomes\\MSA66_LSDV_genomes.fasta","r")


#######################################Calculating the relative position of each cutting position#################################################
ligne=refgenome.readline()
compteur = 0
seq=""

##Extracting the sequences without the /n at the end of each line
while (ligne !=""): 
    ligne = refgenome.readline()
    for i in range(len(ligne)):
        if ligne[i] != "\n":
            seq=seq+(ligne[i])
            
positions=[]

##Keeping each relative position in the list "positions"
for i in range(len(seq)):
    if seq[i] in ["a","t","c","g","-"]:
        compteur+=1
        if compteur%cut == 0 :
            positions.append(i)

refgenome.close()
################################################################################################################################################


############################################Cutting the alignment at each cutting position######################################################
ligne=fasta.readline()
titre = ligne
seq=""

while ligne:
    if ligne[0] == ">": ##Iterate through each aligned genome. You first get the n sections of the first genome, then the n sections of the second genome, etc.
        titre = ligne
        ligne = fasta.readline()
        while ligne and ligne[0]!=">": ##Getting the sequence WITHOUT the /n in the end of each line
            for i in range(len(ligne)):
                if ligne[i] != "\n":
                    seq=seq+(ligne[i])
            ligne=fasta.readline()
        subseqs=[]
        for k in range(0, len(positions)-1): ##Cutting the sequence in equal portions of length = cut
            if seq[positions[k]:positions[k+1]] !="":
                subseqs.append ( [seq[positions[k]:positions[k+1]]] )
        if seq[positions[-1]:] != "": ##Get the last nucleotides if there are some (length < cut)
            subseqs.append([seq[positions[-1]:]])
        subseqs = [seq[0:positions[0]]] + subseqs #Get the first "cut" nucleotides

        ##Actually writing the n subsequences in n different fasta files
        for j in range(len(subseqs)): 
            fastaa = open ("C:\\Users\\frevillet\\Desktop\\1k_66_subseq\\" + str(j*cut)+".fasta", "a")
            fastaa.write(titre + "".join(subseqs[j]) + "\n\n")
            fastaa.close()
        seq = ""
        titre=""
fasta.close()

################################################################################################################################################

