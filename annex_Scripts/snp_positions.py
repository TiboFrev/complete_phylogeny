vcf = open ("C:\\Users\\frevillet\\Desktop\\fasta_files\\vcf\\62_genomes\\62_2334_snp.vcf")

SNP=[]

####This part extracts the positions of each SNP out of the VCF file##########
line = vcf.readline().split()
while line[0][0] == "#":
    line =  vcf.readline().split() ##Ignoring the commentary lines 

SNP = SNP + [line[1]]
while line:
        line = vcf.readline().split()
        if len(line)>1:
            SNP = SNP + [line[1]] ##Every position is stored in a list
##############################################################################


####This part extract the sequence of the aligned reference genome (removing the \n)##########
fasta = open("C:\\Users\\frevillet\\Desktop\\fasta_files\\62_genomes\\referencegenomev62.fasta","r")
ligne=fasta.readline()
compteur = 0
seq=""
while (ligne !=""):
    ligne = fasta.readline()
    for i in range(len(ligne)):
        if ligne[i] != "\n":
            seq=seq+(ligne[i])
##############################################################################################


###This part calculate the relative position using the reference genome########          
pos = []
for i in range(len(seq)):
    if seq[i] in ["a","t","c","g"]:
        compteur+=1
    if str(i) in SNP :
        for i in range(SNP.count(str(i))):
            pos.append(compteur)
###############################################################################

###Each position is now stored in the corresponding 1k subpart ################
SNPnbr = {}
for i in range (0,(max(pos)//1000)+1):
    SNPnbr[i*1000] = len(list(filter(lambda n: i*1000<=n<(i+1)*1000, pos)))
###############################################################################


#######################creating the file#######################################    
output = open("C:\\Users\\frevillet\\Desktop\\SNPpositions.txt","w")
for j in range(len(SNPnbr.keys())):
    output.write(str(list(SNPnbr.keys())[j]) + " : " + str(SNPnbr[list(SNPnbr.keys())[j]]) + "\n")           
output.close()
###############################################################################
