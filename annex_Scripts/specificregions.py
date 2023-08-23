###This script takes in argument a list of positions (integers) and a fasta file
###It returns several fasta files with the sequences of the specified regions of the input fasta file

##Import the necessary package
import os

##Argument collection
#each element is one edge, so one region is two elements (eg 135021-135095)
regs = [135021,135095,143865,144132,116144,118620,133724,134435] #regions to extract
fasta = open("C:\\Users\\frevillet\\Desktop\\fasta_files\\66_genomes\\66_LSDV_genomes.fasta", "r") #fasta file from which to extract

ligne = fasta.readline()
titre = ligne
seq = ""

while ligne:
    if ligne[0] == ">": #enter the loop each time you read a new sequence
        titre = ligne
        ligne = fasta.readline()
        while ligne and ligne[0] != ">": #Iterate through the entire sequence 
            for i in range(len(ligne)):
                if ligne[i] != "\n":
                    seq = seq + (ligne[i]) #extract only the actual sequence (without the \n on each line)
            ligne = fasta.readline()
        subseqs = []
        for k in range(0, len(regs), 2): #for each regions
            start = regs[k]
            end = regs[k + 1]
            if seq[start:end] != "":
                subseqs.append([seq[start:end]]) #extract the specific region
                
        for j in range(len(subseqs)): #write the regions in a fasta file
                fastaa = open("C:\\Users\\frevillet\\Desktop\\blast\\trees66\\" + str(regs[2*j]+30) + ".fasta", "a")
                fastaa.write(titre + "".join(subseqs[j]) + "\n\n")
                fastaa.close()
        seq = ""
        titre = ""
fasta.close()
