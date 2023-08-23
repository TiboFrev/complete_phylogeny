###This script takes in argument a sequence alignment and returns a csv file containing the percentage of gaps in each aligned sequence
###It can be used to find poor quality sequences in an alignment for example

##Import the required packages
import csv

##Initialisation
alignment = open("C:\\Users\\frevillet\\Desktop\\LSDV_GOATPOX_SHEEPOX.fasta","r") #alignment in which you want to count the gaps for each sequence
ligne=alignment.readline()
titre = ligne
length=0
gaps=0
summary={}

while ligne:
    if ligne[0] == ">": #each time you encounter a new sequence (header), enter the loop
        titre = ligne[0:-1]
        titre = titre.replace('"','')
        ligne = alignment.readline()
        while ligne and ligne[0]!=">": #Iterate through the entire sequence
            for i in range(len(ligne)):
                if ligne[i] != "\n":
                    length += 1 #total length
                    if ligne[i] == "-":
                        gaps+=1 #number of gaps
            ligne=alignment.readline()
        summary[titre] = (((gaps*100)/length))  #associate each sequence with its gaps percentage 
        length = 0
        gaps = 0
        titre=""
alignment.close()

with open('gapsPOXVIRUSES.csv', 'w') as csvfile: #summarise the information in a csv file
    filewriter = csv.writer(csvfile, delimiter=',',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['Genome', 'gaps %'])
    for i in range(len(summary)):
        filewriter.writerow([(list(summary.keys())[i]), summary[list(summary.keys())[i]]])
    
   
