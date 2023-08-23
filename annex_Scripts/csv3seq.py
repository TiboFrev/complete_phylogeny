###This script extract information from the output file of the recombination detection tool 3seq.
###It has been created for specific analysis for LSDV, so it might require big changes to be adapted for another work


##Import the required packages
import csv
import re
from pprint import pprint
import sys

##Sort each genome to his belonging clade
Clade13 = set(["OM984486.1","MW355944.1","OP922506.1","OM105589.1","OP508345.1","OP752701.1","OQ267778.1","OQ267777.1","ON152411.1","MW732649.1","OP985536.1","OL752713.2","MZ577075.1","MZ577074.1","OM803091.1","OM803092.1","MZ577073.1","MZ577076.1","OM793603.1","OM793602.1","OM984485.1"])
Clade12 = set(["MN995838.1","KY702007.1","MT643825.1","KY829023.3","KX894508.1","MN642592.1","MW030512.1","MH893760.2","MW699032.1","OQ588787.1","AF409137.1","MW656253.1","OK318001.1","MW631933.1","KX683219.1","OP688129.1","OK422494.1","OP688128.1","OP297402.1","MN072619.1","NC_003027.1"])
Clade11 = set(["OM793605.1","OM793607.1","OM793606.1","OM793608.1", "MK441838.1", "AF409138.1", "OM793609.1", "KX764645.1", "MG972412.1", "KX764643.1", "KX764644.1", "MW656252.1", "MW435866.1", "OM793604.1", "MN636839.1","MN636843.1", "MN636842.1", "MN636838.1"])
Clade2 =set(["OM373209.1", "ON400507.1", "OK422492.1" ,  "OK422493.1"])
R = set(["MT134042.1","OL542833.1","MH646674.1","OM530217.1"])
vaccine = set(["MK441838.1","AF409138.1","OM793609.1","KX764645.1","MG972412.1","KX764643.1","KX764644.1","KX683219.1","MW631933.1"])

L=[]
freq = {}

##Small function which takes an iterator and return TRUE if all the elements are equals
def all_equal(iterator):
    iterator = iter(iterator)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)


##opening the CSV file
with open('3seq.62.csv', mode ='r')as file:
    
    ##reading the CSV file
    csvFile = csv.reader(file)
    summ = {}
    next(csvFile)
    for line in csvFile:
        line = line[0].split()
        #extract the accession numbers of parents (N1 and N2) and child (N3) genomes with regex
        N1 =(re.findall("[A-Za-z]{2}\\d+\\.\\d+|[A-Za-z]{2}_\\d+\\.\\d+|[A-Za-z]\\d{5}\\.\\d+",line[0])[0]) ##N1 is the first parent genome
        N2 =(re.findall("[A-Za-z]{2}\\d+\\.\\d+|[A-Za-z]{2}_\\d+\\.\\d+|[A-Za-z]\\d{5}\\.\\d+",line[1])[0]) ##N2 is the second parent genome
        N3 =(re.findall("[A-Za-z]{2}\\d+\\.\\d+|[A-Za-z]{2}_\\d+\\.\\d+|[A-Za-z]\\d{5}\\.\\d+",line[2])[0]) ##N3 is the child genome
       
        ##Creating a dictionnary with a parent couple as a key, and the child genome and the recombination regions as values
        if (N1,N2) not in summ :
            summ[(N1,N2)] = {"child" : [[N3, line[12],line[14]]]}
        else:
            summ[(N1,N2)]["child"].append([N3, line[12],line[14]])

        ##Creating a list L of regions of recombinations
        ##Creating a dictionnary freq of the frequency of the regions in the results
        if [line[12],line[14]] not in L: #L contains the regions, and freq the frequencies of these regions in recombination results
            L.append([line[12],line[14]])
            freq[str(line[12])+"-"+str(line[14])]=1
        else:
            freq[str(line[12])+"-"+str(line[14])]+=1

for i in summ: ##simplify the dictionnary in case one recombination led to several child genomes
    LISTE=[]
    LISTE2=[]
    for j in summ[i]["child"]:
        LISTE.append((j[1],j[2]))
    if all_equal(LISTE):
        for j in summ[i]["child"]:
            LISTE2.append(j[0])
        summ[i] = {"region":(j[1],j[2]),"child":LISTE2}

##count the numbers of recombination between genome from Clade 1.1 with genome from Clade 1.1, Clade 1.1 with 1.2 etc.
unun=0
undeux=0
untrois=0
unr=0
deuxdeux=0
deuxtrois=0
deuxr=0
troistrois=0
troisr=0
rr=0
c2un=0
c2deux=0
c2trois=0
c2r=0
c2c2=0

print("Vaccine strains implicated :")
for key in summ.keys():
    for i in key:
        if i in vaccine:
            print(key)
            print(str(len(summ[key]["child"]))+" recombination(s)")

print("\n")
for key in summ.keys():
    c11 = sum(elem in Clade11 for elem in key)
    c12 = sum(elem in Clade12 for elem in key)
    c13 = sum(elem in Clade13 for elem in key)
    c2 = sum(elem in Clade2 for elem in key)
    r = sum(elem in R for elem in key)
    
    print(f"Tuple {key}:")
    print(f"clade11: {c11} accession(s)")
    print(f"clade12: {c12} accession(s)")
    print(f"clade13: {c13} accession(s)")
    print(f"recombinants: {r} accession(s)")
    print(f"clade2: {c2} accession(s)")
    a = (len(summ[key]["child"]))
    print("for "+str(a)+" recombinations\n")

    if c11 == 1 and c12==1:
        undeux+=a
    elif c13==1 and r==1:
        troisr+=a
    elif c12 ==2:
        deuxdeux+=a
    elif c11==1 and r==1:
        unr+=a
    elif c11==1 and c13==1:
        untrois+=a
    elif c12==1 and c13==1:
        deuxtrois+=a
    elif c12==1 and r==1:
        deuxr += a
    elif c11 == 2:
        unun+=a
    elif c13==2:
        troistrois+=a
    elif r==2:
        rr+=a
    elif c2 == 2:
        c2c2+=a
    elif c2 == 1 and c11 == 1:
        c2un+=a
    elif c2 == 1 and c12 == 1:
        c2deux+=a
    elif c2 == 1 and c13 == 1:
        c2trois+=a
    elif c2 == 1 and r == 1:
        c2r+=a
print("Clade\t1.1\t1.2\t1.3\tr\t2\n1.1  \t"+str(unun)+"\t"+str(undeux)+"\t"+str(untrois)+"\t"+str(unr)+"\t"+str(c2un)+"\n1.2  \t"+str(undeux)+"\t"+str(deuxdeux)+"\t"+str(deuxtrois)+"\t"+str(deuxr)+"\t"+str(c2deux)+"\n1.3  \t"+str(untrois)+"\t"+str(deuxtrois)+"\t"+str(troistrois)+"\t"+str(troisr)+"\t"+str(c2trois)+"\nr    \t"+str(unr)+"\t"+str(deuxr)+"\t"+str(troisr)+"\t"+str(rr)+"\t"+str(c2r)+"\n2  \t"+str(c2un)+"\t"+str(c2deux)+"\t"+str(c2trois)+"\t"+str(c2r)+"\t"+str(c2c2)+"\n")
print(str(len(L))+"unique regions of recombination\n")

##calculate the real position in the genome (using one aligned genome as reference)
fasta = open("C:\\Users\\frevillet\\Desktop\\fasta_files\\62_genomes\\referencegenomev62.fasta","r")
ligne=fasta.readline()
compteur = 0
seq=""
while (ligne !=""):
    ligne = fasta.readline()
    for i in range(len(ligne)):
        if ligne[i] != "\n":
            seq=seq+(ligne[i])
M=[]
for i in range(len(L)):
    for j in range(len(L[i])):
        M.append(str.split(L[i][j],sep='-')[0])
        M.append(str.split(L[i][j],sep='-')[1])

pos = {}
for i in range(len(seq)):
    if seq[i] in ["a","t","c","g"]:
        compteur+=1
    if str(i) in M :
        pos[str(i)]=compteur    
for i in range(len(M)):
    M[i] = pos[M[i]]

    
Q=[]
for i in range(0,len(M),4):
    Q.append(str(M[i])+"-"+str(M[i+1])+"-"+str(M[i+2])+"-"+str(M[i+3]))

freq = dict(zip(Q,list(freq.values())))
print("Recombination frequency of each unique couple:\n")
print(freq)
print("\nSummary of each unique couple and their child recombinant genome:")
pprint(summ)

dico={}
for i in summ:
    for j in i:
        if j not in dico:
            dico[j]=len(summ[i]["child"])
        else:
            dico[j]=dico[j]+len(summ[i]["child"])

