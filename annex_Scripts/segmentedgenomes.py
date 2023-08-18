refgenome = open("C:\\Users\\frevillet\\Desktop\\fasta_files\\66_genomes\\referencegenomev66.fasta","r")
cut = int(input("Size of the wanted segment of genomes :"))
ligne=refgenome.readline()
compteur = 0
seq=""
while (ligne !=""):
    ligne = refgenome.readline()
    for i in range(len(ligne)):
        if ligne[i] != "\n":
            seq=seq+(ligne[i])

positions=[]

for i in range(len(seq)):
    if seq[i] in ["a","t","c","g","-"]:
        compteur+=1
        if compteur%cut == 0 :
            positions.append(i)

refgenome.close()

fasta=open("C:\\Users\\frevillet\\Desktop\\fasta_files\\66_genomes\\MSA66_LSDV_genomes.fasta","r")

ligne=fasta.readline()
titre = ligne
seq=""

while ligne:
    if ligne[0] == ">":
        titre = ligne
        ligne = fasta.readline()
        while ligne and ligne[0]!=">":
            for i in range(len(ligne)):
                if ligne[i] != "\n":
                    seq=seq+(ligne[i])
            ligne=fasta.readline()
        subseqs=[]
        for k in range(0, len(positions)-1):
            if seq[positions[k]:positions[k+1]] !="":
                subseqs.append ( [seq[positions[k]:positions[k+1]]] )
        if seq[positions[-1]:] != "":
            subseqs.append([seq[positions[-1]:]])
        subseqs = [seq[0:positions[0]]] + subseqs
        for j in range(len(subseqs)):
            fastaa = open ("C:\\Users\\frevillet\\Desktop\\1k_66_subseq\\" + str(j*cut)+".fasta", "a")
            fastaa.write(titre + "".join(subseqs[j]) + "\n\n")
            fastaa.close()
        seq = ""
        titre=""
fasta.close()


