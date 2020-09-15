from itertools import *

f = open("nucleotide_alignment.fa", "r")

ntseqs = [] # nucleotide sequences
ntgaps = [] # gaps in nucleotide sequences
ntnames = []
oldseq = ""
name = ""
coordinate = 0
for x in f.readlines():
    if x[0] == ">":
        name = x[1:-1]
        if oldseq != "":
            ntseqs.append(oldseq)
            ntnames.append(name)
            oldseq = ""
            coordinate = 0
    if x[0] != ">":
        for p in x:
            if p != "\n":
                coordinate = coordinate + 1
                if p == "-":
                    ntgaps.append(coordinate)
                oldseq = oldseq + p

ntseqs.append(oldseq)
ntnames.append(name)

ntgaps.sort()
f.close()

#get rid of sites with > 50% gaps
positions_to_trim = []
for i in count(0):
    if i > len(ntseqs[0]):
        break
    if len(list(compress(range(len(ntgaps)), [x == i for x in ntgaps]))) >= 0.5*len(ntseqs):
        positions_to_trim.append(i)

ntseqs_trim = []

for s in ntseqs:
    l = list(s)
    for i in sorted(positions_to_trim, reverse=True):
        del l[i]
    ntseqs_trim.append("".join(l))

f = open("aminoacid_alignment.fa", "r")

aaseqs = [] # amino acid sequences
aanames = []
name = ""
coordinate = 0
for x in f.readlines():
    if x[0] == ">":
        name = x[1:-1]
        aanames.append(name)
        if oldseq != "":
            aaseqs.append(oldseq)
            oldseq = ""
    if x[0] != ">":
        for p in x:
            if p != "\n":
                oldseq = oldseq + p

f.close()
i = 0
matches = []
for p in range(len(aaseqs[0]) - 5):
    for s,n in zip(aaseqs,aanames):
        query = s[p:p+5]
        # when I added this line I got an error
        middle_aa = query[3]
        # looking for GGDEF or GGEEF motif
        if query[:2] == "GG" and middle_aa in ["D","E"] and query[4:] == "EF":
            matches.append(n + " " + str(p) + " " + middle_aa)
    i = i + 1

f = open("matches.txt", "w")

for match in matches:
    f.write(match + "\n")

f.close()
