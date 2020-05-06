"""Kyle Friend, Washington and Lee University, April 2020"""

import numpy as np
from scipy import stats

#The first section loads a text file with gene names associated with ribosome densities.
#Data are then stored in a dictionary.
RDPs = {}
n = 0
rpfs = []
with open('Ribosome_densities.txt', 'r') as infile:
    for line in infile:
        line = line.rstrip()
        if n % 5 == 0:
            gene = line
            rpfs = []
        else:
            line = line.split()
            line = [float(x) for x in line]
            count = 0
            for val in line:
                if val == 0:
                    count += 1
            if count > 0.25*len(line):
                pass
            else:
                rpfs.append(line)
        if len(rpfs) == 4:
            RDPs[gene] = rpfs
            rpfs = []
        n += 1

#Next RDI values are calculated for the two conditions.
all_genes = list(RDPs.keys())
for gene in all_genes:
    rpfs = RDPs[gene]
    norms = []
    for rpf in rpfs:
        rdi = []
        for i in range(len(rpf)):
            rdi.append(rpf[i]*(i+1))
        rdi = sum(rdi)/(len(rdi)*sum(rpf))
        norms.append(rdi)
    norms = [(norms[0] + norms[1])/2, (norms[2] + norms[3])/2]
    RDPs[gene] = norms

#Then a file is written with the data. Control treated cells in left column.
#Torin1-treated cells in the right column.
#Lastly, a Mann-Whitney test is done to calculate significance.

cons = []
tors = []
outfile = open('All_RDIs.txt', 'w')
for gene in all_genes:
    vals = RDPs[gene]
    cons.append(vals[0])
    tors.append(vals[1])
    outfile.write(str(gene) + '\t' + str(vals[0]) + '\t' + str(vals[1]) + '\n')
outfile.close()

(U_stat, p_value) = stats.mannwhitneyu(cons, tors)
print(p_value)
