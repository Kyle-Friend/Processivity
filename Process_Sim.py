"""Kyle Friend at Washington and Lee University
Two files are written here. One contains columns with RPI values, and the other
contains columns with simulated ribosome densities versus the fraction of the open
reading frame."""

import numpy as np
all_genes = {}
n = 0

"""The Mock mRNAs file contains uniform densities corresponding to every mRNA in
mouse ESCs with densities based on translational status. The first operation
generates a dictionary for analysis."""
with open('Mock_mRNAs.txt', 'r') as infile:
    for line in infile:
        line = line.rstrip()
        if n % 2 == 0:
            gene = line
        else:
            line = line.split()
            line = [float(x) for x in line]
            all_genes[gene] = line
        n += 1

outfile = open('RPIs_Tens.txt', 'w')
probs = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10] #Change these values to change the likelihood that a ribosome will abort translation.
genes = list(all_genes.keys())
total = np.zeros((100, 10))
n = 0
print(len(genes))
for gene in genes:
    dens = all_genes[gene]
    if len(dens) < 100: continue
    for j in range(len(probs)):
        temp = []
        for i in range(len(dens)):
            if i == 0:
                temp.append(dens[i])
                continue
            num = np.random.randint(0, probs[j])
            if num == 1:
                val = (10)/len(dens)
                val = temp[-1] - val
                if val <= 0:
                    temp.append(0)
                else:
                    temp.append(val)
            else:
                temp.append(temp[-1])
        RPI = 0 #This is where RPI is calculted for every gene.
        for i in range(len(temp)):
            RPI = RPI + (i+1)*temp[i]
        RPI = RPI/(sum(temp)*len(temp))
        outfile.write(str(RPI) + '\t')
        block = len(temp)//100 #This operation is used to chunk all codons into 100 separate bins.
        temp2 = []
        for i in range(0, len(temp), block):
            if i + block > len(temp) - 1:
                temp2.append(sum(temp[i:]))
            else:
                temp2.append(sum(temp[i:i+block]))
        while(len(temp2)) > 100:
            num = np.random.randint(len(temp2))
            temp2.pop(num)
        for i in range(len(temp2)):
            total[i][j] = total[i][j] + temp2[i]
    outfile.write('\n')
    n += 1
    if n % 100 == 0:
        print(n)

outfile.close()
outfile2 = open('Tens_Proc.txt', 'w')
np.savetxt(outfile2, total, delimiter = '\t', newline = '\n')
outfile2.close()
