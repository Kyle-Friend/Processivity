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

#Next biological replicates are normalized to the maximum value.
all_genes = list(RDPs.keys())
for gene in all_genes:
    rpfs = RDPs[gene]
    norms = []
    for rpf in rpfs:
        top = max(rpf)
        rpf = [x/top for x in rpf]
        norms.append(rpf)
    RDPs[gene] = norms

#Then arrays are created with data for the two treatments.
#The array split operation creates 100 bins with ribosomal densities.

Con_Codons = []
Tor_Codons = []
for gene in all_genes:
    rpfs = RDPs[gene]
    Con_avg = []
    Tor_avg = []
    if len(rpfs[0]) < 100: continue
    for i in range(len(rpfs[0])):
        Con_avg.append((rpfs[0][i] + rpfs[1][i])/2)
        Tor_avg.append((rpfs[2][i] + rpfs[3][i])/2)
    Con_avg = np.array_split(Con_avg, 100)
    Tor_avg = np.array_split(Tor_avg, 100)
    bin_con = []
    bin_tor = []
    for val in Con_avg:
        bin_con.append(np.mean(val))
    for val in Tor_avg:
        bin_tor.append(np.mean(val))
    Con_Codons.append(bin_con)
    Tor_Codons.append(bin_tor)

Con_Codons = np.array(Con_Codons)
Tor_Codons = np.array(Tor_Codons)

#Maximum values are created to set the final output to 1.
Con_dense = np.mean(Con_Codons, axis = 0)
Tor_dense = np.mean(Tor_Codons, axis = 0)
Con_dense = list(Con_dense)
top_cons = max(Con_dense)
Tor_dense = list(Tor_dense)
top_tors = max(Tor_dense)

#Average ribosome density is calculated on a per codon basis, and confidence intervals are also generated.
#Statistical analysis is also performed using a Mann-Whitney test.
(dim_x, dim_y) = Con_Codons.shape

outfile = open('Normalized_RPFs_Percent_ORF.txt', 'w')

for i in range(dim_y):
    cons = Con_Codons[:, i]
    cons = [x/top_cons for x in cons]
    tors = Tor_Codons[:, i]
    tors = [x/top_tors for x in tors]
    n = len(cons)
    con_avg, se_cons = np.mean(cons), stats.sem(cons)
    tors_avg, se_tors = np.mean(tors), stats.sem(tors)
    h_cons = se_cons * stats.t.ppf((1 + 0.99)/2.0, n - 1)
    h_tors = se_tors * stats.t.ppf((1 + 0.99)/2.0, n - 1)
    (U_stat, p_value) = stats.mannwhitneyu(cons, tors)
    summary = [con_avg, tors_avg, h_cons, h_tors, p_value]
    summary = [str(x) for x in summary]

    outfile.write('\t'.join(summary) + '\n')

outfile.close()
