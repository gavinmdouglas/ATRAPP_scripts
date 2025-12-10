#!/usr/bin/python3

import os
import sys
import gzip
from collections import defaultdict

# Get breakdown of number of HGT calls of each type per pairwise genome comparison.
cluster_hitfile = sys.argv[1]

if not os.path.exists(cluster_hitfile):
    print(f"Error: {cluster_hitfile} does not exist.")
    sys.exit(1)

hits_95 = defaultdict(int)
hits_99 = defaultdict(int)
taxa_diff = {}

with gzip.open(cluster_hitfile, 'rt') as f:
    header = f.readline().strip().split("\t")
    col_to_i = {x:i for i,x in enumerate(header)}
    for line in f:
        line = line.strip().split("\t")
        genome1 = line[col_to_i["gene1_genome"]]
        genome2 = line[col_to_i["gene2_genome"]]
        genome_pair = ','.join(sorted([genome1, genome2]))
        identity = float(line[col_to_i["identity"]])
        highest_tax_diff = line[col_to_i["highest_tax_diff"]]

        if genome_pair in taxa_diff:
            if highest_tax_diff != taxa_diff[genome_pair]:
                print("Taxa diff mismatch")
                sys.exit()
        else:
            taxa_diff[genome_pair] = highest_tax_diff
    
        if identity >= 95 and identity < 99:
            hits_95[genome_pair] += 1
        elif identity >= 99:
            hits_99[genome_pair] += 1
        else:
            print("Identity not in range")
            sys.exit()

print("\t".join(["Taxa_combo", "Highest_tax_diff", "95_gene_count", "99_gene_count", "both_gene_count"]))
for genome_pair in sorted(list(taxa_diff.keys())):
    print("\t".join([genome_pair, taxa_diff[genome_pair], str(hits_95[genome_pair]), str(hits_99[genome_pair]), str(hits_95[genome_pair] + hits_99[genome_pair])]))
