#!/usr/bin/python3

import gzip
import sys

# Get cluster IDs of each pairwise best hits in table produced by putative_hgt_per_cluster.py.
# Also get taxa pair ID, to help speed up analyses.
taxa_file = '/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz'
genome_to_taxon = {}
with gzip.open(taxa_file, 'rt') as taxa_fh:
    taxa_header = taxa_fh.readline()
    for line in taxa_fh:
        line = line.strip().split('\t')
        genome_to_taxon[line[0]] = line[1]

representative_file = '/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene-catalog-membership.tsv.gz'
gene_to_rep = {}
with gzip.open(representative_file, 'rt') as rep_fh:
    rep_header = rep_fh.readline()
    for line in rep_fh:
        gene, rep = line.strip().split('\t')
        gene_to_rep[gene] = rep

with gzip.open('/mfs/gdouglas/projects/ocean_mags/clusters/all_best_hits.tsv.gz', 'rt') as clusters_fh:
    clusters_header = clusters_fh.readline().strip()
    print(clusters_header + '\tsorted_genome_pair\tgene1_rep\tgene2_rep')
    for line in clusters_fh:
        outline = line.strip().split('\t')

        # Sort genome IDs and combine.
        gene1_taxid = genome_to_taxon[outline[2]]
        gene2_taxid = genome_to_taxon[outline[3]]
        sorted_genome_pair = ','.join(sorted([gene1_taxid, gene2_taxid]))
        outline.append(sorted_genome_pair)

        gene1 = outline[0]
        if gene1 in gene_to_rep:
            outline.append(gene_to_rep[gene1])
        else:
            outline.append(gene_to_rep['NA'])
            print('No rep found for gene1: ' + gene2, file=sys.stderr)

        gene2 = outline[1]
        if gene2 in gene_to_rep:
            outline.append(gene_to_rep[gene2])
        else:
            outline.append(gene_to_rep['NA'])
            print('No rep found for gene2: ' + gene2, file=sys.stderr)

        if outline[-1] != outline[-2]:
            print('Rep genes differ for: ' + gene1 + ' and ' + gene2, file=sys.stderr)

        print('\t'.join(outline))
