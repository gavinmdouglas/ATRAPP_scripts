#!/usr/bin/python3

import argparse
import gzip
from collections import defaultdict


def main():

    parser = argparse.ArgumentParser(

        description='''
Compute number of filtered clusters shared by each pairwise genome, and include the phylum level of each genome.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="INPUT", type=str,
                        help="Path to (gzipped) table with filtered genes to clusters (from map_filt_multi_genes_to_cluster.py).",
                        required=True)

    parser.add_argument('-t', '--tax',
                        metavar="TAX-TAB", type=str,
                        help="Path to (gzipped) table with taxonomic information per genome.",
                        required=True)

    args = parser.parse_args()

    # Read in genome taxonomy.
    tax = defaultdict(dict)
    tax_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "MAG"]
    tax_col_to_i = {}
    with gzip.open(args.tax, "rb") as tax_fh:
        tax_header = tax_fh.readline().decode("utf-8").strip().split("\t")
        for i, col in enumerate(tax_header):
            tax_col_to_i[col] = i

        for tax_line in tax_fh:
            tax_line = tax_line.decode("utf-8").strip().split("\t")
            mag_id = tax_line[tax_col_to_i["MAG"]]
            for tax_level in tax_levels:
                tax[mag_id][tax_level] = tax_line[tax_col_to_i[tax_level]]

    # Read in cluster to genomes.
    cluster_to_genomes = defaultdict(set)
    with gzip.open(args.input, "rb") as input_fh:
        input_header = input_fh.readline().decode("utf-8").strip().split("\t")
        input_col_to_i = {col: i for i, col in enumerate(input_header)}

        for input_line in input_fh:
            input_line = input_line.decode("utf-8").strip().split("\t")
            cluster_id = input_line[input_col_to_i["cluster_id"]]
            genome_id = input_line[input_col_to_i["genome_id"]]

            cluster_to_genomes[cluster_id].add(genome_id)

    print('\t'.join(["genome1", "genome2", "num_shared_clusters", "genome1_phylum", "genome2_phylum"]))

    genome_pairs = defaultdict(int)
    for cluster_id, genomes in cluster_to_genomes.items():
        genomes = sorted(list(genomes))
        for i in range(len(genomes)):
            for j in range(i + 1, len(genomes)):
                genome1 = genomes[i]
                genome2 = genomes[j]
                genome_pairs[(genome1, genome2)] += 1

    for (genome1, genome2), shared_count in genome_pairs.items():
        genome1_phylum = tax[genome1]["Phylum"] if genome1 in tax else "Unknown"
        genome2_phylum = tax[genome2]["Phylum"] if genome2 in tax else "Unknown"
        print('\t'.join([genome1, genome2, str(shared_count), genome1_phylum, genome2_phylum]))


if __name__ == '__main__':
    main()
