#!/usr/bin/python3

import argparse
import sys
import gzip
from collections import defaultdict


def main():

    parser = argparse.ArgumentParser(

        description='''
Compute prevalence across taxa per each gene cluster found in at least two different genera, with genes on scaffolds >=5kbp and alignments >=500bp.
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

    print('\t'.join(['Cluster_ID', 'Num_Genomes',
                     'Num_Genus', 'Num_Genus_MultiGenome',
                     'Num_Family', 'Num_Family_MultiGenome',
                     'Num_Order', 'Num_Order_MultiGenome',
                     'Num_Class', 'Num_Class_MultiGenome',
                     'Num_Phylum', 'Num_Phylum_MultiGenome',
                     'Num_Domain', 'Num_Domain_MultiGenome']))

    for cluster_id, genomes in cluster_to_genomes.items():
        outline = [cluster_id, str(len(genomes))]
        
        higher_tax_levels = ['Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']

        for tax_level in higher_tax_levels:
            taxa_counts = defaultdict(int)
            for genome_id in genomes:
                if genome_id not in tax:
                    sys.exit(f"Genome {genome_id} not found in taxonomy.")
                taxa = tax[genome_id][tax_level]
                taxa_counts[taxa] += 1
            
            num_taxa = len(taxa_counts)
            num_multi_genome_taxa = sum(1 for count in taxa_counts.values() if count > 1)

            outline.append(str(num_taxa))
            outline.append(str(num_multi_genome_taxa))
        print('\t'.join(outline))


if __name__ == '__main__':
    main()
