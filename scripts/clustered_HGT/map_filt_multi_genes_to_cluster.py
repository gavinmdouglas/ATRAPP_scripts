#!/usr/bin/python3

import argparse
import os
import sys
import gzip
from collections import defaultdict

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from functions import read_fasta

def main():

    parser = argparse.ArgumentParser(

        description='''
Get mapfile of genes to clusters, based on the following conditions:
1) Genes must be on scaffolds >= 5kbp.
2) Genes must be across genomes from at least two different genera.
3) Cluster alignment must be at least 500 bp.
4) At least two genes must be in the >95% identity cluster, after these above conditions.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-o', '--output',
                        metavar="OUTPUT", type=str,
                        help="Path to output file.",
                        required=True)

    parser.add_argument('-c', '--cdhit_output',
                        metavar="CDHIT-OUTPUT", type=str,
                        help="Path to cdhit output file with gene clusters (gzipped).",
                        required=False,
                        default='/mfs/gdouglas/projects/ATRAPP/cdhit_working/atrapp_genes_cdhit_clusters.txt.clstr.gz')

    parser.add_argument('-a', '--alignment_folder',
                        metavar="ALIGN-FOLDER", type=str,
                        help="Path to folder with gene cluster alignments.",
                        required=False,
                        default='/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/cluster_seqs_aligned')

    parser.add_argument('-t', '--tax',
                        metavar="TAX-TAB", type=str,
                        help="Path to (gzipped) table with taxonomic information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz')

    parser.add_argument('-g', '--geneinfo',
                        metavar="GENE-INFO", type=str,
                        help="Path to (gzipped) table with gene information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ATRAPP/genome_process/gene_info.tsv.gz')

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

    gene_to_genome = dict()
    gene_to_scaffold = dict()
    gene_tab_col = {}
    passing_scaffolds = set()
    with gzip.open(args.geneinfo, "rb") as gene_fh:
        # gene_id gene_type       gene_length     genome_id       scaffold_id     scaffold_length
        gene_header = gene_fh.readline().decode("utf-8").strip().split("\t")
        for i, col in enumerate(gene_header):
            gene_tab_col[col] = i

        for gene_line in gene_fh:
            gene_line = gene_line.decode("utf-8").strip().split("\t")

            gene_to_genome[gene_line[gene_tab_col['gene_id']]] = gene_line[gene_tab_col['genome_id']]
            gene_to_scaffold[gene_line[gene_tab_col['gene_id']]] = gene_line[gene_tab_col['scaffold_id']]

            scaffold_length = int(gene_line[gene_tab_col['scaffold_length']])
            if scaffold_length >= 5000:
                passing_scaffolds.add(gene_line[gene_tab_col['scaffold_id']])

    pass_ids = set()
    cluster_to_genes = {}
    cluster_id = None
    current_cluster_genes = []

    with gzip.open(args.cdhit_output, 'rt') as cluster_fh:
        for line in cluster_fh:

            if line[0] == '>':
                if len(current_cluster_genes) > 1:
                    cluster_to_genes[cluster_id] = current_cluster_genes

                cluster_id = 'c' + line.split()[1]
                current_cluster_genes = []
            else:
                seq_id = line.split()[2][1:].split('...')[0]
                if seq_id in pass_ids:
                    sys.exit('Sequence ' + seq_id + ' appears in multiple clusters.')
                else:
                    pass_ids.add(seq_id)
                current_cluster_genes.append(seq_id)

    # And for final cluster if it has at least two sequences.
    if len(current_cluster_genes) > 1:
        cluster_to_genes[cluster_id] = current_cluster_genes

    with open(args.output, 'w') as out_fh:
        out_fh.write("gene_id\tcluster_id\tgenome_id\tscaffold_id\n")

        for cluster_id, gene_list in cluster_to_genes.items():

            # Filter genes based on scaffold length.
            filtered_genes = [g for g in gene_list if gene_to_scaffold[g] in passing_scaffolds]

            if len(filtered_genes) < 2:
                continue

            # Check genera represented.
            genera = set()
            for gene_id in filtered_genes:
                genome_id = gene_to_genome[gene_id]
                genus = tax[genome_id]["Genus"]
                genera.add(genus)

            if len(genera) < 2:
                continue

            # Check alignment length.
            align_path = os.path.join(args.alignment_folder, cluster_id + '.fa')
            if not os.path.exists(align_path):
                sys.exit(f'Alignment file {align_path} does not exist.')

            seqs = read_fasta(align_path)
            align_length = len(next(iter(seqs.values())))
            if align_length < 500:
                continue

            # If all conditions are met, write to output.
            for gene_id in filtered_genes:
                out_fh.write(f"{gene_id}\t{cluster_id}\t{gene_to_genome[gene_id]}\t{gene_to_scaffold[gene_id]}\n")

if __name__ == '__main__':
    main()
