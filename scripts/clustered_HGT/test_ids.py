#!/usr/bin/python3

import argparse
import os
import sys
import gzip
from pprint import pprint

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def main():

    parser = argparse.ArgumentParser(

        description='''
Parse aligned FASTAs of sequences from the same cluster.
Identify reciprocal best-hits above 95% (that are at least between genomes of different genera or above).
Also, consider all hits initially, but only output reciprocal best hits that are on scaffolds >= 5000 bp.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="FILES", type=str,
                        help="Path to file with paths to input FASTAs. One per line. Done this way to make it trivial to parallelize.",
                        required=True)

    parser.add_argument('-g', '--geneinfo',
                        metavar="GENEINFO", type=str,
                        help="Path to (gzipped) table with gene information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ATRAPP/genome_process/gene_info.tsv.gz')

    args = parser.parse_args()

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

    with open(args.input, 'r') as input_fh:
        for filepath in input_fh:
            filepath = filepath.strip()
            if not filepath or len(filepath) == 0:
                continue
        
            seq_ids = []
            with open(filepath, 'r') as fasta_in:
                for line in fasta_in:
                    line = line.rstrip()
                    if len(line) == 0:
                        continue
                    # If header-line then split by whitespace, take the first element,
                    # and define the sequence name as everything after the ">".
                    if line[0] == ">":
                        name = line[1:]
                        seq_ids.append(name)

            matched_id = 0
            for raw_seq_id in seq_ids:
                seq_id_split = raw_seq_id.split("_")

                # Remove any elements of seq_id_split that matched ['fa', '.fa', 'cleaned', '.cleaned', 'cl', '.cl', 'prokk', '.prokk']
                seq_id_split = [x for x in seq_id_split if x not in ['fa', '.fa', 'cleaned', '.cleaned', 'cl', '.cl', 'prokk', '.prokk']]

                clean_id = seq_id_split[2] + '_' + seq_id_split[3]
                if clean_id not in gene_to_genome:
                    sys.stderr.write(f"Warning: {clean_id}, (raw: {raw_seq_id}) not found in gene info table.\n")
                else:
                    matched_id += 1

            print(f"Matched {matched_id} out of {len(seq_ids)} sequence IDs to gene info table.")

if __name__ == '__main__':
    main()
