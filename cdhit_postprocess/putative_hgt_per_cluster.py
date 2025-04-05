#!/usr/bin/python3

import argparse
import os
import sys
import gzip
from collections import defaultdict
import hammingdist

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

def main():

    parser = argparse.ArgumentParser(

        description='''
Parse aligned FASTAs of sequences from the same cluster.
Identify reciprocal best-hits above 95% (that are at least between genomes of different genera or above).
Also, consider all hits initially, but only output reciprocal best hits that are on passing scaffolds.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="FILES", type=str,
                        help="Path to file with paths to input FASTAs. One per line. Done this way to make it trivial to parallelize.",
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar="OUTPUT", type=str,
                        help="Path to output file.",
                        required=True)

    parser.add_argument('-t', '--tax',
                        metavar="TAXTAB", type=str,
                        help="Path to (gzipped) table with taxonomic information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/mapfiles/MAG_taxa_breakdown.tsv.gz')

    parser.add_argument('-p', '--passing',
                        metavar="PASS", type=str,
                        help="Path to file with passing scaffolds.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt')

    parser.add_argument('-g', '--geneinfo',
                        metavar="GENEINFO", type=str,
                        help="Path to (gzipped) table with gene information per genome.",
                        required=False,
                        default='/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/gene_info_allscaffolds.tsv.gz')

    args = parser.parse_args()

    passing_scaffolds = set()
    with open(args.passing, 'r') as pass_fh:
        for line in pass_fh:
            passing_scaffolds.add(line.strip())

    # Read in genome taxonomy.
    tax = defaultdict(dict)
    tax_levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Strain"]
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
    with gzip.open(args.geneinfo, "rb") as gene_fh:
        gene_header = gene_fh.readline().decode("utf-8").strip().split("\t")
        for i, col in enumerate(gene_header):
            gene_tab_col[col] = i

        for gene_line in gene_fh:
            gene_line = gene_line.decode("utf-8").strip().split("\t")
            gene_to_genome[gene_line[gene_tab_col['gene']]] = gene_line[gene_tab_col['genome']]
            gene_to_scaffold[gene_line[gene_tab_col['gene']]] = gene_line[gene_tab_col['scaffold']]

    with open(args.output, 'w') as out_fh:

        print("\t".join(['gene1', 'gene2', 'gene1_genome', 'gene2_genome', 'highest_tax_diff', 'identity']), file=out_fh)

        with open(args.input, 'r') as input_fh:
            for filepath in input_fh:
                filepath = filepath.strip()
                if not filepath or len(filepath) == 0:
                    continue
            
                seq_ids = []
                seqs = {}
                name = None
                with open(filepath, 'r') as fasta_in:
                    for line in fasta_in:
                        line = line.rstrip()
                        if len(line) == 0:
                            continue
                        # If header-line then split by whitespace, take the first element,
                        # and define the sequence name as everything after the ">".
                        if line[0] == ">":
                            name = line[1:]
                            if name in seqs:
                                sys.stderr("Stopping due to duplicated id in file: " + name)
                            seqs[name] = ''
                            seq_ids.append(name)
                        else:
                            seqs[name] += line

                # Make sure all sequences are same length.
                seq_len = len(seqs[list(seqs.keys())[0]])
                for seq in seqs.values():
                    if len(seq) != seq_len:
                        sys.exit("Error: Not all sequences are the same length.")

                # Skip if seq_len <= 500.
                if seq_len <= 500:
                    continue

                # Get pairwise distances of aligned sequences.
                dist = hammingdist.from_fasta(filepath)                    

                top_matches = {}
                top_identities = {}
                ambig_seqs = set()
                for i in range(len(seq_ids) - 1):
                    for j in range(i + 1, len(seq_ids)):
                        i_id = seq_ids[i]
                        j_id = seq_ids[j]
                        iden = (1 - (dist[i, j] / seq_len)) * 100.0
                        if iden >= 95.0:
                            if i_id not in top_matches or iden > top_identities[i_id]:
                                top_matches[i_id] = j_id
                                top_identities[i_id] = iden
                            elif iden == top_identities[i_id]:
                                ambig_seqs.add(i_id)

                            if j_id not in top_matches or iden > top_identities[j_id]:
                                top_matches[j_id] = i_id
                                top_identities[j_id] = iden
                            elif iden == top_identities[j_id]:
                                ambig_seqs.add(j_id)
                
                parsed_seqs = set()
                for seq_id in seq_ids:
                    if seq_id not in top_matches or seq_id in ambig_seqs or seq_id in parsed_seqs:
                        continue
                        
                    best_match = top_matches[seq_id]
                    reverse_best_match = top_matches[best_match]

                    if seq_id == reverse_best_match and best_match not in ambig_seqs and best_match not in parsed_seqs:
                        parsed_seqs.add(best_match)
                        parsed_seqs.add(seq_id)

                        if seq_id not in gene_to_scaffold or best_match not in gene_to_scaffold:
                            continue

                        gene1_scaffold = gene_to_scaffold[seq_id]
                        gene2_scaffold = gene_to_scaffold[best_match]
                        if gene1_scaffold not in passing_scaffolds or gene2_scaffold not in passing_scaffolds:
                            continue

                        gene1_genome = gene_to_genome[seq_id]
                        gene2_genome = gene_to_genome[best_match]
                        if gene1_genome == gene2_genome:
                            continue

                        diff_found = False
                        for tax_level in tax_levels:
                            if tax[gene1_genome][tax_level] != tax[gene2_genome][tax_level]:
                                highest_tax_diff = tax_level
                                diff_found = True
                                break
                        if not diff_found:
                            sys.exit("Error: No taxonomic difference found between genomes " + gene1_genome + " and " + gene2_genome + ".")

                        genome_order = sorted([])

                        if highest_tax_diff in  ['Species', 'Strain']:
                            continue
                        
                        print("\t".join([seq_id, best_match, gene1_genome, gene2_genome, highest_tax_diff, str(top_identities[seq_id])]), file=out_fh)

if __name__ == '__main__':
    main()
