#!/usr/bin/python3

import argparse
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from functions import write_fasta, read_fasta


def main():

    parser = argparse.ArgumentParser(

        description='''
Parse CD-HIT output file of cluster definitions, and original input FASTA. Output FASTA files for each cluster with at least two sequences.''',
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-c', '--cluster',
                        metavar="CLUSTER", type=str,
                        help="Path to CD-HIT cluster output file.",
                        required=True)

    parser.add_argument('-f', '--fasta',
                        metavar="FASTA", type=str,
                        help="Path to input FASTA file for clustering.",
                        required=True)

    parser.add_argument('-o', '--outdir',
                        metavar="OUTDIR", type=str,
                        help="Path to output directory.",
                        required=True)
    
    parser.add_argument('-t', '--trim',
                        action='store_true',
                        help="Trim header of sequences in FASTA to first whitespace-delimited field.")

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header=args.trim)

    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir)

    # Keep track of previously parsed sequence IDs as a sanity check.
    pass_ids = set()
    cluster_seqs = {}
    cluster_id = None

    with open(args.cluster, 'r') as cluster_fh:
        for line in cluster_fh:

            if line[0] == '>':
                if len(cluster_seqs) > 1:
                    write_fasta(cluster_seqs, os.path.join(args.outdir, cluster_id + '.fa'))

                cluster_id = 'c' + line.split()[1]
                cluster_seqs = {}
            else:
                seq_id = line.split()[2][1:].split('...')[0]
                if seq_id in pass_ids:
                    sys.exit('Sequence ' + seq_id + ' appears in multiple clusters.')
                else:
                    pass_ids.add(seq_id)
                cluster_seqs[seq_id] = seqs[seq_id]

    # And for final cluster if it has at least two sequences.
    if len(cluster_seqs) > 1:
        write_fasta(cluster_seqs, os.path.join(args.outdir, cluster_id + '.fa'))


if __name__ == '__main__':
    main()
