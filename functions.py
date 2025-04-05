#!/usr/bin/python3

import gzip
import sys
import textwrap

def read_fasta(filename, cut_header=False):
    seq = {}
    name = None

    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:
        line = line.rstrip()
        if len(line) == 0:
            continue
        if line[0] == ">":
            if cut_header:
                name = line.split()[0][1:]
            else:
                name = line[1:]
            name = name.rstrip()
            if name in seq:
                sys.exit("Stopping due to duplicated id in file: " + name)
            seq[name] = ''
        else:
            line = line.rstrip()
            seq[name] += line

    fasta_in.close()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")
    out_fasta.close()
