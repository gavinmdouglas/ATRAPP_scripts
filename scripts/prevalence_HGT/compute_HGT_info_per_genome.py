#!/usr/bin/python3

import argparse
import sys
import gzip
from collections import defaultdict
import pandas as pd


def main():

    parser = argparse.ArgumentParser(

        description='''
Compute breakdown of number of clusters involved in HGT, found across other taxa, etc., per genome.
''',

formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input',
                        metavar="INPUT", type=str,
                        help="Path to (gzipped) table with filtered genes to clusters (from map_filt_multi_genes_to_cluster.py).",
                        required=True)

    parser.add_argument('-p', '--prevalence',
                        metavar="PREVALENCE-TAB", type=str,
                        help="Path to (gzipped) table with cluster taxon prevalence information.",
                        required=True)

    parser.add_argument('-b', '--besthits',
                        metavar="BESTHITS-TAB", type=str,
                        help="Path to (gzipped) table with best hits for clusters.",
                        required=True)

    parser.add_argument('-o', '--output',
                        metavar="OUTPUT", type=str,
                        help="Path to output table with per-genome HGT info.",
                        required=True)

    args = parser.parse_args()

    # Read in cluster to genomes.
    genome_to_clusters = defaultdict(set)
    cluster_to_genomes = defaultdict(set)
    gene_to_cluster = dict()
    with gzip.open(args.input, "rb") as input_fh:
        input_header = input_fh.readline().decode("utf-8").strip().split("\t")
        input_col_to_i = {col: i for i, col in enumerate(input_header)}

        for input_line in input_fh:
            input_line = input_line.decode("utf-8").strip().split("\t")
            gene_id = input_line[input_col_to_i["gene_id"]]
            cluster_id = input_line[input_col_to_i["cluster_id"]]
            genome_id = input_line[input_col_to_i["genome_id"]]

            cluster_to_genomes[cluster_id].add(genome_id)
            genome_to_clusters[genome_id].add(cluster_id)
            gene_to_cluster[gene_id] = cluster_id

    # Read in prevalence info.
    prev_tab = pd.read_csv(args.prevalence, sep="\t", compression='gzip', index_col=0)
    
    # Identify columns containing "Multi".
    col2rm = ["Num_Genomes"]
    for col in prev_tab.columns:
        if "Multi" in col:
            col2rm.append(col)
    prev_tab = prev_tab.drop(columns=col2rm)

    # Convert prev_tab to 1 for all values > 1.
    # Do this so that the values indicate there being multiple taxa with the cluster per column.
    prev_tab = (prev_tab > 1).astype(int)
    
    # Replace "Num_" in all column names with "Cross_"
    prev_tab.columns = [col.replace("Num_", "Cross_") for col in prev_tab.columns]

    # Loop through genomes, get subset of rows in prev_tab corresponding to clusters in that genome.
    # Then sum the columns.
    genome_prev_df = pd.DataFrame(index=sorted(list(genome_to_clusters.keys())),
                                  columns=prev_tab.columns,
                                  dtype=int)

    for genome_id, cluster_ids in genome_to_clusters.items():
        genome_clusters = prev_tab.loc[prev_tab.index.isin(cluster_ids), :]
        genome_prev_df.loc[genome_id] = genome_clusters.sum(axis=0).astype(int)

    genome_prev_df.insert(0, "HGT_num", 0)

    # Read in HGT best-hits.
    with gzip.open(args.besthits, "rb") as besthits_fh:
        besthits_header = besthits_fh.readline().decode("utf-8").strip().split("\t")
        besthits_col_to_i = {col: i for i, col in enumerate(besthits_header)}

        for line in besthits_fh:
            line = line.decode("utf-8").strip().split("\t")
            gene1 = line[besthits_col_to_i["gene1"]]
            gene2 = line[besthits_col_to_i["gene2"]]
            cluster1 = gene_to_cluster[gene1]
            cluster2 = gene_to_cluster[gene2]

            if cluster1 != cluster2:
                sys.exit(f"Gene {gene1} and {gene2} map to different clusters: {cluster1}, {cluster2}")

            cluster_id = cluster1

            # Make sure cluster is in index of prev_tab.
            if cluster_id not in prev_tab.index:
                sys.exit(f"Cluster {cluster_id}, which is in best-hits table, not found in prevalence table.")

            # Make sure both genomes are in genome_prev_df.
            genome1 = line[besthits_col_to_i["gene1_genome"]]
            genome2 = line[besthits_col_to_i["gene2_genome"]]

            # Throw error if genomes are the same.
            if genome1 == genome2:
                sys.exit(f"Genes {gene1} and {gene2} are from the same genome {genome1} in best-hits table.")

            # Make sure genomes are in cluster_to_genomes for this cluster.
            if genome1 not in cluster_to_genomes[cluster_id]:
                sys.exit(f"Genome {genome1} not found in cluster {cluster_id}, which is in best-hits table.")
            if genome2 not in cluster_to_genomes[cluster_id]:
                sys.exit(f"Genome {genome2} not found in cluster {cluster_id}, which is in best-hits table.")

            if genome1 not in genome_prev_df.index:
                sys.exit(f"Genome {genome1}, which is in best-hits table, not found in genome prevalence dataframe.")
            if genome2 not in genome_prev_df.index:
                sys.exit(f"Genome {genome2}, which is in best-hits table, not found in genome prevalence dataframe.")

            genome_prev_df.at[genome1, "HGT_num"] += 1
            genome_prev_df.at[genome2, "HGT_num"] += 1

    # Write out.
    genome_prev_df = genome_prev_df.astype(int)
    genome_prev_df.to_csv(args.output, sep="\t", index=True, header=True, index_label="genome_id")


if __name__ == '__main__':
    main()
