# For each genome, determine the total number of genes,
# and the breakdown of genes with HGT hits, that are across multiple species, genera, families, orders, classes, phyla, domains.
# Return the tally of the number of different genes found across different species, strains, etc.
# Also return the number of genes in different MGE categories.

import pandas as pd
from collections import defaultdict
import sys

# Read in pairwise HGT hits.
hgt_hits = pd.read_csv("/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/all_best_hits.tsv.gz", sep="\t", compression="gzip")

# Read in mappings of gene ID to cluster ID and genome ID,
# as well as genome ID to taxa and cluster ID to COG info.
genome_to_taxa = pd.read_csv("/mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz", sep="\t", index_col=0)

gene_info = pd.read_csv("/mfs/gdouglas/projects/ATRAPP/genome_process/gene_info.tsv.gz", sep="\t", index_col=0)

gene_to_cluster = pd.read_csv("/mfs/gdouglas/projects/ocean_mags/clusters/gene_to_cluster.tsv.gz", sep="\t", index_col=0)
scaffolds_to_keep = set(pd.read_csv("/mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/scaffolds_5000bp.txt", header=None)[0].tolist())

# Subset hgt_hits to where gene1_genome and gene2_genome are in genome_to_taxa.index.
hgt_hits = hgt_hits[hgt_hits['gene1_genome'].isin(genome_to_taxa.index) & hgt_hits['gene2_genome'].isin(genome_to_taxa.index)]

# Get breakdown of HGT hits per cluster.
# Loop through each row of hgt_hits.
genome_hgt_counts = defaultdict(int)

for idx, row in hgt_hits.iterrows():
    gene1 = row['gene1']
    gene2 = row['gene2']

    if gene1 not in gene_to_cluster.index:
        sys.exit(f"Gene {gene1} not found in gene_to_cluster.")
    if gene2 not in gene_to_cluster.index:
        sys.exit(f"Gene {gene2} not found in gene_to_cluster.")
    cluster1 = gene_to_cluster.loc[gene1, 'cluster_id']
    cluster2 = gene_to_cluster.loc[gene2, 'cluster_id']

    if cluster1 != cluster2:
        raise ValueError(f"Genes {gene1} and {gene2} are in different clusters: {cluster1} and {cluster2}.")

    genome_hgt_counts[row['gene1_genome']] += 1
    genome_hgt_counts[row['gene2_genome']] += 1

# Get clusters to genes.
cluster_to_genes = defaultdict(set)
for gene, row in gene_to_cluster.iterrows():
    if gene in gene_info_w_mge.index:
        cluster_to_genes[row['cluster_id']].add(gene)

tax_levels = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']
genomes_genes_across_mult_levels = dict()
for level in tax_levels:
    genomes_genes_across_mult_levels[level] = defaultdict(int)

total_genes = defaultdict(int)
genes_on_scaffold_w_provirus = defaultdict(int)
proMGE_genes = defaultdict(int)
genomad_plasmid_genes = defaultdict(int)
genomad_virus_genes = defaultdict(int)

for cluster_id, genes in cluster_to_genes.items():
    genomes_in_cluster = set()

    for gene in genes:
        if not gene in gene_info_w_mge.index:
            sys.exit(f"Gene {gene} not found in gene_info_w_mge.")
        gene_info = gene_info_w_mge.loc[gene]
        genomes_in_cluster.add(gene_info['genome'])
        total_genes[gene_info['genome']] += 1
        if gene_info['scaffold_contains_provirus'] == 'Yes':
            genes_on_scaffold_w_provirus[gene_info['genome']] += 1
        if gene_info['proMGE'] == 'Yes':
            proMGE_genes[gene_info['genome']] += 1
        if gene_info['genomad'] == 'Plasmid':
            genomad_plasmid_genes[gene_info['genome']] += 1
        if gene_info['genomad'] == 'Virus':
            genomad_virus_genes[gene_info['genome']] += 1

    if len(genomes_in_cluster) == 0:
        continue

    genome_taxa_subset = genome_to_taxa.loc[list(genomes_in_cluster)]

    if len(genomes_in_cluster) != len(genome_taxa_subset):
        raise ValueError(f"Some genomes in cluster {cluster_id} are not in genome_to_taxa.")

    for level in tax_levels:
        if genome_taxa_subset[level].nunique() > 1:
            for genome in genomes_in_cluster:
                genomes_genes_across_mult_levels[level][genome] += 1

print("genome\ttotal_genes\thgt_gene_calls\tgenes_on_scaffold_w_provirus\tproMGE_genes\tgenomad_plasmid_genes\tgenomad_virus_genes\tgenes_across_mult_species\tgenes_across_mult_genera\tgenes_across_mult_families\tgenes_across_mult_orders\tgenes_across_mult_classes\tgenes_across_mult_phyla\tgenes_across_mult_domains")

all_genomes = sorted(list(total_genes.keys()))

for genome in all_genomes:
    total = total_genes[genome]
    hgt = genome_hgt_counts[genome]
    provirus = genes_on_scaffold_w_provirus[genome]
    proMGE = proMGE_genes[genome]
    plasmid = genomad_plasmid_genes[genome]
    virus = genomad_virus_genes[genome]
    mult_species = genomes_genes_across_mult_levels['Species'][genome]
    mult_genera = genomes_genes_across_mult_levels['Genus'][genome]
    mult_families = genomes_genes_across_mult_levels['Family'][genome]
    mult_orders = genomes_genes_across_mult_levels['Order'][genome]
    mult_classes = genomes_genes_across_mult_levels['Class'][genome]
    mult_phyla = genomes_genes_across_mult_levels['Phylum'][genome]
    mult_domains = genomes_genes_across_mult_levels['Domain'][genome]

    print(f"{genome}\t{total}\t{hgt}\t{provirus}\t{proMGE}\t{plasmid}\t{virus}\t{mult_species}\t{mult_genera}\t{mult_families}\t{mult_orders}\t{mult_classes}\t{mult_phyla}\t{mult_domains}")
