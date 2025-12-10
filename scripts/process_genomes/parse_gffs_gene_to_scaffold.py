#!/usr/bin/python3

import os
import sys

# Parse GFFs to get mapping of gene IDs to scaffold and genome IDs. Also include the scaffold lengths.

# Folder with Prokka output folders. 
indir = "/mfs/nicot/mybook2/MAGs_ATRAPP/dRep_out_MAGs_last_98/dereplicated_genomes/prok_out"

prokka_folders = []
for item in os.listdir(indir):
    item_path = os.path.join(indir, item)
    if os.path.isdir(item_path):
        prokka_folders.append(item)

scaffold_lengths = {}

print("gene_id\tgene_type\tgene_length\tgenome_id\tscaffold_id\tscaffold_length")

parsed_genomes = set()
for prokka_folder in prokka_folders:
    if prokka_folder == 'gff' or prokka_folder == "ffn_file":
        continue

    full_folder = os.path.join(indir, prokka_folder)
    expected_gff = os.path.join(full_folder, prokka_folder.replace('.prokk', '') + ".gff")
    if not os.path.exists(expected_gff):
        sys.exit(f"Expected GFF file {expected_gff} does not exist.")
    genome_id = prokka_folder.replace('.fa.cleaned.fa', '')
    genome_id = genome_id.replace('.fa.cleaned', '')
    genome_id = genome_id.replace('.cleaned.fa', '')
    genome_id = genome_id.replace('.cl.fa', '')
    genome_id = genome_id.replace('.fa.cl', '')
    genome_id = genome_id.replace('.fa.prokk', '')
    genome_id = genome_id.replace('.prokk', '')
    genome_id = genome_id.replace('.cl', '')

    if '.fa' in genome_id or 'eaned' in genome_id or '.cl' in genome_id or 'prokk' in genome_id:
        sys.exit(f"Unexpected genome ID format: {genome_id}")
    
    if genome_id in parsed_genomes:
        sys.exit(f"Duplicate genome ID found: {genome_id}.")
    parsed_genomes.add(genome_id)

    scaffold_lengths = dict()

    with open(expected_gff, 'r') as gff_file:
        for line in gff_file:
            if line.startswith("#"):
                if line.startswith("##sequence-region"):
                    line_split = line.split()
                    scaffold_id = line_split[1]
                    if line_split[2] != "1":
                        sys.exit(f"Unexpected sequence region format: {line}.")
                    scaffold_lengths[scaffold_id] = line_split[3]
            else:
                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue
                scaffold_id = fields[0]
                feature_type = fields[2]
                if feature_type == 'repeat_region':
                    continue
                feature_info = fields[8]
                feature_id = feature_info.split(";")[0].split("=")[1]
                feature_length = str(int(fields[4]) - int(fields[3]) + 1)
                if scaffold_id not in scaffold_lengths:
                    sys.exit(f"Scaffold ID {scaffold_id} not found in scaffold lengths.")
                print(f"{feature_id}\t{feature_type}\t{feature_length}\t{genome_id}\t{scaffold_id}\t{scaffold_lengths[scaffold_id]}")