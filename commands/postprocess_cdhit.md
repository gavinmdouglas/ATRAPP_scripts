0. CD-HIT was run on all genes across MAGs, to create sequence clusters

This was done by first concatenating all gene sequence FASTAs:
cd /mfs/gdouglas/projects/ATRAPP/cdhit_working

```
for genefasta in /mfs/nicot/mybook2/MAGs_ATRAPP/dRep_out_MAGs_last_98/dereplicated_genomes/prok_out/ffn_file/*.ffn; do
	awk '{ print $1 }' >> atrapp_genes.ffn
done
```

And then running cd-hit:

```
conda activate cd-hit

cd-hit-est -i ./atrapp_genes.ffn \
		   -o ./atrapp_genes_cdhit_clusters.txt \
           -c 0.95 \
           -M 100000 \
           -T 64 \
           -g 1 \
           -d 0 \
           -s 0.95 \
           -aL 0.95 \
           -aS 0.95 \
           2> cdhit_stderr.txt \
           > cdhit_stdout.txt

# Remove FASTA of representative sequences, to save space:
rm atrapp_genes_cdhit_clusters.txt

conda deactivate
```

We then **post-processed the clusters to identify reciprocal best-hits of genes across genomes**, using the following steps. Note that several of these scripts are very minor modifications of the same ones
we used for the same workflow for this manuscript: https://github.com/gavinmdouglas/ocean_cooccur_hgt (under GPL 3 license).

1. Get FASTAs of all clusters with more than 1 sequence.
```
cd /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess

python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/cdhit_postprocess/clusters_to_fastas.py \
	-c /mfs/gdouglas/projects/ATRAPP/cdhit_working/atrapp_genes_cdhit_clusters.txt.clstr \
	-f /mfs/gdouglas/projects/ATRAPP/cdhit_working/atrapp_genes.ffn \
	-o /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/cluster_seqs

# Then gzipped large CD-HIT output files.
gzip /mfs/gdouglas/projects/ATRAPP/cdhit_working/*
```


2. Align with MAFFT.

Run MAFFT with GNU parallel. Generated commands use a quick Python script.
```
mkdir /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/cluster_seqs_aligned

# (Note that this is just a quick script and requires hard-coded variable names to be changed at the top)
python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/cdhit_postprocess/prep_mafft_cmds.py
```

Then to run MAFFT (in this case which was installed in the "phylo" conda environment).
```
conda activate phylo
cat mafft_cmds.sh | parallel -j 200 --joblog mafft_cmds.log '{}'
```

And then afterwards ran this command to make sure that all commands finished:
```
python ~/local/parallel_joblog_summary/joblog_summary.py --cmds mafft_cmds.sh --log mafft_cmds.log
```

3. Get mapping of gene IDs to genomes and scaffolds.

Also capture whether it is a CDS or other element type (tRNA), element length, and scaffold lengths.

```
cd /mfs/gdouglas/projects/ATRAPP/genome_process

python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/process_genomes/parse_gffs_gene_to_scaffold.py > gene_info.tsv

gzip gene_info.tsv
```

Genes on scaffolds < 5000 bp will be included in the analysis, but not retained for the final set of best hits, since we are concerned they are more likely to be contaminants.


4. Get reciprocal best hits, on passing scaffolds, and between different genera and above.

This will be with the `/mfs/gdouglas/scripts/ocean_cooccur_hgt/scripts/clustered_HGT/putative_hgt_per_cluster.py` script.

First, get sets of files to process (in 100 random splits, to make it easily to run in parallel), which is done by this hard-coded script.

```
cd /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess

conda deactivate
python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/clustered_HGT/prep_aligned_paths.py
```

Then, to prep commands:
```
conda activate phylo
mkdir best_hits_raw
for PATHFILE in aligned_paths/*.txt; do
	CHUNK=$( basename $PATHFILE .txt )
	echo "python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/clustered_HGT/putative_hgt_per_cluster.py \
			-i $PATHFILE \
			-o best_hits_raw/$CHUNK.tsv \
			> /dev/null" >> putative_hgt_per_cluster_cmds.sh
done
```

Then to run all these commands in parallel:
```
cat putative_hgt_per_cluster_cmds.sh | parallel -j 10 --eta --joblog putative_hgt_per_cluster_cmds.log {}
```

All output saved to "all_best_hits.tsv.gz"

```
cat best_hits_raw/*tsv >> tmp1
head -n 1 tmp1  > all_best_hits.tsv
grep -v "highest_tax_diff" tmp1 >> all_best_hits.tsv
gzip all_best_hits.tsv
rm tmp1
gzip all_best_hits.tsv
```

Get summary of HGT counts by genome pair:
```
python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/clustered_HGT/pairwise_genome_summary.py \
	/mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/all_best_hits.tsv.gz \
	> clusterbased_hgt_tallies.tsv

gzip clusterbased_hgt_tallies.tsv
```
