## 1. Get mapfile of genes to clusters after filtering steps.

These steps include:
- Genes must be on scaffolds >= 5000 bp.
- Genes must be across genomes from at least two different genera.
- Cluster alignment must be >= 500 bp.
- At least two genes must be in the >95% identity cluster, after these above conditions.

That is done with the following script:
```
python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/clustered_HGT/map_filt_multi_genes_to_cluster.py \
	-o /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/gene_to_cluster_multi_filt.tsv

gzip /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/gene_to_cluster_multi_filt.tsv
```

## 2. Get number of higher-order taxa that encode each cluster

Broken down by the total number of genera, families, etc., as well as the numbers of cases where higher order taxa have at least two genomes with the cluster (to give a sense of more widely spread clusters).

```
cd /mfs/gdouglas/projects/ATRAPP/HGT_prevalence/

python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/prevalence_HGT/compute_taxon_prevalence_per_cluster.py \
	-i /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/gene_to_cluster_multi_filt.tsv.gz \
	--tax /mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz \
	> cluster_taxon_prevalence.tsv

gzip cluster_taxon_prevalence.tsv
```

## 3. Get number of filtered clusters shared between pairwise genomes.

The motivation for this data type is that it could be used to make a quick network. I also included the phylum level of each genome, to make it easier to see which are Cyanobacteria vs. others.

```
python /mfs/gdouglas/scripts/ATRAPP_scripts/scripts/prevalence_HGT/compute_pairwise_clusters_shared.py \
	-i /mfs/gdouglas/projects/ATRAPP/cdhit_postprocess/gene_to_cluster_multi_filt.tsv.gz \
	--tax /mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz \
	> pairwise_genomes_clusters_shared.tsv

gzip pairwise_genomes_clusters_shared.tsv
```

