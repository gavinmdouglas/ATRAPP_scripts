CD-HIT was run on all genes across MAGs, to create sequence clusters. This was run with the options `-c 0.95 -M 100000 -T 64 -g 1  -d 0  -s 0.95  -aL 0.95 -aS 0.95`.

We then post-processed the clusters to identify reciprocal best-hits of genes across genomes

1. Get FASTAs of all clusters with more than 1 sequence.

```
python /mfs/gdouglas/scripts/ATRAPP_scripts/cdhit_postprocess/clusters_to_fastas.py \
	-c /mfs/gdouglas/projects/ocean_mags/clusters/gene-catalog.output.clstr \
	-f /mfs/gdouglas/projects/ocean_mags/Sunagawa_dataset/all_cds_sequences.fasta \
	-o /mfs/gdouglas/projects/ocean_mags/clusters/cluster_seqs
```


# Align with MAFFT.
mkdir /mfs/gdouglas/projects/ocean_mags/clusters/cluster_seqs_aligned

python prep_mafft_cmds.py

conda activate phylo

cat mafft_cmds.sh | parallel -j 200 --joblog mafft_cmds.log '{}'

python ~/local/parallel_joblog_summary/joblog_summary.py --cmds mafft_cmds.sh --log mafft_cmds.log --cmds_to_run mafft_cmds_remain1.sh --failed_cmds mafft_cmds_failed1.sh

# To get full paths for input files.
sed 's/--preservecase cluster_seqs/--preservecase \/mfs\/gdouglas\/projects\/ocean_mags\/clusters\/cluster_seqs/g' mafft_cmds_remain1.sh > mafft_cmds_remain1_fullpath.sh

# To avoid having to type password every time (which would be impossible for all these connections):
eval $(ssh-agent -s)
ssh-add ~/.ssh/id_rsa

# Shared server command:
cat mafft_cmds_remain1_fullpath.sh | parallel --joblog mafft_cmds_remain1_fullpath.log -S 256/sp5000,96/sp2000,30/sp64,30/sp63,15/sp32,15/sp31,15/sp28 "source /mfs/gdouglas/local/miniforge3/etc/profile.d/conda.sh; conda activate phylo; '{}'"

# For the extra quotes around the command caused issues when comparing commands later, so removed them from log file:
sed "s/'//g" mafft_cmds_remain1_fullpath.log > mafft_cmds_remain1_fullpath_noquote.log

# Also had to remove the source part of the command:
sed "s/source \/mfs\/gdouglas\/local\/miniforge3\/etc\/profile.d\/conda.sh; conda activate phylo; //g" mafft_cmds_remain1_fullpath_noquote.log > mafft_cmds_remain1_fullpath_noquote_cmdonly.log

python ~/local/parallel_joblog_summary/joblog_summary.py --cmds mafft_cmds_remain1_fullpath.sh --log mafft_cmds_remain1_fullpath_noquote_cmdonly.log --cmds_to_run mafft_cmds_fullpath_noquote_cmdonly_remain2.sh --failed_cmds mafft_cmds_remain1_fullpath_noquote_cmdonly_failed1.sh


# Then ran final set:
cat mafft_cmds_remain1_fullpath_noquote_cmdonly_failed1.sh | parallel -j 100 --joblog mafft_cmds_remain1_fullpath_noquote_cmdonly_failed1_rerun.log '{}'

python ~/local/parallel_joblog_summary/joblog_summary.py --cmds mafft_cmds_remain1_fullpath_noquote_cmdonly_failed1.sh --log mafft_cmds_remain1_fullpath_noquote_cmdonly_failed1_rerun.log

# Get reciprocal best hits, on passing scaffolds, and between different genera and above.

python /mfs/gdouglas/scripts/ocean_mag_hgt/scripts/clustered_HGT/putative_hgt_per_cluster.py \
	-i ./all_aligned_paths.txt \
	-o all_best_hits.tsv \
	2> all_best_hits.log \
	> /dev/null


# First, get sets of files to process (in 100 random splits, to make it easily to run in parallel).
mkdir aligned_paths
conda deactivate
python /mfs/gdouglas/scripts/ocean_mag_hgt/scripts/clustered_HGT/prep_aligned_paths.py

conda activate phylo
mkdir best_hits_raw
for PATHFILE in /mfs/gdouglas/projects/ocean_mags/clusters/aligned_paths/*.txt; do
	CHUNK=$( basename $PATHFILE .txt )
	echo "python /mfs/gdouglas/scripts/ocean_mag_hgt/scripts/clustered_HGT/putative_hgt_per_cluster.py \
			-i $PATHFILE \
			-o /mfs/gdouglas/projects/ocean_mags/clusters/best_hits_raw/$CHUNK.tsv \
			> /dev/null" >> putative_hgt_per_cluster_cmds.sh
done

cat putative_hgt_per_cluster_cmds.sh | parallel -j 100 --eta --joblog putative_hgt_per_cluster_cmds.log {}

# All output saved to "all_best_hits.tsv.gz"

# Afterwards realized that it would be convenient to have the representative gene in this table too.
# Did so with:

python /mfs/gdouglas/scripts/ocean_mag_hgt/scripts/clustered_HGT/get_besthit_clusters.py \
	> all_best_hits_w_repID.tsv

# Then to get table of tallies of HGT counts per pairwise genome comparison.

python /mfs/gdouglas/scripts/ocean_mag_hgt/scripts/clustered_HGT/pairwise_genome_summary.py > clusterbased_hgt_tallies.tsv
