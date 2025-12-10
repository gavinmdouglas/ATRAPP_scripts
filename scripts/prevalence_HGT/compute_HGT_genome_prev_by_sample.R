rm(list = ls(all.names = TRUE))

# Create tables breaking down proportion of genomes per sample with at least 1 HGT event
# (or one genome with a gene found across multiple genera, etc.)
# Do so for:
# (1) Tara oceans samples matched across lower and upper size fractions (i.e. those from same biological sample, but filtered differently).
# (2) Also do so for Tara oceans samples (for lower fraction, as this has much higher sample size) for which there is matched environmental data.

compute_summary_vals <- function(intab, sampletype) {
  total_genomes <- nrow(intab)
  mean_genes <- mean(intab$total_genes)
  prop_hgt <- length(which(intab$hgt_gene_calls > 0)) / total_genomes
  prop_provirus <- length(which(intab$genes_on_scaffold_w_provirus > 0)) / total_genomes
  prop_proMGE <- length(which(intab$proMGE_genes > 0)) / total_genomes
  prop_plasmid <- length(which(intab$genomad_plasmid_genes > 0)) / total_genomes
  prop_virus <- length(which(intab$genomad_virus_genes > 0)) / total_genomes
  prop_mult_species <- length(which(intab$genes_across_mult_species > 0)) / total_genomes
  prop_mult_genera <- length(which(intab$genes_across_mult_genera > 0)) / total_genomes
  prop_mult_families <- length(which(intab$genes_across_mult_families > 0)) / total_genomes
  prop_mult_orders <- length(which(intab$genes_across_mult_orders > 0)) / total_genomes
  prop_mult_classes <- length(which(intab$genes_across_mult_classes > 0)) / total_genomes
  prop_mult_phyla <- length(which(intab$genes_across_mult_phyla > 0)) / total_genomes
  prop_mult_domains <- length(which(intab$genes_across_mult_domains > 0)) / total_genomes

  col_names <- paste0(sampletype, "_", c("total_genomes", "mean_genes", "prop_hgt", "prop_provirus",
                                         "prop_proMGE", "prop_plasmid", "prop_virus", "prop_mult_species",
                                         "prop_mult_genera", "prop_mult_families", "prop_mult_orders",
                                         "prop_mult_classes", "prop_mult_phyla", "prop_mult_domains"))

  result <- data.frame(total_genomes, mean_genes, prop_hgt, prop_provirus, prop_proMGE, prop_plasmid,
                       prop_virus, prop_mult_species, prop_mult_genera, prop_mult_families,
                       prop_mult_orders, prop_mult_classes, prop_mult_phyla, prop_mult_domains)

  names(result) <- col_names
  return(result)
}

info <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/genome_gene_prevalence.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# So for the sample-matched upper and lower samples.
prev <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_tara_w_protist_fraction.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

env_data <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGEA_env_data.tsv.gz",
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

genome_tax <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                         header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)

prev <- prev[rownames(env_data), ]
prev <- prev[rowSums(prev) >= 10, colSums(prev) > 0]

colnames(prev) <- genome_tax[colnames(prev), "MAG"]

sample_filt_groups <- env_data$sample_w_diff_fraction[duplicated(env_data$sample_w_diff_fraction)]

raw_combined_summary <- list()

for (filt_group in sample_filt_groups) {
  metadata_subset <- env_data[env_data$sample_w_diff_fraction == filt_group, ]

  lower_sample <- rownames(metadata_subset)[which(metadata_subset$Fraction_lower_µm_used_on_board_to_prepare_samp == 0.22 & metadata_subset$Fraction_upper_µm_used_on_board_to_prepare_samp <= 3)]
  upper_sample <- rownames(metadata_subset)[which(metadata_subset$Fraction_lower_µm_used_on_board_to_prepare_samp >= 3 & metadata_subset$Fraction_upper_µm_used_on_board_to_prepare_samp == 20)]

  if ((length(lower_sample) != 1) || (length(upper_sample) != 1)) {
    stop('Not one matching lower and upper each!')
  }

  lower_genomes <- colnames(prev)[which(as.integer(prev[lower_sample, ]) > 0)]
  lower_summary <- compute_summary_vals(info[lower_genomes, ], "lower")

  upper_genomes <- colnames(prev)[which(as.integer(prev[upper_sample, ]) > 0)]
  upper_summary <- compute_summary_vals(info[upper_genomes, ], "upper")

  raw_combined_summary[[filt_group]] <- data.frame(filt_group = filt_group, lower_sample = lower_sample, upper_sample = upper_sample,
                                                   cbind(lower_summary, upper_summary))

}

combined_summary <- do.call(rbind, raw_combined_summary)

write.table(x = combined_summary,
            file = gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_fraction_sample_matched_HGT_prev.tsv.gz"),
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)



# Then for #2 -- Tara samples with more detailed environmental profiles. Split by whether the size fraction was the lower or upper fraction.
rm(list = ls(all.names = TRUE))

compute_summary_vals_simple <- function(intab) {
  total_genomes <- nrow(intab)
  mean_genes <- mean(intab$total_genes)
  prop_hgt <- length(which(intab$hgt_gene_calls > 0)) / total_genomes
  prop_provirus <- length(which(intab$genes_on_scaffold_w_provirus > 0)) / total_genomes
  prop_proMGE <- length(which(intab$proMGE_genes > 0)) / total_genomes
  prop_plasmid <- length(which(intab$genomad_plasmid_genes > 0)) / total_genomes
  prop_virus <- length(which(intab$genomad_virus_genes > 0)) / total_genomes
  prop_mult_species <- length(which(intab$genes_across_mult_species > 0)) / total_genomes
  prop_mult_genera <- length(which(intab$genes_across_mult_genera > 0)) / total_genomes
  prop_mult_families <- length(which(intab$genes_across_mult_families > 0)) / total_genomes
  prop_mult_orders <- length(which(intab$genes_across_mult_orders > 0)) / total_genomes
  prop_mult_classes <- length(which(intab$genes_across_mult_classes > 0)) / total_genomes
  prop_mult_phyla <- length(which(intab$genes_across_mult_phyla > 0)) / total_genomes
  prop_mult_domains <- length(which(intab$genes_across_mult_domains > 0)) / total_genomes

  col_names <- c("total_genomes", "mean_genes", "prop_hgt", "prop_provirus",
                 "prop_proMGE", "prop_plasmid", "prop_virus", "prop_mult_species",
                 "prop_mult_genera", "prop_mult_families", "prop_mult_orders",
                 "prop_mult_classes", "prop_mult_phyla", "prop_mult_domains")

  result <- data.frame(total_genomes, mean_genes, prop_hgt, prop_provirus, prop_proMGE, prop_plasmid,
                       prop_virus, prop_mult_species, prop_mult_genera, prop_mult_families,
                       prop_mult_orders, prop_mult_classes, prop_mult_phyla, prop_mult_domains)

  names(result) <- col_names
  return(result)
}

info <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/genome_gene_prevalence.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# So for the sample-matched upper and lower samples.
prev <- read.table('/mfs/gdouglas/projects/ocean_mags/networks/combined_tables/metaG_presence_tara_w_protist_fraction.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

genome_tax <- read.table('/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/MAG_taxa_breakdown.tsv.gz',
                         header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 2)

env_data <- read.table("/mfs/gdouglas/projects/ocean_hgt_zenodo/mapfiles/Tara_PANGEA_env_data.tsv.gz",
                       header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

env_data <- env_data[, -which(colnames(env_data) %in% c("Sample_ID_BioSamples_accession_number", "Sample_ID_ENA_sample_accession_number"))]

env_data_lower <- env_data[which(env_data$Fraction_lower_µm_used_on_board_to_prepare_samp == 0.22 & env_data$Fraction_upper_µm_used_on_board_to_prepare_samp <= 3), ]

prev_lower <- prev[rownames(env_data_lower), ]
prev_lower <- prev_lower[rowSums(prev_lower) >= 10, colSums(prev_lower) > 0]
colnames(prev_lower) <- genome_tax[colnames(prev_lower), "MAG"]

raw_lower <- list()
for (sampleid in rownames(prev_lower)) {
  sample_genomes <- colnames(prev_lower)[which(prev_lower[sampleid, ] > 0)]
  sample_genome_summary <- compute_summary_vals_simple(info[sample_genomes, ])
  raw_lower[[sampleid]] <- data.frame(sampleid = sampleid, sample_genome_summary)
}

combined_lower <- do.call(rbind, raw_lower)

env_data_upper <- env_data[which(env_data$Fraction_lower_µm_used_on_board_to_prepare_samp >= 3 & env_data$Fraction_upper_µm_used_on_board_to_prepare_samp == 20), ]

prev_upper <- prev[rownames(env_data_upper), ]
prev_upper <- prev_upper[rowSums(prev_upper) >= 10, colSums(prev_upper) > 0]
colnames(prev_upper) <- genome_tax[colnames(prev_upper), "MAG"]

raw_upper <- list()
for (sampleid in rownames(prev_upper)) {
  sample_genomes <- colnames(prev_upper)[which(prev_upper[sampleid, ] > 0)]
  sample_genome_summary <- compute_summary_vals_simple(info[sample_genomes, ])
  raw_upper[[sampleid]] <- data.frame(sampleid = sampleid, sample_genome_summary)
}

combined_upper <- do.call(rbind, raw_upper)

write.table(x = combined_lower,
            file = gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_lower_fraction_env_matched.tsv.gz"),
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

write.table(x = combined_upper,
            file = gzfile("/mfs/gdouglas/projects/ocean_hgt_zenodo/hgt_prev_analyses/Tara_upper_fraction_env_matched.tsv.gz"),
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
