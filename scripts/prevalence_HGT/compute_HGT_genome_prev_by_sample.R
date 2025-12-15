rm(list = ls(all.names = TRUE))

# Create tables breaking down proportion of genomes per sample with at least 1 HGT event
# (or one genome with a gene found across multiple genera, etc.)

prev <- read.table('/mfs/gdouglas/projects/ATRAPP/cooccur/df_cv.tsv.gz',
                   header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1,
                   check.names = FALSE)

rownames(prev) <- gsub(".fa.cleaned", "", rownames(prev))
rownames(prev) <- gsub(".fa.cl", "", rownames(prev))
rownames(prev) <- gsub(".cl", "", rownames(prev))

# Convert to prevalence (0's and 1's)
prev[prev > 0] <- 1

genome_tax <- read.table('/mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz',
                         header=TRUE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

# Make sure that all prev rownames match expected format.
if (length(setdiff(rownames(prev), rownames(genome_tax))) > 0) {
  stop('Some genomes in prev not found in genome_tax!')
}

HGT_info <- read.table("/mfs/gdouglas/projects/ATRAPP/HGT_prevalence/genome_HGT_info.tsv.gz",
                       header=TRUE, sep="\t", row.names=1)

raw <- list()

label_map <- c(
  "HGT_num" = "HGT_prop",
  "Cross_Genus" = "Cross_Genus_prop",
  "Cross_Family" = "Cross_Family_prop",
  "Cross_Order" = "Cross_Order_prop",
  "Cross_Class" = "Cross_Class_prop",
  "Cross_Phylum" = "Cross_Phylum_prop",
  "Cross_Domain" = "Cross_Domain_prop"
)

for (samp in colnames(prev)) {
  samp_genomes <- rownames(prev)[prev[, samp] > 0]
  num_genomes <- length(samp_genomes)

  # Get subset of HGT_info
  samp_genomes_intersect <- intersect(samp_genomes, rownames(HGT_info))

  if (length(samp_genomes_intersect) > 0) {
    samp_HGT_info <- HGT_info[samp_genomes_intersect, , drop = FALSE]
  } else {
    samp_HGT_info <- HGT_info[0, ]
  }

  samp_HGT_info_prop <- colSums(samp_HGT_info > 0) / num_genomes

  raw[[samp]] <- data.frame(sample=samp,
                            num_genomes=num_genomes)

  for (coln in names(samp_HGT_info_prop)) {
    raw[[samp]][, label_map[[coln]]] <- samp_HGT_info_prop[[coln]]
  }

}

combined <- do.call(rbind, raw)

write.table(x = combined,
            file = gzfile("/mfs/gdouglas/projects/ATRAPP/HGT_prevalence/mgs_sample_hgt_prev.tsv.gz"),
            sep = "\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
