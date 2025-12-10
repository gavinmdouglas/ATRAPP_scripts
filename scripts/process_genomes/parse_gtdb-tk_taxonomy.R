rm(list = ls(all.names = TRUE))

bac120_taxonomy_file <- '/mfs/nicot/mybook2/MAGs_ATRAPP/dRep_out_MAGs_last_98/gtdbtk.bac120.summary.tsv'
ar53_taxonomy_file <- '/mfs/nicot/mybook2/MAGs_ATRAPP/dRep_out_MAGs_last_98/gtdbtk.ar53.summary.tsv'

bac120_taxonomy <- read.table(bac120_taxonomy_file, stringsAsFactors = FALSE, sep = '\t', comment.char = '', quote = '', header=TRUE)
ar53_taxonomy <- read.table(ar53_taxonomy_file, stringsAsFactors = FALSE, sep = '\t', comment.char = '', quote = '', header=TRUE)

gtdb_classifications <- rbind(ar53_taxonomy, bac120_taxonomy)

taxonomy_tab <- data.frame(matrix(NA, nrow = nrow(gtdb_classifications), ncol=8))
colnames(taxonomy_tab) <- c('MAG', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

taxonomy_tab$MAG <- sub('\\.fa\\.cleaned$', '', gtdb_classifications$user_genome)
taxonomy_tab$MAG <- sub('\\.cleaned$', '', taxonomy_tab$MAG)
taxonomy_tab$MAG <- sub('\\.fa\\.cl$', '', taxonomy_tab$MAG)
taxonomy_tab$MAG <- sub('\\.cl$', '', taxonomy_tab$MAG)

taxonomy_tab$Domain <- sub(';p__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Phylum <- sub(';c__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Class <- sub(';o__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Order <- sub(';f__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Family <- sub(';g__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Genus <- sub(';s__.*$', '', gtdb_classifications$classification)
taxonomy_tab$Species <- gtdb_classifications$classification

write.table(x = taxonomy_tab,
            file = gzfile('/mfs/gdouglas/projects/ATRAPP/genome_process/ATRAPP_MAG_GTDBtk_taxonomy.tsv.gz'),
            sep = '\t', quote=FALSE, row.names = FALSE, col.names = TRUE)
