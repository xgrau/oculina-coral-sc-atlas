# libraries
library("rtracklayer")
library("igraph")
source("../scripts/helper.R")
graphics.off()


# list of species
sps_list = c("Ocupat","Amil","Fspp","Gasp","Spis")

# reference species (Tadh, it's got the best assembly)
sps_ref = "Ocupat"

# files
ogs_fn = "../data/orthology_Anthozoa/dir_step4/orthologous_pairs.txt"
out_fn = "results_wga/"

# load orthology
ogs = read.table(ogs_fn, col.names = c("gene1", "gene2"))

ogs_g   = igraph::graph_from_data_frame(ogs, directed = FALSE)
ogs_g_c = igraph::components(ogs_g)
ogs_g_d = data.frame(orthogroup = ogs_g_c$membership, transcript_id = names(ogs_g_c$membership))

# load gffs
list_gff = list()
list_fas = list()
for (spi in sps_list) {
	
	message(sprintf("pairs | load %s", spi))
	gfi = rtracklayer::readGFF(sprintf("../data/reference/%s_long.annot.gtf", spi))
	list_gff[[spi]] = gfi [ gfi$type == "transcript", ]
	list_fas[[spi]] = ape::read.FASTA(sprintf("../data/reference/%s_gDNA.fasta", spi))
	
}

# loop and crosstabulate
for (spi in sps_ref) {

	list_hms = list()
	list_xts = list()
	for (spj in sps_list [ sps_list != spi ]) {

		message(sprintf("pairs | find homologous chromosomes | %s-%s", spi, spj))

		# retrieve pairs
		ci = list_gff[[spi]][,c("seqid","transcript_id")]
		cj = list_gff[[spj]][,c("seqid","transcript_id")]
		
		# for reference species, keep large scaffolds only
		ci_large = rev(names(which(sort(table(ci$seqid)) >= 50)))
		ci_large = stringr::str_sort(ci_large, numeric = TRUE)
		
		# ignore chromosomes with <X genes
		ci = ci [ ci$seqid %in% names(which(table(ci$seqid) >= 100)) , ]
		cj = cj [ cj$seqid %in% names(which(table(cj$seqid) >= 100)) , ]
		
		# add fine-grained orthology
		ci = merge(ci, ogs_g_d, by = "transcript_id", all.x = FALSE, all.y = FALSE)
		cj = merge(cj, ogs_g_d, by = "transcript_id", all.x = FALSE, all.y = FALSE)
		
		# merge sps i and j
		cm = merge(ci, cj, by = "orthogroup", suffixes = c("_spi","_spj"), all.x = FALSE, all.y = FALSE)
		
		# crosstabulate
		xt = table(cm$seqid_spi, cm$seqid_spj)
		xt_f = t(t(xt) / colSums(xt))
		xt_f = xt_f [ ci_large, ]
		xt_f [ is.na(xt_f) ] = 0
		xt_f = xt_f [ , apply(xt_f, 2, function(c) max(c) > 0.05) ]
		
		# save crosstabulation matrix 
		list_xts[[spj]] = xt_f
		
		ordered_chr_spj = select_top_markers(t(xt_f), n_top_markers = 1e6, n_markers_rollmean = 1)
		
		# heatmap
		list_hms[[spj]] = plot_complex_heatmap(
			t(xt_f)[ ordered_chr_spj, ],
			name = "fraction",
			title_row = sprintf("%s, n=%i", spj, nrow(t(xt_f))),
			title_col = sprintf("%s, n=%i", spi, ncol(t(xt_f))),
			fontsize = 5,
			color_mat = c("gray95","#d6e72e","#6fb600","#003f4d"),
			cell_border = gpar(col = "white", lwd = 1, lty = 1),
			heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
			color_min = 0.05,
			color_max = 1,
			use_raster = FALSE,
			cluster_row = FALSE,
			cluster_col = FALSE)
			
		# list of homologs for each chrom in spi
		xt_d = as.data.frame(xt_f)
		colnames(xt_d) = c("spi","spj","fraction_orthologs_j_in_i")
		xt_c = as.data.frame(xt)
		colnames(xt_c) = c("spi","spj","num_orthologs_j_in_i")
		xt_d$merge_id = paste(xt_d$spi, xt_d$spj)
		xt_c$merge_id = paste(xt_c$spi, xt_c$spj)
		xt_m = merge(xt_d, xt_c[,c("merge_id","num_orthologs_j_in_i")], by = "merge_id", all.x = TRUE, all.y = FALSE)
		xt_m$genes_in_i = table(cm$seqid_spi) [ as.character(xt_m$spi) ]
		xt_m$genes_in_j = table(cm$seqid_spj) [ as.character(xt_m$spj) ]
		xt_m = xt_m [ xt_m$fraction_orthologs_j_in_i > 0.05, ]
		xt_m = xt_m [ order(xt_m$spi, -xt_m$num_orthologs_j_in_i), ]
		xt_m$merge_id = NULL
		
		# rename chromosomes to add species name
		xt_m$spi = paste(spi, xt_m$spi, sep = ".")
		xt_m$spj = paste(spj, xt_m$spj, sep = ".")
		
		# save
		xt_m = xt_m [ order(xt_m$spi), ]
		write.table(xt_m, sprintf("%s/chrom_homology.%s.to_%s.csv",out_fn, spi, spj), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
		
	}
	
	message(sprintf("pairs | find homologous chromosomes | plot"))
	pdf(sprintf("%s/chrom_homology.%s.pdf",out_fn, spi), width = 4, height = 120)
	print(list_hms[[1]] %v% list_hms[[2]] %v% list_hms[[3]] %v% list_hms[[4]])
	dev.off()

}
