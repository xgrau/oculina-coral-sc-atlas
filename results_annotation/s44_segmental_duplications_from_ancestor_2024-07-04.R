# load libs
suppressMessages(require("ape"))
suppressMessages(require("scales"))
suppressMessages(require("igraph"))
suppressMessages(require("rtracklayer"))
suppressMessages(require("UpSetR"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))

## Scope of analysis ##

prob_thr = 0.95
distance_thr = 1e5

list_tests = list(
	"allsclerplusplus" = list(sps_list = c("Ocupat","Ocuarb","Gasp","Fspp","Spis","Pocdam","Amil","Adig","Gfas"), anc_list = c("Scleractinia","Hexacorallia","Anthozoa"))
)

## Define input ##

# input
dat_fn = "../data/orthology_Metazoa_plus/orthogroup_conservation.csv"
daa_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"
ort_fn = "results_gene_family_evolution_anthozoa_plus/orthogroup_counts.posteriors.G2.csv"
phy_fn = "../data/species_tree.Anthozoa_plus.newick"
dol_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_summary.csv"

# output
out_fn = "results_segmental_duplications_plus/"
dir.create(out_fn, showWarnings = FALSE)


# read gene classification, taxonomy, etc
message("posteriors | read orthology")
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
daa = read.table(daa_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ort_names_v = dic_from_vecs(dat$orthogroup, dat$orthogroup_name)
ort_genes_v = dic_from_vecs(dat$gene, dat$orthogroup)
ort_genon_v = dic_from_vecs(dat$gene, dat$orthogroup_name)
ora_names_v = dic_from_vecs(daa$orthogroup, daa$orthogroup_name)
ora_genes_v = dic_from_vecs(daa$gene, daa$orthogroup)
ora_genon_v = dic_from_vecs(daa$gene, daa$orthogroup_name)
oga_gtv = dic_from_vecs(daa$transcript, daa$gene)
gg_v = dic_from_vecs(gsub("_","-", daa$gene), daa$gene)

# read species tree
message("posteriors | read tree")
phy = ape::read.tree(phy_fn)
sps_list = phy$tip.label
nod_list = c(phy$tip.label, phy$node.label)
root_label = phy$node.label[1]

# read posteriors matrices
message("posteriors | read posteriors processed matrices")
mat_pres = as.matrix(read.table(sprintf("%s.pres.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_pre1 = as.matrix(read.table(sprintf("%s.pres_singlecopy.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_prem = as.matrix(read.table(sprintf("%s.pres_multicopy.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_loss = as.matrix(read.table(sprintf("%s.loss.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_gain = as.matrix(read.table(sprintf("%s.gain.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_redu = as.matrix(read.table(sprintf("%s.redu.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_expa = as.matrix(read.table(sprintf("%s.expa.csv", gsub("\\.csv$","",ort_fn)), sep = "\t"))
mat_summ = read.table(sprintf("%s.summ.csv", gsub("\\.csv$","",ort_fn)), sep = "\t")

# heatmap col
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))


## Identify ancestral single-copy genes ##

for (nn in 1:length(list_tests)) {
	
	tid = names(list_tests)[nn]
	anc_list = list_tests[[nn]]$anc_list
	sps_list = list_tests[[nn]]$sps_list
	message(sprintf("posteriors | %s | %i species in %i ancestors", tid, length(sps_list), length(anc_list)))

	ogs_singlecopy_d = data.frame()
	for (ani in anc_list) {
		
		ogs_singlecopy_a = names(which(mat_pre1 [ , ani ] >= prob_thr))
		ogs_singlecopy_i = data.frame(orthogroup = ogs_singlecopy_a, ancestor = ani)
		ogs_singlecopy_d = rbind(ogs_singlecopy_d, ogs_singlecopy_i)
		
	}

	ogs_singlecopy_a = aggregate(ancestor ~ orthogroup, ogs_singlecopy_d, function(vv) paste(vv, collapse = ","))


	## Identify extant multi-copy genes (candidate SDs) ##

	gr_list = vector(mode = "list", length = length(sps_list))
	names(gr_list) = sps_list
	for (spi in sps_list) {
		
		# load gene coordinates, add orthogroups
		gti = rtracklayer::import(sprintf("../data/reference/%s_long.annot.gtf", spi))
		ggi = gti [ gti$type == "transcript" ]
		ggi$orthogroup = ora_genes_v [ ggi$gene_id ]
		ggi$orthogroup = factor(ggi$orthogroup , levels = unique(ora_genes_v))
		gr_list[[spi]] = ggi
		
		count_ogs = table(ggi$orthogroup)
		ogs_singlecopy_a[,spi] = count_ogs [ ogs_singlecopy_a$orthogroup ]

	}

	ogs_singlecopy_a_f = ogs_singlecopy_a [ apply(ogs_singlecopy_a[,sps_list], 1, function(vv) any(vv > 1)), ]
	ogs_singlecopy_a_f_v = dic_from_vecs(ogs_singlecopy_a_f$orthogroup, ogs_singlecopy_a_f$ancestor)


	## Identify segmental duplications ##

	segments_d = data.frame()
	for (spi in sps_list) {
		
		ogs_expanded_i = ogs_singlecopy_a_f$orthogroup [ ogs_singlecopy_a_f[,spi] > 1 ]
		ggi = gr_list[[spi]]
		message(sprintf("posteriors | %s | n=%i ancestral single-copy families expanded in %s", tid, length(ogs_expanded_i), spi))
		
		segments_i = data.frame()
		for (ogi in ogs_expanded_i) {
			
			ggi_i = ggi [ ggi$orthogroup == ogi ]
			
			# pairwise distances between all paralogs
			pairwise_i = sapply(1:length(ggi_i), function(vv) { GenomicRanges::distance(ggi_i[vv], ggi_i, ignore.strand = TRUE) })
			diag(pairwise_i) = NA
			
			# which paralogs are close enough? (below threshold)
			pairwise_i_m = apply(pairwise_i, 1, function(vv) vv < distance_thr )
			pairwise_i_m = pairwise_i_m * 1
			pairwise_i_m [ is.na(pairwise_i_m) ] = 0
			rownames(pairwise_i_m) = colnames(pairwise_i_m) = ggi_i$gene_id
			pairwise_i_g = igraph::graph_from_adjacency_matrix(pairwise_i_m)
			pairwise_i_g_components = igraph::components(pairwise_i_g)
			pairwise_i_d = data.frame(
				duplicate_set = pairwise_i_g_components$membership,
				gene = names(pairwise_i_g_components$membership),
				chromosome = as.character(seqnames(ggi_i)),
				start = start(ggi_i),
				end = end(ggi_i),
				orthogroup = ogi,
				species = spi,
				gene_name = ort_genon_v [ names(pairwise_i_g_components$membership) ]
			)
			
			rownames(pairwise_i_d) = NULL
			pairwise_i_d = pairwise_i_d [ pairwise_i_d$duplicate_set %in% which(table(pairwise_i_d$duplicate_set) > 1), ]
			pairwise_i_d$duplicate_set = sprintf("%s:%s:set%03d", ogi,spi, as.numeric(as.factor(pairwise_i_d$duplicate_set)))
			# pairwise_i_num_segmental_dups = length(unique(pairwise_i_d$duplicate_set))
			
			if (nrow(pairwise_i_d)>0) {
				segments_i = rbind(segments_i, pairwise_i_d)
			}
			
		}
		
		segments_d = rbind(segments_d, segments_i)
		message(sprintf("posteriors | %s | n=%i ancestral single-copy families expanded in %s | n=%i families duplicated in %i segments", tid, length(ogs_expanded_i), spi, length(unique(segments_i$orthogroup)), length(unique(segments_i$duplicate_set))))
		
	}


	# summarise overlaps
	pdf(sprintf("%s/overlap.%s.venn.OcuSpiAmi.pdf", out_fn, tid), width = 5, height = 5)
	vvi = venn.three(
		unique(segments_d[segments_d$species=="Ocupat","orthogroup"]), 
		unique(segments_d[segments_d$species=="Spis","orthogroup"]), 
		unique(segments_d[segments_d$species=="Amil","orthogroup"]), 
		catname1 = "Ocupat",
		catname2 = "Spis",
		catname3 = "Amil",
		main = "")
	dev.off()

	# summarise overlaps
	pdf(sprintf("%s/overlap.%s.venn.OcuOcu.pdf", out_fn, tid), width = 5, height = 5)
	vvo = venn.two(
		unique(segments_d[segments_d$species=="Ocupat","orthogroup"]), 
		unique(segments_d[segments_d$species=="Ocuarb","orthogroup"]), 
		catname1 = "Ocupat",
		catname2 = "Ocuarb",
		main = "")
	dev.off()

	# save list
	segments_d$species = factor(segments_d$species, levels = sps_list)
	segments_d = segments_d [ order(segments_d$orthogroup,segments_d$species), ]
	segments_d$single_copy_in = ogs_singlecopy_a_f_v [ segments_d$orthogroup ]
	write.table(segments_d, file = sprintf("%s/summary_segmental_duplications.%s.csv", out_fn, tid), sep = "\t", quote = FALSE, row.names = FALSE)

	# # reload
	# segments_d = read.table(sprintf("%s/summary_segmental_duplications.%s.csv", out_fn, tid), sep = "\t", header = TRUE)
	# segments_d$species = factor(segments_d$species, levels = list_tests[[tid]]$sps_list)
	
	# aggregate at the duplicate set level
	segments_d_a = segments_d [ , c("duplicate_set","species","orthogroup","single_copy_in") ]
	segments_d_a = unique(segments_d_a)

	# aggregate at the orthogroup level
	segments_d_aa = segments_d [ , c("species","orthogroup","single_copy_in") ]
	segments_d_aa = unique(segments_d_aa)

	# some denominators for counts
	num_extant_expanded = apply(ogs_singlecopy_a_f[,sps_list], 2, function(vv) sum(vv > 1))
	num_extant_families = table(unique(daa [ c("orthogroup","species") ])$species)
	num_extant_genes = table(daa$species)

	# summary counts
	pdf(sprintf("%s/summary_segmental_duplications.%s.per_sps.pdf", out_fn, tid), width = 16, height = 8)
	layout(matrix(1:9, nrow = 3, byrow = TRUE))

	hh = table(segments_d$species)
	b = barplot(
		hh,
		col = "lightblue3", las = 2,
		main = "Num genes in segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, hh, pos = 1, cex = 0.8, col = "darkblue")	
		
	hh = table(segments_d_a$species)
	b = barplot(
		hh,
		col = "lightblue3", las = 2,
		main = "Num segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, hh, pos = 1, cex = 0.8, col = "darkblue")		
		
	hh = table(segments_d_aa$species)
	b = barplot(
		hh,
		col = "lightblue3", las = 2,
		main = "Num expanded single-copy ancestral families in segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, hh, pos = 1, cex = 0.8, col = "darkblue")	

	hh = num_extant_expanded 
	b = barplot(
		num_extant_expanded,
		col = "lightblue3", las = 2,
		main = "Num expanded single-copy ancestral families with any kind of duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, hh, pos = 1, cex = 0.8, col = "darkblue")

	hh = table(segments_d_aa$species) / num_extant_expanded
	b = barplot(
		hh,
		col = "lightblue3", las = 2, ylim = c(0,1),
		main = "Fraction expanded single-copy ancestral families in segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, sprintf("%.3f", hh), pos = 1, cex = 0.8, col = "darkblue")

	hh = table(segments_d_aa$species) / num_extant_families [ sps_list ]
	b = barplot(
		hh,
		col = "lightblue3", las = 2, 
		main = "Fraction of all extant families in segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, sprintf("%.3f", hh), pos = 1, cex = 0.8, col = "darkblue")

	hh = table(segments_d$species) / num_extant_genes [ sps_list ]
	b = barplot(
		hh,
		col = "lightblue3", las = 2, 
		main = "Fraction of all extant genes in segmental duplications", cex.main = 0.7, font.main = 1)
	text(b, hh, sprintf("%.3f", hh), pos = 1, cex = 0.8, col = "darkblue")

	# close
	dev.off()


	## Set overlap analysis ##
	
	# # setup upset dataframe, for all species in test
	# xtt = table(segments_d$orthogroup, segments_d$species)
	# xtm = as.data.frame(table_to_matrix(xtt))
	# xtl = apply(xtm, 2, function(v) names(which(v>0)))
	# xtd = UpSetR::fromList(xtl)
	# # plot
	# uu1 = UpSetR::upset(
	# 	xtd,
	# 	nsets = ncol(xtd),
	# 	keep.order = TRUE,
	# 	sets = colnames(xtd),
	# 	cutoff = 0,
	# 	order.by = "freq",
	# 	mainbar.y.label = "shared orthologs with segmental duplications",
	# 	sets.x.label = "num orthologs with segmental duplications",
	# 	matrix.color = "black",
	# 	main.bar.color = "lightblue3", 
	# 	sets.bar.color = "lightblue4"
	# )

	# # setup upset dataframe, for a subset of species
	# sps_focus = c("Ocupat","Spis","Amil","Nvec","Xesp")
	# sps_focus_f = sps_focus [ sps_focus %in% sps_list ]
	# xtl_f = xtl [ sps_focus_f ]
	# xtd = UpSetR::fromList(xtl_f)
	# # plot
	# uu2 = UpSetR::upset(
	# 	xtd,
	# 	nsets = ncol(xtd),
	# 	keep.order = TRUE,
	# 	sets = colnames(xtd),
	# 	cutoff = 0,
	# 	order.by = "freq",
	# 	mainbar.y.label = "shared families with segmental duplications",
	# 	sets.x.label = "num families with segmental duplications",
	# 	matrix.color = "black",
	# 	main.bar.color = "lightblue3", 
	# 	sets.bar.color = "lightblue4"
	# )
	
	# pdf(sprintf("%s/summary_segmental_duplications.%s.overlaps.pdf", out_fn, tid), width = 10, height = 5)
	# print(uu1)
	# print(uu2)
	# dev.off()


	## Functional enrichments ##

	dir.create(sprintf("%s/tests_%s_nov24", out_fn, tid), showWarnings = FALSE)


	for (spi in sps_list [ sps_list %in% c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp") ]) {

		# functional enrichments
		go_annot = gsa_topgo_load_emapper(emapper_fn = sprintf("../data/reference/%s_ensembl.GO.csv", spi), index_col_GOs = 2)
		pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi))
		pf_annot_fullarch = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi), do_unique = FALSE)
		names(go_annot) = oga_gtv [ names(go_annot) ]
		names(pf_annot) = oga_gtv [ names(pf_annot) ]
		names(pf_annot_fullarch) = oga_gtv [ names(pf_annot_fullarch) ]

		# load KO info in species of interest
		keg_lll = gsa_enrichment_load_pfam_list(sprintf("../results_scatlas/results_pathways/kegg_gene_to_map_%s.csv", spi), architecture_sep = ",")
		names(keg_lll) = gg_v [ names(keg_lll) ]


		# functional enrichments, first group
		fgi_genes = segments_d [ segments_d$species == spi, "gene" ]
		fgi_genes = fgi_genes [ !is.na(fgi_genes) ]
		bgi_genes = unique(daa [ daa$species == spi, "gene" ])
		
		message(sprintf("functional enrichments | %s | %s n=%i genes", tid, spi, length(fgi_genes)))
		if (length(fgi_genes)>0) {
			
			# KEGG enrichment
			kegg_enrichment = gsa_enrichment_hypergeometric_v2(
				annotation = keg_lll,
				genes_fg = fgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/all.%s.kegg.csv", out_fn, tid, spi)
			)

			# test pfam
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/all.%s.pfam.csv", out_fn, tid, spi),
			)
			
			list_arch = pf_annot_fullarch [ fgi_genes [ fgi_genes %in% names(pf_annot_fullarch) ] ]
            list_arch_v = sapply(list_arch, function(v) paste(v, collapse = "/"))
            list_arch_v = sort(list_arch_v)
            write.table(list_arch_v, sprintf("%s/tests_%s_nov24/all.%s.pfam_architectures.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)
			
			# test GOs
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/all.%s.topgo.csv", out_fn, tid, spi),
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)

			# integrate GOs and pfam in a single table
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("%s/tests_%s_nov24/all.%s.summary.csv", out_fn, tid, spi),
				extra_annot_dict = ort_genon_v
			) 
			
		}
		
		# enrichment in the intersect (using this sps as ref)
		fgi_genes = segments_d [ segments_d$species == spi & segments_d$orthogroup %in% vvi$intersect_123, "gene" ]
		fgi_genes = fgi_genes [ !is.na(fgi_genes) ]
		bgi_genes = unique(daa [ daa$species == spi, "gene" ])

		message(sprintf("functional enrichments | %s | %s n=%i genes in intersect", tid, spi, length(fgi_genes)))
		if (length(fgi_genes)>0) {
			
			# KEGG enrichment
			kegg_enrichment = gsa_enrichment_hypergeometric_v2(
				annotation = keg_lll,
				genes_fg = fgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/intersect.%s.kegg.csv", out_fn, tid, spi)
			)

			# test pfam
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/intersect.%s.pfam.csv", out_fn, tid, spi)
			)
			
			list_arch = pf_annot_fullarch [ fgi_genes [ fgi_genes %in% names(pf_annot_fullarch) ] ]
            list_arch_v = sapply(list_arch, function(v) paste(v, collapse = "/"))
            list_arch_v = sort(list_arch_v)
            write.table(list_arch_v, sprintf("%s/tests_%s_nov24/intersect.%s.pfam_architectures.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

			# test GOs
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/intersect.%s.topgo.csv", out_fn, tid, spi),
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)

			# integrate GOs and pfam in a single table
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("%s/tests_%s_nov24/intersect.%s.summary.csv", out_fn, tid, spi),
				extra_annot_dict = ort_genon_v
			) 
			
		}
		
		
		if (spi %in% c("Ocupat","Ocuarb")) { 
			
			# enrichment in the oculina intersect (using this sps as ref)
			fgi_genes = segments_d [ segments_d$species == spi & segments_d$orthogroup %in% vvo$list_intersect, "gene" ]
			fgi_genes = fgi_genes [ !is.na(fgi_genes) ]
			bgi_genes = unique(daa [ daa$species == spi, "gene" ])

			message(sprintf("functional enrichments | %s | %s n=%i genes in oculina intersect", tid, spi, length(fgi_genes)))
			if (length(fgi_genes)>0) {
				
				# KEGG enrichment
				kegg_enrichment = gsa_enrichment_hypergeometric_v2(
					annotation = keg_lll,
					genes_fg = fgi_genes,
					out_fn = sprintf("%s/tests_%s_nov24/intersectocuocu.%s.kegg.csv", out_fn, tid, spi)
				)

				# test pfam
				pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
					annotation = pf_annot,
					genes_fg = fgi_genes,
					genes_bg = bgi_genes,
					out_fn = sprintf("%s/tests_%s_nov24/intersectocuocu.%s.pfam.csv", out_fn, tid, spi)
				)
				
				list_arch = pf_annot_fullarch [ fgi_genes [ fgi_genes %in% names(pf_annot_fullarch) ] ]
				list_arch_v = sapply(list_arch, function(v) paste(v, collapse = "/"))
				list_arch_v = sort(list_arch_v)
				write.table(list_arch_v, sprintf("%s/tests_%s_nov24/intersectocuocu.%s.pfam_architectures.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

				# test GOs
				go_enrichment_l = gsa_topgo_enrichment_v2(
					annotation = go_annot,
					genes_fg = fgi_genes,
					genes_bg = bgi_genes,
					out_fn = sprintf("%s/tests_%s_nov24/intersectocuocu.%s.topgo.csv", out_fn, tid, spi),
					ontologyset = c("BP","MF","CC"),
					tg_test = "fisher",
					tg_algorithm = "elim",
					top_markers = 30,
					nodesize = 10
				)

				# integrate GOs and pfam in a single table
				gsa_merge_tables_v2(
					main_enrichment = go_enrichment_l[[1]],
					main_mapping = go_enrichment_l[[2]],
					pfam_enrichment = pf_enrichment_l[[1]],
					pfam_mapping = pf_enrichment_l[[2]],
					output_fn = sprintf("%s/tests_%s_nov24/intersectocuocu.%s.summary.csv", out_fn, tid, spi),
					extra_annot_dict = ort_genon_v
				) 
				
			}
		}
		
		# enrichment in the non-intersect (using this sps as ref)
		fgi_genes = segments_d [ segments_d$species == spi & !segments_d$orthogroup %in% vvi$intersect_123, "gene" ]
		fgi_genes = fgi_genes [ !is.na(fgi_genes) ]
		bgi_genes = unique(daa [ daa$species == spi, "gene" ])

		message(sprintf("functional enrichments | %s | %s n=%i genes in non-intersect", tid, spi, length(fgi_genes)))
		if (length(fgi_genes)>0) {
			
			# KEGG enrichment
			kegg_enrichment = gsa_enrichment_hypergeometric_v2(
				annotation = keg_lll,
				genes_fg = fgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/nonintersect.%s.kegg.csv", out_fn, tid, spi)
			)

			# test pfam
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/nonintersect.%s.pfam.csv", out_fn, tid, spi)
			)
			
			list_arch = pf_annot_fullarch [ fgi_genes [ fgi_genes %in% names(pf_annot_fullarch) ] ]
            list_arch_v = sapply(list_arch, function(v) paste(v, collapse = "/"))
            list_arch_v = sort(list_arch_v)
            write.table(list_arch_v, sprintf("%s/tests_%s_nov24/nonintersect.%s.pfam_architectures.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = TRUE, col.names = FALSE)

			# test GOs
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = fgi_genes,
				genes_bg = bgi_genes,
				out_fn = sprintf("%s/tests_%s_nov24/nonintersect.%s.topgo.csv", out_fn, tid, spi),
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)

			# integrate GOs and pfam in a single table
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("%s/tests_%s_nov24/nonintersect.%s.summary.csv", out_fn, tid, spi),
				extra_annot_dict = ort_genon_v
			) 
			
			
		}
		
	}


}
