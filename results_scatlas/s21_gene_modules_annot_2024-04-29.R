# libraries
suppressMessages(library("scales"))
suppressMessages(source("../scripts/Gene_module_functions.R"))
suppressMessages(source("../scripts/Cross_species_functions.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))

# data index
sps_list = c("Ocupat","Ocuarb","Spin","Amil","Nvec","Xesp")

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)
oga_v = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)

# load gene names
gna_fn = "../data/orthology_Metazoa_plus/orthogroup_conservation.csv"
gna = read.table(gna_fn, sep = "\t", header = TRUE)
ogm_gv = dic_from_vecs(gna$gene, gna$gene_name)
ogm_gvo = dic_from_vecs(gna$gene, gna$orthogroup_name)
# dictionaries transcript to gene
gna_gtv = dic_from_vecs(names = gna$transcript, terms = gna$gene)
gna_tgv = dic_from_vecs(names = gna$gene, terms = gna$transcript)


# dataframes to carry
mc_wgcna_gmods_annot_t = data.frame()
module_annot_t = data.frame()
for (spi in sps_list) {

	# set working species
	if (spi == "Spio" | spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}
	
	# output
	out_fn = sprintf("results_metacell_%s_filt/", spi)
	dir.create(out_fn, showWarnings = FALSE)

	## Load annotations ##

	# load gene annotations for this species
	gna_i = gna [ gna$species == spi_w, ]
	gna_p = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi_w), sep = "\t", col.names = c("transcript","blast","pfam"))
	gna_p$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi_w), vector_to_fix = gna_p$transcript)
	gna_i = merge(gna_i, gna_p[,c("gene","pfam")], by.x = "gene", by.y = "gene", all.x = TRUE, all.y = TRUE)
	gene_annot = data.frame(gene = gna_i$gene, name = gna_i$orthogroup_name, pfam = gna_i$pfam)
	gene_annot = gene_annot [ !duplicated(gene_annot$gene), ]
	rownames(gene_annot) = gene_annot$gene
	
	# load TF annotations
	gene_annot_tfs = ogm [ grepl(":", ogm$orthogroup_name) & ogm$species == spi_w, ]
	gene_annot_tfs_v = dic_from_vecs(gene_annot_tfs$gene, gene_annot_tfs$orthogroup_name)
	gene_annot$name [ rownames(gene_annot) %in% gene_annot_tfs$gene ] = paste(gene_annot_tfs_v [ rownames(gene_annot) [ rownames(gene_annot) %in% gene_annot_tfs$gene ] ] , gene_annot$name [ rownames(gene_annot) %in% gene_annot_tfs$gene ], sep = ":")
	gene_annot$is_tf = gene_annot$gene %in% gene_annot_tfs$gene
	
	# seurat gene dict
	gg_v = dic_from_vecs(gsub("_","-", gna_p$gene), gna_p$gene)

	message(sprintf("wgcna | %s | load footprints...", spi))
	mc_fp = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", out_fn, spi))
	mc_uf = readRDS(sprintf("%s/dat.%s.expression.mcs_umifrac.rds", out_fn, spi))
	
	message(sprintf("wgcna | %s | load annotations...", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	ctm = read.table(sprintf("%s/annot.%s.mcs.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
	ctt_bct_v = dic_from_vecs(ctt$cell_type, ctt$broad_cell_type)

	# reload
	message(sprintf("wgcna | %s | reload WGCNA...", spi))
	# mc_wgcna = readRDS(sprintf("%s/gmod.%s.wgcna.wgcna.rds", out_fn, spi))
	mc_wgcna_gmods = readRDS(sprintf("%s/gmod.%s.wgcna.memberships.rds", out_fn, spi))
	mc_wgcna_gmods_annot = read.table(sprintf("%s/gmod.%s.wgcna.memberships.csv", out_fn, spi), sep = "\t", header = TRUE)
	mc_wgcna_me = readRDS(sprintf("%s/gmod.%s.wgcna.ME.rds", out_fn, spi))
	mc_wgcna_me = mc_wgcna_me [ ctm$metacell, ]
	
	# plot mc fp
	gene_annot_gm = gene_annot
	gene_annot_gm$gene = NULL
	rownames(gene_annot_gm) = gsub("_","-", rownames(gene_annot_gm))
	gmod_plotMCheatmap_annotate_modules(
		expr_matrix = mc_fp,
		gmods = mc_wgcna_gmods,
		me = mc_wgcna_me,
		expr_matrix_colors = ctm_mcs_col_v,
		ex_output_file = "/dev/null",
		an_output_file = "/dev/null",
		me_output_file = sprintf("%s/gmod.%s.gmod_eigengenes.pdf", out_fn, spi),
		ex_width = 200, ex_height = 100,
		me_width = (10 + nrow(mc_wgcna_me) / 20), me_height = (10 + ncol(mc_wgcna_me) / 20),
		resolution_rate = 1, 
		do_expression = FALSE,
		eigen_min = 0,
		cor_cutoff_max = NULL,
		annotation = gene_annot_gm,
		highlight_genes = rownames(gene_annot[gene_annot$is_tf,]), 
		heatmap_colors = c("gray95","orange","orangered2","#520c52"),
		heatmap_colors_cor = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		highlight_genes_annot = gene_annot_gm[gene_annot_gm$is_tf,"name"]
	)
	
	
	# unlist_gmods = as.character(unlist(mc_wgcna_gmods))
	# mc_wgcna_me_signed = WGCNA::signedKME(t(mc_fp[unlist_gmods, ]), datME = mc_wgcna_me)

	
	# annotate modules to cell types
	colnames(mc_wgcna_me) = gsub("^ME","",colnames(mc_wgcna_me))
	module_annot = data.frame(row.names = colnames(mc_wgcna_me), module = colnames(mc_wgcna_me))
	pdf(sprintf("%s/gmod.%s.wgcna.gmod_annotation.activity_per_mc.pdf", out_fn, spi), width = 8, height = 5)
	layout(matrix(1:2, ncol = 2, byrow = TRUE))
	for (moi in colnames(mc_wgcna_me)) {

		# sort metacells by eigenvalue
		mc_eigen_d = data.frame(metacell = rownames(mc_wgcna_me), eigen = mc_wgcna_me[,moi])
		mc_eigen_d = merge(mc_eigen_d, ctm[,c("metacell","metacell_annotation")], by.x = "metacell", by.y = "metacell", all.x = TRUE, all.y = FALSE)
		mc_eigen_d = mc_eigen_d [ order(mc_eigen_d$eigen, decreasing = TRUE), ]
		
		# mc_eigen_d_f = head(mc_eigen_d,100)
		# mc_eigen_d_f_a = aggregate(eigen ~ metacell_annotation, mc_eigen_d_f, sum)
		# mc_eigen_d_f_a = mc_eigen_d_f_a [ order(mc_eigen_d_f_a[,2], decreasing = TRUE), ]
		# mc_eigen_d_f_a = mc_eigen_d_f_a [ mc_eigen_d_f_a[,2]>0, ]
		# moi_tag_full = paste(sprintf("%s(%.2f)",mc_eigen_d_a[,1],mc_eigen_d_a[,2]), collapse = ";")
		# moi_tag_top  = head(mc_eigen_d_a[,1], 1)
		
		# scale eigenvalue
		mc_eigen_d$eigen_scaled = as.numeric(scale(mc_eigen_d$eigen))
		# select top mcs
		if (any(mc_eigen_d$eigen_scaled >= qnorm(0.99))) {
			qqp = 0.99
		} else if (any(mc_eigen_d$eigen_scaled >= qnorm(0.95))) {
			qqp = 0.95
		} else {
			qqp = 0.75
		}
		index_thr = max( max(which(mc_eigen_d$eigen_scaled >= qnorm(qqp))) )
		
		# find tag based on top metacell
		mc_eigen_d_f = mc_eigen_d[ 1:index_thr, ]
		mc_eigen_d_a = aggregate(eigen ~ metacell_annotation, data = mc_eigen_d_f, sum)
		mc_eigen_d_a = mc_eigen_d_a [ order(mc_eigen_d_a[,2], decreasing = TRUE), ]
		mc_eigen_d_a = mc_eigen_d_a [ mc_eigen_d_a[,2] / mc_eigen_d_a[1,2] > 0.1, ]
		moi_tag_full = paste(sprintf("%s(%.2f)",mc_eigen_d_a[,1],mc_eigen_d_a[,2]), collapse = ";")
		moi_tag_top  = head(mc_eigen_d_a[,1], 1)
		
		plot(
			head(mc_eigen_d$eigen, 30),
			col = ctm_mcs_col_v[head(mc_eigen_d$metacell, 30)],
			las = 1, ylab = "eigenvalue", pch = 19,
			ylim = c(0,1))
		title(main = sprintf("top metacells at p=%.2f in %s module\n%s", qqp, moi, moi_tag_full), cex.main = 0.7, font.main = 1, adj = 0)
		abline(v = index_thr+0.5, col = "darkred", lty = 2)
		text(1:30, head(mc_eigen_d$eigen, 30)+0.01, head(mc_eigen_d$metacell_annotation, 30), srt = 90, cex = 0.5, adj = 0, col = "thistle4")
		plot(
			head(mc_eigen_d$eigen_scaled, 30),
			col = ctm_mcs_col_v[head(mc_eigen_d$metacell,30)],
			las = 1, ylab = "z-score", pch = 19,
			cex.main = 0.7, font.main = 1
		)
		abline(h = qnorm(c(0.95,0.99)), lty = 2, col = "darkgray")
		abline(v = index_thr+0.5, col = "darkred", lty = 2)
		text(1:30, head(mc_eigen_d$eigen_scaled, 30), head(mc_eigen_d$metacell_annotation, 30), srt = 90, cex = 0.5, adj = 0, col = "thistle4")
		
		
		# annotate module expression pattern
		module_annot[moi,"annotation_top"] = moi_tag_top
		module_annot[moi,"annotation_string"] = moi_tag_full
		module_annot[moi,"metacell_top"] = mc_eigen_d$metacell[1]
		module_annot[moi,"num_active_mcs_total"] = as.numeric(index_thr)
		module_annot[moi,"num_active_mcs_in_top"] = sum(mc_eigen_d_f$metacell_annotation == moi_tag_top)
		module_annot[moi,"frac_active_mcs_in_top"] = sprintf("%.3f",sum(mc_eigen_d_f$metacell_annotation == moi_tag_top) / sum(ctm$metacell_annotation == moi_tag_top))
		module_annot[moi,"is_transversal"] = as.numeric(length(unique(ctt_bct_v [ mc_eigen_d_a$metacell_annotation ])) > 1)
		# is general?		
		median_eigen_all = median(mc_eigen_d$eigen)
		median_eigen_top = median(mc_eigen_d$eigen[mc_eigen_d$metacell_annotation==mc_eigen_d$metacell_annotation[1]])
		module_annot[moi,"is_general"] = as.numeric((median_eigen_top > median_eigen_all * 1.2) > 1)
		
		
	}
	dev.off()
	
	# add top TFs
	module_genes = read.table(sprintf("%s/gmod.%s.wgcna.memberships.csv", out_fn, spi), header = TRUE, sep = "\t")
	module_genes$is_tf [is.na(module_genes$is_tf)] = FALSE
	module_genes_tfs = module_genes[ module_genes$is_tf, ]
	module_genes_tfs_a = aggregate(name ~ module, module_genes_tfs, function(v) { paste(stringr::str_trunc(head(unique(v), 20),60), collapse = "; ") })
	module_genes_tfs_v = dic_from_vecs(module_genes_tfs_a[,1], module_genes_tfs_a[,2])
	module_annot$top_tfs = module_genes_tfs_v[module_annot$module]
	
	# expression tag
	module_annot$expression_tag = "ct-specific"
	module_annot$expression_tag [ module_annot$num_active_mcs_total > module_annot$num_active_mcs_in_top * 1.2 & module_annot$num_active_mcs_total > 6 ] = "ct-general"
	module_annot$expression_tag [ module_annot$is_general == 1 & module_annot$num_active_mcs_total > 1 ] = "ct-general"
	module_annot$expression_tag [ module_annot$is_transversal == 1 ] = "transversal"
	module_annot$expression_tag [ module_annot$is_transversal == 1 & grepl("Foxj1", module_annot$top_tfs)] = "ciliary"
	module_annot$expression_tag [ module_annot$is_transversal == 1 & grepl("Myc", module_annot$top_tfs)] = "myc+"
	
	# sort and write
	module_annot$annotation_top = factor(module_annot$annotation_top, levels = unique(ctt$cell_type))
	module_annot = module_annot [ order(module_annot$annotation_top, -module_annot$num_active_mcs_in_top, -as.numeric(module_annot$frac_active_mcs_in_top)), ]
	module_annot_top_v = dic_from_vecs(module_annot$module, as.character(module_annot$annotation_top))	
	module_annot$module_name = sprintf("%s_%s | %s | %s", spi, module_annot$module, module_annot$annotation_top, module_annot$expression_tag)
	module_annot$module_color = ctt_cts_col_v [ module_annot$annotation_top ]
	module_annot_module_name_v = dic_from_vecs(module_annot$module, module_annot$module_name)
	
	# add orthology info to gene-per-module list
	mc_wgcna_gmods_annot$module_name = module_annot_module_name_v [ mc_wgcna_gmods_annot$module ]
	mc_wgcna_gmods_annot$orthogroup = oga_v [ mc_wgcna_gmods_annot$gene ]

	# concatenate...
	module_annot$module = paste(spi, module_annot$module, sep = "_")
	mc_wgcna_gmods_annot$module = paste(spi, mc_wgcna_gmods_annot$module, sep = "_")
	module_annot_t = rbind(module_annot_t, module_annot)
	mc_wgcna_gmods_annot_t = rbind(mc_wgcna_gmods_annot_t, mc_wgcna_gmods_annot)
	genes_per_gmod = table(mc_wgcna_gmods_annot$gene, mc_wgcna_gmods_annot$module)
	genes_per_gmod_jaccard = jaccard(genes_per_gmod, genes_per_gmod)
	write.table(module_annot, sprintf("%s/gmod.%s.wgcna.gmod_annotation.activity_per_mc.csv", out_fn, spi), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	write.table(mc_wgcna_gmods_annot, sprintf("%s/gmod.%s.wgcna.gmod_annotation.annot_table.csv", out_fn, spi), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
	
	# add barplots of eigenvalues, UMIfrac of contents, and FP of top TFs
	pdf(sprintf("%s/gmod.%s.gmod_eigengenes.barplots.pdf", out_fn, spi), width = 12, height = 32)
	for (moi in colnames(mc_wgcna_me)) {
		layout(1:16)
		# eigenvalue of module
		barplot(mc_wgcna_me[,moi], main = sprintf("%s_%s eigenvalues", spi, moi), space = 0, ylim = c(-0.1,0.5), col = ctm_mcs_col_v, las = 1, ylab = "eigenvalue", border = NA)

		# genes in module
		moi_genes = mc_wgcna_gmods_annot [ mc_wgcna_gmods_annot$module == sprintf("%s_%s", spi, gsub("^ME","",moi)) , "gene" ]
		mc_uf_i = mc_uf [ gsub("_","-", moi_genes), ]
		mc_uf_i_s = colSums(mc_uf_i)
		barplot(mc_uf_i_s, main = sprintf("%s_%s UMIfrac", spi, moi), space = 0, col = ctm_mcs_col_v, las = 1, ylab = "UMI/10^4", border = NA)
		
		# tfs in module
		moi_tfs = mc_wgcna_gmods_annot [ mc_wgcna_gmods_annot$module == sprintf("%s_%s", spi, gsub("^ME","",moi)) & mc_wgcna_gmods_annot$is_tf, "gene" ]
		for (tfi in moi_tfs) {
			vv = mc_fp[gsub("_","-",tfi),]
			vix = which(vv > 4)
			vv [ vix ] = 4
			b = barplot(
				vv,
				cex = 0.7,
				xlab = "",
				ylab = "FC",
				ylim = c(0,4),
				space = 0,
				border = NA,
				pch = 19, col = ctm_mcs_col_v[ctm$metacell], las = 2, cex.names = 0.6)
			text(b[vix], y = 4, "^", col = "gray10", cex = 0.5)
			title(main = sprintf("%s_%s\n%s (%s)\n%s", spi, moi, tfi, oga_v[tfi], ogm_gvo[tfi]), cex.main = 1, font.main = 1, col.sub = "gray10", col.main = "gray20")
			abline(h=1, lty = 2, col = "gray5")
		}
		
	}
	dev.off()

	
	# annots
	go_annot = gsa_topgo_load_emapper(emapper_fn = sprintf("../data/reference/%s_ensembl.GO.csv", spi_w), index_col_GOs = 2)
	pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi_w))
	names(pf_annot) = oga_gtv[names(pf_annot)]
	names(go_annot) = oga_gtv[names(go_annot)]

	# load KO info in species of interest
	keg_ll = read.table(sprintf("results_pathways/kegg_map_to_%s.csv",  spi), header = TRUE, sep = "\t")
	
	# load KO-to-KEGG pathway mappings
	kg2name_d = read.table("../data/reference/kegg_pathway2name.tsv", sep = "\t", header = FALSE, col.names = c("kegg", "kegg_name"))
	kg2name_v = dic_from_vecs(kg2name_d$kegg, kg2name_d$kegg_name)

	ko2kg_l = gsa_enrichment_load_pfam_list("../data/reference/kegg_pathway2ko.tsv", architecture_sep = ",")
	ko2kg_d = data.frame(
		kegg = unlist(sapply(1:length(ko2kg_l), function(nn) rep(names(ko2kg_l)[nn], length(ko2kg_l[[nn]])))),
		ko = unlist(ko2kg_l))
	keg_ld = merge(keg_ll[,c("gene","name","ko")], ko2kg_d, by.x = "ko", by.y = "ko", all.x = TRUE, all.y = FALSE)
	keg_ld = unique(keg_ld)
	keg_ld = keg_ld [ order(keg_ld$gene, keg_ld$kegg, keg_ld$ko), ]
	keg_ld$kegg_name = sprintf("%s: %s", keg_ld$kegg, stringr::str_trunc(kg2name_v[keg_ld$kegg], 40))
	keg_ld = keg_ld [ !is.na(keg_ld$kegg), ]
	keg_ldl = aggregate(kegg_name ~ gene, data = keg_ld, function(vv) { paste(vv, collapse = ",") })
	write.table(keg_ldl, sprintf("results_pathways/kegg_gene_to_map_%s.csv", spi), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	keg_lll = gsa_enrichment_load_pfam_list(sprintf("results_pathways/kegg_gene_to_map_%s.csv", spi), architecture_sep = ",")
	names(keg_lll) = gg_v [ names(keg_lll) ]


	dir.create(sprintf("%s/tests_gmod_nov24/", out_fn), showWarnings = FALSE)
	for (nn in 1:length(mc_wgcna_gmods)) {
		
		# foreground and backgound
		mki_module = names(mc_wgcna_gmods)[nn]
		mki_genes = unique(gg_v[mc_wgcna_gmods[[nn]]])
		mki_genes = as.character(mki_genes[ !is.na(mki_genes) ] )
		bgi_genes = unique(gg_v [unlist(mc_wgcna_gmods)])
		bgi_genes = as.character(bgi_genes [ !is.na(bgi_genes) ])
		
		message(sprintf("functional enrichments | %s module, n=%i genes", mki_module, length(mki_genes)))
		if (length(mki_genes)>0) {

			# # KEGG enrichment
			# kegg_enrichment = gsa_enrichment_hypergeometric_v2(
			# 	annotation = keg_lll,
			# 	genes_fg = mki_genes,
			# 	genes_bg = bgi_genes [ bgi_genes %in% names(keg_lll) ],
			# 	out_fn = sprintf("results_metacell_%s_filt/tests_gmod_nov24/test.%s.%s.kegg.csv", spi, spi, make.names(mki_module))
			# )
			# keg_ld_i = keg_ld [ keg_ld$gene %in% mki_genes, ]
			# keg_ld_i$kegg_name = factor(keg_ld_i$kegg_name, levels = kegg_enrichment[[1]]$annot)
			# keg_ld_i = keg_ld_i [ order(keg_ld_i$kegg_name, keg_ld_i$ko, keg_ld_i$gene), ]
			
			# # KEGG enrichment pvalue
			# kegg_enrichment_pv_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$pval_adj)			
			# kegg_enrichment_fg_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$freq_in_fg)	
			# keg_ld_i$kegg_pvalue = kegg_enrichment_pv_v[keg_ld_i$kegg_name]		
			# keg_ld_i$kegg_n_in_fg = kegg_enrichment_fg_v[keg_ld_i$kegg_name]		
			# write.table(keg_ld_i, sprintf("results_metacell_%s_filt/tests_gmod_nov24/test.%s.%s.keggsum.csv", spi, spi, make.names(mki_module)), sep = "\t", quote = FALSE, row.names = FALSE)
			
			# test domain enrichment
			message(sprintf("pathways %s | %s test pfam...", spi, mki_module))
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = mki_genes,
				genes_bg = bgi_genes [ bgi_genes %in% names(pf_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_gmod_nov24/test.%s.%s.pfam.csv", spi, spi, make.names(mki_module)),
				name_fg = make.names(mki_module)
			)

			# test GOs
			message(sprintf("pathways %s | %s test GOs...", spi, mki_module))
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = mki_genes,
				genes_bg = bgi_genes [ bgi_genes %in% names(go_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_gmod_nov24/test.%s.%s.topgo.csv", spi, spi, make.names(mki_module)),
				name_fg = make.names(mki_module),
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)
			
			# integrate GOs and pfam in a single table
			message(sprintf("pathways %s | merge GO+pfam...", spi, mki_module))
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("results_metacell_%s_filt/tests_gmod_nov24/test.%s.%s.summary.csv", spi, spi, make.names(mki_module)),
				extra_annot_dict = ogm_gv
			) 


		}

	}
	
}


message("done!")
