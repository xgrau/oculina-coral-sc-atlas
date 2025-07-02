#### Input ####

# libraries
source("../scripts/Downstream_functions.R")
source("../scripts/Cross_species_functions.R")
source("../scripts/helper.R")
library("umap")
library("phangorn")
library("ape")
graphics.off()

# input
out_fn = "results_csps_cell_types/"
dir.create(out_fn, showWarnings = FALSE)

# species
# data index
list_comparisons = list(
	alpoculi = list(c("Ocupat","Ocuarb"), ref = "Ocupat"),
	alpscler = list(c("Ocupat","Ocuarb","Amil","Spin"), ref = "Ocupat")
)

# lateral modules
lateral_modules = list(
	"Ocupat" = c(ciliary = "magenta", stress = c("maroon","white","bisque4"), proliferation = c("lightgreen", "pink"))
)

# orthogroups
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
ogm_gv_seu = dic_from_vecs(names = gsub("_", "-", ogm$gene), terms = ogm$orthogroup_name)


# fc threshold
fc_thr = 2

# fraction of variance required in PCA based dendrogram
# fraction_required = 0.4 # unset; determined automatically based on elbow rule!

# completeness
fraction_genes = 0.7

# identify icc markers
min_ec = 0

# bootstrap
num_bs = 100
collapse_bs = num_bs * 0.2

focus = "cell_type"
focid = "cts"

for (nn in 1:length(list_comparisons)) {

	# set id
	set_id = names(list_comparisons)[nn]
	
	sps_list = list_comparisons[[nn]][[1]]
	sps_ref = list_comparisons[[nn]][["ref"]]
	sps_query = sps_list [ !sps_list %in% sps_ref ]

	# load cell type data for all species
	ann_cts = c()
	for (spi in sps_list) {
		ctt_spi_fn = sprintf("results_metacell_%s_filt/annot.%s.leiden.csv", spi, spi)
		ctt_spi = read.table(ctt_spi_fn, header = TRUE, comment.char = "", sep = "\t")
		# get color annotations for each species
		ann_spi_v = dic_from_vecs(sprintf("%s|%s", spi, ctt_spi$cell_type), ctt_spi$color)
		ann_cts = c(ann_cts, ann_spi_v)
	}

	# loop query species, load icc ortholog pairs
	for (spi in sps_query) {

		# load cross-species marker data (from `s01`)
		message(sprintf("csps %s | load %s-%s ICC data...", set_id, sps_ref, spi))
		icc_genes_fn = sprintf("results_metacell_%s_filt/csps/dat.icc.%s-%s.ec_scores.csv", sps_ref, sps_ref, spi)
		if (file.exists(icc_genes_fn)) {
			icc_genes_i = read.table(icc_genes_fn, header = TRUE)
			icc_genes_i = icc_genes_i [ icc_genes_i$ec_value >= min_ec, ]
		} else {
			message(sprintf("MISSING! -> %s", icc_genes_fn))
		}
		
		if (spi == sps_query[1]) {
			icc_genes_f = icc_genes_i[,c(1,2)]
		} else {
			icc_genes_f = merge(icc_genes_f, icc_genes_i[,c(1,2)], by = "sp1", all.x = FALSE, all.y = FALSE)
		}

	}
	colnames(icc_genes_f) = c(sps_ref, sps_query)

	# loop species, load mcfp data
	for (spi in sps_list) {

		# load cross-species marker data (from `s01`)
		message(sprintf("csps %s | load %s expression data...", set_id, spi))
		fps_genes_fn = sprintf("results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi, spi)
		if (file.exists(fps_genes_fn)) {
			fps_i = readRDS(fps_genes_fn)
			fps_i_f = quantile_normalisation(fps_i)
			# fps_i_f = fps_i
			fps_i_f = fps_i_f [ icc_genes_f[,spi] , ]
			colnames(fps_i_f) = sprintf("%s|%s", spi, colnames(fps_i_f))
		} else {
			message(sprintf("MISSING! -> %s", fps_genes_fn))
		}
		
		if (spi == sps_list[1]) {
			csps_m = fps_i_f
		} else {
			csps_m = cbind(csps_m, fps_i_f)
		}

	}

	# drop trans
	csps_m = csps_m [ ,!grepl("unknown",colnames(csps_m)) ]

	# remove markers with an excess of absences (keep genes present in X% of tissues)
	csps_m = csps_m [ rowSums(!is.na(csps_m)) >= ncol(csps_m) * fraction_genes , ]

	# remove markers from lateral modules (stress, ciliary, translation/transcription-related modules)
	gmo = read.table(sprintf("results_metacell_%s_filt/gmod.%s.wgcna.gmod_annotation.annot_table.csv", sps_ref, sps_ref), sep = "\t", header = TRUE)
	gmo$moduleid = gsub(sprintf("^%s_", sps_ref), "", gmo$module)
	latmod_genes_v = c()
	for (latmoi in lateral_modules[[sps_ref]]) { 
		latmod_genes_v = c(latmod_genes_v, gmo [ gmo$moduleid == latmoi, "gene" ])
	}
	latmod_genes_v = gsub("_","-", latmod_genes_v)
	csps_m = csps_m [ !rownames(csps_m) %in% latmod_genes_v , ]

	# na to zero
	csps_m [ is.na(csps_m) ] = 0

	# # restrict to tfs?	
	# list_tfs = unique(read.table(sprintf("../results_annotation/gene_annotations/tfs.%s_genes.txt", sps_ref))[,1])
	# list_tfs = dictionary_t2g(sprintf("../data/reference/%s_long.annot.gtf", sps_ref), list_tfs)
	# list_tfs = list_tfs [ list_tfs %in% rownames(csps_m) ]
	# csps_m = csps_m [ list_tfs, ]

	# footprint discretisation
	message(sprintf("csps %s | %s discretise fps at fc>=%.2f...", set_id, focid, fc_thr))

	# binarise
	csps_m_b = (csps_m >= fc_thr) * 1

	# trinarise?
	# csps_m_t = csps_m
	# csps_m_t [ csps_m <= 1.0                ] = "A"
	# csps_m_t [ csps_m > 1 & csps_m < fc_thr ] = "C"
	# csps_m_t [ csps_m >= fc_thr             ] = "G"
	# csps_m_t [ csps_m >= fc_thr * 2         ] = "T"

	# drop non-enriched markers?
	# csps_m_b = csps_m_b [ apply(csps_m_b, 1, function(r) max(r) == max(csps_m_b)), ]


	# correlation
	message(sprintf("csps | weighted correlation %s...", set_id))
	icc_ecv = read.table(sprintf("results_metacell_%s_filt/csps/dat.icc.%s-%s.ec_scores.csv", sps_ref, sps_ref, sps_list[length(sps_list)]), header = TRUE, sep = "\t")
	icc_ecv_v = dic_from_vecs(icc_ecv$sp1, icc_ecv$ec_value)
	com = WGCNA::cor(csps_m, method = "pearson", weights.x = icc_ecv_v[rownames(csps_m)])

	# load annotations from both species
	ctt = data.frame()
	cmd = data.frame()
	for (spi in sps_list) {
		ctt_i = read.table(sprintf("results_metacell_%s_filt/annot.%s.leiden.csv", spi, spi), sep = "\t", header = TRUE, comment.char = "")
		ctt_i$cell_type_sps = paste(spi, make.names(ctt_i$cell_type), sep = "|")
		ctt = rbind(ctt, ctt_i)
		cmd_i = read.table(sprintf("results_metacell_%s_filt/dat.%s.cell_metadata.csv", spi, spi), sep = "\t", header = TRUE, comment.char = "")
		cmd_i$cell_type_sps = paste(spi, make.names(cmd_i$cell_type), sep = "|")
		cmd = rbind(cmd, cmd_i[,c("cell_type","cell_type_sps","color")])
	}
	dic_cts_col = dic_from_vecs(ctt$cell_type_sps, ctt$color)

	# print, original order	
	message(sprintf("csps | heatmap %s, original...", set_id))
	hm1 = plot_complex_heatmap(
		com,
		name = "wpearson",
		color_min = 0,
		color_max = 0.8,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		do_dotplot = FALSE,
		use_raster = FALSE,
		categories_row = rownames(com),
		categories_col = colnames(com),
		colors_row = dic_cts_col,
		colors_col = dic_cts_col,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		cex_dotplot = 0.03
	)
	# print, original order	
	message(sprintf("csps | heatmap %s, clustering...", set_id))
	hm2 = plot_complex_heatmap(
		com,
		name = "wpearson",
		color_min = 0,
		color_max = 0.8,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = TRUE,
		cluster_col = TRUE,
		do_dotplot = FALSE,
		use_raster = FALSE,
		categories_row = rownames(com),
		categories_col = colnames(com),
		colors_row = dic_cts_col,
		colors_col = dic_cts_col,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		cex_dotplot = 0.03
	)
	pdf(sprintf("%s/csps.%s.wpearson.hm.%s.pdf", out_fn, set_id, focid), width = 8+dim(com)[2]/20, height = 8+dim(com)[1]/20)
	print(hm1)
	print(hm2)
	dev.off()
	
	# sorted barplots
	message(sprintf("csps | sorted barplots %s...", set_id))
	pdf(sprintf("%s/csps.%s.wpearson.bp.%s.pdf", out_fn, set_id, focid), width = 8, height = 4)
	par(mar = c(8.1, 4.1, 4.1, 2.1))
	for (nnn in 1:nrow(com)) {
		
		sps_nn = gsub("\\|.*", "", rownames(com)[nnn])
		com_dd_f_v = com[nnn,]
		com_dd_f_v = com_dd_f_v [ 1:length(com_dd_f_v) != nnn ]
		# com_dd_f_v = com_dd_f_v [ !grepl(sps_nn, names(com_dd_f_v)) ]
		com_dd_f_v = sort(com_dd_f_v, decreasing = TRUE)
		com_dd_f_v = com_dd_f_v [ com_dd_f_v >= 0.05 ]
		b = barplot(com_dd_f_v, col = dic_cts_col [ names(com_dd_f_v) ] , las = 2, cex.names = 0.5, ylim = c(0,1), space = 0, ylab = "wpearson", xlim = c(0,40))
		text(b, com_dd_f_v, sprintf("%.2f", com_dd_f_v), srt = 90, cex = 0.7, col = "gray5")
		title(main = sprintf("%s\n%s", rownames(com)[nnn], focus), adj = 0, cex.main = 1, font.main = 1)
		
	}
	dev.off()

	message(sprintf("csps %s | %s n=%i markers in %i %ss...", set_id, focid, nrow(csps_m_b), ncol(csps_m_b), focus))


	# reduce
	pdf(sprintf("%s/csps.%s.dimred.%s.pdf", out_fn, set_id, focid))
	# umap
	message(sprintf("csps %s | %s reduce UMAP...", set_id, focid))
	csps_u = umap::umap(t(csps_m_b))
	plot(csps_u$layout, col = ann_cts [ rownames(csps_u$layout) ], pch = c(15,17,18,19,1,2,5,6) [ factor(gsub("\\|.*","", rownames(csps_u$layout))) ] )
	title(main = "umap")
	text(csps_u$layout, rownames(csps_u$layout), col = scales::alpha(ann_cts [ rownames(csps_u$layout) ], 0.6), cex = 0.4)

	# pca
	message(sprintf("csps %s | %s reduce PCA...", set_id, focid))
	if (ncol(csps_m_b) < nrow(csps_m_b)) {

		csps_p = princomp(csps_m_b)
		plot(csps_p$loadings[,c(1,2)], col = ann_cts [ rownames(csps_p$loadings) ], pch = c(15,17,18,19,1,2,5,6) [ factor(gsub("\\|.*","", rownames(csps_p$loadings))) ] )
		title(main = "pca")
		text(csps_p$loadings[,c(1,2)], rownames(csps_p$loadings), col = scales::alpha(ann_cts [ rownames(csps_p$loadings) ], 0.6), cex = 0.4)

		# top?
		pov = csps_p$sdev^2/sum(csps_p$sdev^2)
		
		num_pcs_v = find_pca_elbow(pov)
		n_pcs = num_pcs_v["Perpendicular line"]
		fraction_required = cumsum(pov) [ n_pcs ]
		
		plot(pov, col = "blue", main = "fraction variance", sub = sprintf(">=%.1fpp variance explained at PC %i", fraction_required * 100, n_pcs))
		abline(v = n_pcs, lty = 2, col = "red")

	}

	dev.off()


	# plot height
	plot_height = ceiling(ncol(csps_m_b) / 6 + 4)

	# dendrogram based on first X PCs
	if (ncol(csps_m_b) < nrow(csps_m_b)) {
		
		message(sprintf("csps %s | %s reduce PCA, keep n=%i PCs...", set_id, focid, n_pcs))
		csps_p_f = csps_p$loadings[,1:n_pcs]

		message(sprintf("csps %s | %s plot PCA dendrogram...", set_id, focid))
		csps_p_f_h = stats::hclust(dist(csps_p_f, method = "manhattan"), method = "average")
		csps_p_f_t = as.phylo(csps_p_f_h)

		# root?
		ali_tips = csps_p_f_t$tip.label
		if (any(grepl("progenitors", ali_tips))) {
			ali_tips = ali_tips [ grepl("progenitor", ali_tips) ]
			ali_tips = ali_tips [ !grepl("Tadh", ali_tips) ]
			ali_root = ape::getMRCA(csps_p_f_t, ali_tips)
			csps_p_f_t = ape::root(csps_p_f_t, node = ali_root)
		}
		
		pdf(sprintf("%s/csps.%s.dendrogram.%s.PCA.pdf", out_fn, set_id, focid), height = plot_height, width = 12)
		ape::plot.phylo(ladderize(csps_p_f_t), tip.color = ann_cts [ csps_p_f_t$tip.label ], font = 1, underscore = TRUE, type = "phylo", main = sprintf("upgma from %i PCs", n_pcs))
		ape::add.scale.bar()
		dev.off()
		ape::write.tree(ladderize(csps_p_f_t), sprintf("%s/csps.%s.dendrogram.%s.PCA.newick", out_fn, set_id, focid))
		
	}

	# maximum likelihood trees
	# binary matrix
	ali_f = phangorn::as.phyDat(t(csps_m_b), type = "USER", levels = names(table(csps_m_b)), names = rownames(csps_m_b), return.index = TRUE)

	# distance and initial parsimony tree
	message(sprintf("csps %s | %s parsimony tree...", set_id, focid))
	parsimony_fun =  function(x) { reorder(ape::as.phylo(fastcluster::hclust(phangorn::dist.logDet(x), method = "average")), "postorder") }
	ali_f_phy = parsimony_fun(ali_f)
	
	# root at progenitors
	ali_tips = ali_f_phy$tip.label
	if (any(grepl("progenitors", ali_tips))) {
		ali_tips = ali_tips [ grepl("neurosecretory_progenitors", ali_tips) ]
		ali_root = ape::getMRCA(ali_f_phy, ali_tips)
		ali_f_phy = ape::root(ali_f_phy, node = ali_root, resolve.root = TRUE)
	}
	
	# avoid inf lengths
	ali_f_phy$edge.length [ is.infinite(ali_f_phy$edge.length) ] = max(ali_f_phy$edge.length [ !is.infinite(ali_f_phy$edge.length) ]) * 1.2

	# bootstrap in UPGMA (may fail, retry often)
	boo_iter = 1
	boo_total = num_bs
	boo_list = NULL
	
	message(sprintf("csps %s | %s parsimony tree, bootstrap...", set_id, focid))
	boo_list = phangorn::bootstrap.phyDat(ali_f, parsimony_fun, multicore = TRUE, bs = 100)
	
	# collapse poorly supported nodes
	collapse_nodes = which(prop.clades(ali_f_phy, boo_list) < collapse_bs) + length(ali_f_phy$tip.label)
	rownames(ali_f_phy$edge) = as.character(1:nrow(ali_f_phy$edge))
	collapse_edges = sort(as.numeric(rownames(ali_f_phy$edge) [ ali_f_phy$edge[,2] %in% collapse_nodes ]))
	ali_f_phy_col = ali_f_phy
	ali_f_phy_col$edge.length[ collapse_edges ] = 0
	ali_f_phy_col = ape::di2multi(ali_f_phy_col)

	# add node names
	ali_f_phy = ape::makeNodeLabel(ali_f_phy)

	message(sprintf("csps %s | %s plot...", set_id, focid))

	pdf(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.pdf", out_fn, set_id, focid), height = plot_height, width = 12)
	phangorn::plotBS(ali_f_phy, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP")
	ape::add.scale.bar()
	phangorn::plotBS(ali_f_phy_col, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with FBP", fc_thr), root.edge = TRUE, method = "FBP", use.edge.length = FALSE)
	phangorn::plotBS(ali_f_phy, boo_list, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f with TBE", fc_thr), root.edge = TRUE, method = "TBE", digits = 1)
	ape::add.scale.bar()
	ape::plot.phylo(ali_f_phy, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, underscore = TRUE, type = "phylo", main = sprintf("UPGMA binarised at fp>%.2f", fc_thr), root.edge = TRUE, use.edge.length = TRUE, show.node.label = TRUE)
	# ape::edgelabels()
	# ape::nodelabels()
	ape::add.scale.bar()
	dev.off()
	ape::write.tree(ali_f_phy, sprintf("%s/csps.%s.dendrogram.%s.UPGMA.newick", out_fn, set_id, focid))


	## Ancestral reconstruction ##

	# loop species, load markers, as table
	list_markers_tfs = data.frame()
	list_markers_tfs_v = c()
	for (spi in sps_list) {

		# load cross-species marker data (from `s01`)
		message(sprintf("csps %s | load %s markers...", set_id, spi))
		mks_genes_fn = sprintf("results_metacell_%s_filt/seu.%s.markers.cts.csv", spi, spi)
		mks_i = read.table(mks_genes_fn, sep = "\t", header = TRUE)
		mks_i = mks_i [ !is.na(mks_i$gene) & mks_i$p_val_adj < 0.01 & mks_i$avg_log2FC > 0, ]
		mks_i$orthogroup = oga_gv [ mks_i$gene ]
		mks_i$orthogroup_name = ogm_gv [ mks_i$gene ]
		
		focus_nodes = unique(mks_i$focus_node)
		for (foi in focus_nodes) {
			mks_ii = mks_i [ mks_i$focus_node == foi & mks_i$is_tf, c("gene", "orthogroup") ]
			mks_ii = unique(mks_ii)
			# mks_ii = mks_ii [ mks_ii != "" ]
			# mks_ii = mks_ii [ !is.na(mks_ii) ]
			if (nrow(mks_ii) > 0) {
				list_markers_tfs = rbind(list_markers_tfs, data.frame(focus_node = sprintf("%s|%s", spi, foi), orthogroup = mks_ii$orthogroup, gene = mks_ii$gene))
			} else {
				list_markers_tfs = rbind(list_markers_tfs, data.frame(focus_node = sprintf("%s|%s", spi, foi), orthogroup = NA, gene = NA))
			}
		}

	}
	list_markers_tfs_m = lapply(apply(combn(names(list_markers_tfs[,c(1,2)]), 2), 2, function(i) list_markers_tfs[,c(1,2)][i]), table)[[1]]
	list_markers_tfs_names_d = data.frame(gene = list_markers_tfs$gene, orthogroup = list_markers_tfs$orthogroup)
	list_markers_tfs_names_d = unique(list_markers_tfs_names_d)
	list_markers_tfs_names_d = list_markers_tfs_names_d [ order(list_markers_tfs_names_d$orthogroup, list_markers_tfs_names_d$gene), ]
	list_markers_tfs_names_d$orthogroup_name = ogm_gv [ list_markers_tfs_names_d$gene ]
	list_markers_tfs_names_d = list_markers_tfs_names_d [ !is.na(list_markers_tfs_names_d$orthogroup_name), ]
	list_markers_tfs_names_v = dic_from_vecs(list_markers_tfs_names_d$orthogroup, list_markers_tfs_names_d$orthogroup_name)

	# wrangle and select markers to reconstruct
	message(sprintf("csps %s | %s init ancestral reconstruction...", set_id, focid))
	ali_f_phy_sansnodes = ali_f_phy
	ali_f_phy_sansnodes$node.label = NULL
	list_markers_tfs_m_f = list_markers_tfs_m [ , apply(list_markers_tfs_m, 2, function(vv) sum( vv > 0 ) >= 2 ) ]
	# rownames(list_markers_tfs_m_f) = make.names(rownames(list_markers_tfs_m_f))
	rownames(list_markers_tfs_m_f) = gsub("/", ".", rownames(list_markers_tfs_m_f))

	# descendant nodes
	descendants_from_nodes = adephylo::listTips(ali_f_phy)
	list_labels = c(ali_f_phy$tip.label, ali_f_phy$node.label)
	nodes_to_test = list_labels[!list_labels %in% ali_f_phy$node.label[1]]

	# open pdf and loop
	message(sprintf("csps %s | %s loop over n=%i markers...", set_id, focid, ncol(list_markers_tfs_m_f)))
	pdf(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.markers.pdf", out_fn, set_id, focid), height = plot_height, width = 24)
	# init tables
	markers_per_node = data.frame()
	for (nnn in 1:ncol(list_markers_tfs_m_f)) {
		
		ggo = colnames(list_markers_tfs_m_f)[nnn]
		vv_i = factor(as.character(list_markers_tfs_m_f[ali_f_phy$tip.label, ggo]), levels = c("0","1"))
		ace_i = ape::ace(vv_i, ali_f_phy_sansnodes, type = "discrete", method = "ML", marginal = TRUE, model = "ER")
		ace_i_d = ace_i$lik.anc
		rownames(ace_i_d) = ali_f_phy$node.label
		vv_i_m = data.frame(row.names = ali_f_phy$tip.label, "0" = as.numeric(vv_i == 0), "1" = as.numeric(vv_i == 1))
		
		ape::plot.phylo(ali_f_phy, tip.color = ann_cts [ ali_f_phy$tip.label ], font = 1, underscore = TRUE, type = "phylo", root.edge = TRUE, use.edge.length = TRUE, show.node.label = TRUE)
		ape::nodelabels(node = 1:ali_f_phy$Nnode+Ntip(ali_f_phy), pie = ace_i_d, piecol = c("honeydew2","blue3"), cex = 0.2)
		ape::tiplabels(pie = vv_i_m, piecol = c("honeydew2","blue3"), cex = 0.1)
		title(main = sprintf("%s\n%s", ggo, list_markers_tfs_names_v[ggo]))
		
		# ace_i_d_f = ace_i_d [ ace_i_d[,"1"] >= 0.9, ]
		# save table
		ace_i_d_i = data.frame(
			orthogroup = ggo,
			orthogroup_name = list_markers_tfs_names_v [ ggo ],
			probability = ace_i_d[,"1"],
			present = ace_i_d[,"1"] >= 0.7,
			node = rownames(ace_i_d),
			tips_in_node = unlist(lapply(descendants_from_nodes[rownames(ace_i_d)], function(vv) paste(sort(names(vv)), collapse = ",")  ))
		)
		markers_per_node = rbind(markers_per_node, ace_i_d_i)
		
	}
	dev.off()

	# log
	message(sprintf("csps %s | %s ancestral reconstruction complete...", set_id, focid))

	# sort and out
	markers_per_node = markers_per_node [ order(markers_per_node$node, -markers_per_node$probability, markers_per_node$orthogroup_name) ,]
	write.table(markers_per_node, sprintf("%s/csps.%s.dendrogram.%s.UPGMA.markers.csv", out_fn, set_id, focid), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

	message(sprintf("csps %s | %s trees done!", set_id, focid))



	## Plot TF expression in ordered heatmap ##

	ali_f_phy_is_tip = ali_f_phy$edge[,2] <= length(ali_f_phy$tip.label)
	ali_f_phy_ordered_tips = ali_f_phy$tip.label[ali_f_phy$edge[ali_f_phy_is_tip, 2]]

	ogs_for_plot = unique(markers_per_node$orthogroup)
	genes_for_plot = intersect(gsub("_", "-", list_markers_tfs$gene), rownames(csps_m))
	csps_m_order_f = csps_m [ genes_for_plot , ali_f_phy_ordered_tips ]
	
	# top per metacell
	expr_m = csps_m_order_f 
	top_markers = select_top_markers(expr_m, matrix_thr = 1.5, n_top_markers = 1e6, n_markers_rollmean = 4)
	expr_m = expr_m [ top_markers, ]
	top_cl_per_gene = colnames(expr_m) [ apply(expr_m, 1, which.max) ]
	rownames(expr_m) = sprintf("%s | %s", rownames(expr_m), ogm_gv_seu[rownames(expr_m)]) 
	rownames(expr_m) = stringr::str_trunc(rownames(expr_m), 100)
	pp1 = plot_complex_heatmap(
		expr_m,
		name = "Footprint",
		cluster_row = FALSE,
		cluster_col = FALSE,
		use_raster = FALSE,
		color_min = 1,
		color_max = 4,
		color_mat = c("gray95","orange","orangered2","#520c52"),
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cl_per_gene,
		categories_col = colnames(expr_m),
		colors_row = ann_cts,
		colors_col = ann_cts,
		fontsize = 6
	)
	pdfheight = max(8, (nrow(expr_m)/8 + 6))
	pp1@matrix_param$width  = unit(ncol(expr_m)/10, "in")
	pp1@matrix_param$height = unit(nrow(expr_m)/10, "in")
	pdf(sprintf("%s/csps.%s.dendrogram.%s.UPGMA.ordered_heatmap_TFs.pdf", out_fn, set_id, focid), width = 20, height = pdfheight)
	print(pp1)
	dev.off()

	## Save matrix ##
	saveRDS(csps_m, sprintf("%s/csps.%s.cspsmatrix.%s.rds", out_fn, set_id, focid))

}