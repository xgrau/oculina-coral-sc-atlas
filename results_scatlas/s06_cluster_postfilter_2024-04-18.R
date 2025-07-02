# libraries
suppressMessages(source("../scripts/Seurat_functions.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))
graphics.off()


# data index
sps_list = c("Ocupat","Ocuarb","Spin","Amil","Nvec","Xesp")

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
ogm_tv = dic_from_vecs(names = ogm$transcript, terms = ogm$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)


for (spi in sps_list) {
	
	# set working species
	if (spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}

	# output
	out_fn = sprintf("results_metacell_%s_filt2/", spi)
	dir.create(out_fn, showWarnings = FALSE)

	## Load gene annotations ##

	# load gene annotations for this species
	gna_i = ogm [ ogm$species == spi_w, ]
	gna_p = read.table(sprintf("../data/reference/%s_long.pep.annotations.csv", spi_w), sep = "\t", col.names = c("transcript","blast","pfam"))
	gna_p$gene = dictionary_t2g(gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi_w), vector_to_fix = gna_p$transcript)
	gna_i = merge(gna_i, gna_p[,c("gene","pfam")], by.x = "gene", by.y = "gene", all.x = TRUE, all.y = TRUE)
	gene_annot = data.frame(gene = gna_i$gene, name = gna_i$orthogroup_name, pfam = gna_i$pfam)
	gene_annot = gene_annot [ !duplicated(gene_annot$gene), ]
	rownames(gene_annot) = gene_annot$gene
	pfa_v = dic_from_vecs(gene_annot$gene, unlist(lapply(gene_annot$pfam, function(v) { paste(unique(unlist(stringr::str_split(v, "/"))), collapse = "/") })))
	
	# load TF annotations
	gene_annot_tfs = ogm [ grepl(":", ogm$orthogroup_name) & ogm$species == spi_w, ]
	gene_annot_tfs_v = dic_from_vecs(gene_annot_tfs$gene, gene_annot_tfs$orthogroup_name)
	gene_annot$name [ rownames(gene_annot) %in% gene_annot_tfs$gene ] = paste(gene_annot_tfs_v [ rownames(gene_annot) [ rownames(gene_annot) %in% gene_annot_tfs$gene ] ] , gene_annot$name [ rownames(gene_annot) %in% gene_annot_tfs$gene ], sep = ":")
	gene_annot$is_tf = gene_annot$gene %in% gene_annot_tfs$gene
	
	# seurat gene dict
	gg_v = dic_from_vecs(gsub("_","-", gna_p$gene), gna_p$gene)
	
	# gene names from spis, for reference
	if (spi != "Spis") {
		og_pairs = clean_og_pairs(
			og_pairs_fn = "../data/orthology_Anthozoa_plus/orthogroup_pairs.genes.csv", 
			sp1 = "Spis", 
			sp2 = spi, 
			t2g = FALSE,
			header = TRUE
		)
		og_pairs_v = dic_from_vecs(og_pairs$sp1, og_pairs$sp2)
		gnspis = read.table("../data/reference/Spis_Genes_ncbi_dataset.tsv", header = TRUE, sep = "\t")
		gnspis = gnspis [ !is.na(gnspis$Protein.accession), ]
		gnspis = gnspis [ gnspis$Protein.accession != "", ]
		gnspis$Name = gsub(" ", "_", gnspis$Name)
		gnspis$gene_focus_species = og_pairs_v [ gnspis$Protein.accession ]
		gnspis = gnspis [ !is.na(gnspis$gene_focus_species),  ]
		gnspis_v = dic_from_vecs(gnspis$gene_focus_species, gnspis$Name)
	} else {
		gnspis = read.table("../data/reference/Spis_Genes_ncbi_dataset.tsv", header = TRUE, sep = "\t")
		gnspis = gnspis [ !is.na(gnspis$Protein.accession), ]
		gnspis = gnspis [ gnspis$Protein.accession != "", ]
		gnspis$Name = gsub(" ", "_", gnspis$Name)
		gnspis_v = dic_from_vecs(gnspis$Protein.accession, gnspis$Name)
	}


	## Reload without altering ##
	
	# object
	message(sprintf("metacell | %s | load Seurat...", spi))
	seu = readRDS(sprintf("%s/dat.%s.seurat_final.rds", out_fn, spi))

	## Load Seurat ##

	message(sprintf("metacell | %s | load annotated clusters...", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt_leiden    = ctt$cluster [ ctt$cluster_class == "leiden" ]
	ctt_lou_cts_v = dic_from_vecs(ctt$cluster, ctt$cell_type)
	ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	
	# update metadata
	message(sprintf("metacell | %s | update metadata...", spi))
	seu@meta.data$cell_type = ctt_lou_cts_v [ as.character(seu$leiden) ]
	seu@meta.data$cell_type [ as.character(seu$metacell) %in% ctt$cluster ] = ctt_lou_cts_v  [ as.character(seu$metacell) ] [ as.character(seu$metacell) %in% ctt$cluster ]
	seu@meta.data$color = ctt_cts_col_v [ as.character(seu$cell_type) ]
	
	# cross-tabulation of metacells and leiden clusters
	message(sprintf("metacell | %s | cross-tabulate metacells and leiden...", spi))
	xtb = table(seu$metacell, seu$leiden)
	xtf_lei = xtb / rowSums(xtb)
	order_cols = names(sort(apply(xtf_lei,1,function(x) which.max(x))))
	xtf_lei = xtf_lei[order_cols,]
	top_cluster_per_mc = levels(seu$leiden) [ apply(xtf_lei, 1, which.max) ]
	pp = plot_complex_heatmap(
		xtf_lei,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cluster_per_mc,
		colors_row = ctt_lou_col_v,
		categories_col = colnames(xtf_lei),
		colors_col = ctt_lou_col_v
	)
	pdf(sprintf("%s/seu.%s.crosstab.lei-metacell.pdf", out_fn, spi), height = (dim(xtf_lei)[1] / 10) + 5, width = (dim(xtf_lei)[2] / 10) + 5)
	print(pp)
	dev.off()

	# cross-tabulation of metacells and cell types
	message(sprintf("metacell | %s | cross-tabulate metacells and cell types...", spi))
	xtb = table(seu$metacell, seu$cell_type)
	xtf_cts = xtb / rowSums(xtb)
	order_cols = names(sort(apply(xtf_cts,1,function(x) which.max(x))))
	xtf_cts = xtf_cts[order_cols,]
	top_cluster_per_mc = levels(seu$cell_type) [ apply(xtf_cts, 1, which.max) ]
	pp = plot_complex_heatmap(
		xtf_cts,
		color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
		cluster_row = FALSE,
		cluster_col = FALSE,
		cell_border = gpar(col = "white", lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cluster_per_mc,
		colors_row = ctt_cts_col_v,
		categories_col = colnames(xtf_cts),
		colors_col = ctt_cts_col_v
	)
	pdf(sprintf("%s/seu.%s.crosstab.cts-metacell.pdf", out_fn, spi), height = (dim(xtf_cts)[1] / 10) + 5, width = (dim(xtf_cts)[2] / 10) + 5)
	print(pp)
	dev.off()


	# setup metacell annotation table, leiden
	message(sprintf("metacell | %s | metacell annotation table, leiden...", spi))
	mc_annot_lei_l = apply(xtf_lei, 1, function(v) { 
		vs = sort(v, decreasing = TRUE)
		vs = vs [ vs > 0.05 ]
		vn = names(vs)
		vt = vn[1]
		vv = paste(sprintf("%s|%.2f", vn, vs), collapse = ";")
		return(list(top = vt, others = vv))
	})
	mc_annot_lei_top = unlist(lapply(mc_annot_lei_l, function(i) i$top))
	mc_annot_lei_ann = unlist(lapply(mc_annot_lei_l, function(i) i$others))

	# setup metacell annotation table, cell types
	message(sprintf("metacell | %s | metacell annotation table, cell types...", spi))
	mc_annot_cts_l = apply(xtf_cts, 1, function(v) { 
		vs = sort(v, decreasing = TRUE)
		vs = vs [ vs > 0.05 ]
		vn = names(vs)
		vt = vn[1]
		vv = paste(sprintf("%s|%.2f", vn, vs), collapse = ";")
		return(list(top = vt, others = vv))
	})
	mc_annot_cts_top = unlist(lapply(mc_annot_cts_l, function(i) i$top))
	mc_annot_cts_ann = unlist(lapply(mc_annot_cts_l, function(i) i$others))
	
	# create metacell annotation table
	ctm = data.frame(
		metacell = names(mc_annot_cts_top),
		metacell_annotation = mc_annot_cts_top,
		color = ctt_cts_col_v [ mc_annot_cts_top ],
		metacell_top_cell_type = mc_annot_cts_top [ names(mc_annot_cts_top) ],
		metacell_ann_cell_type = mc_annot_cts_ann [ names(mc_annot_cts_top) ],
		metacell_top_leiden    = mc_annot_lei_top [ names(mc_annot_cts_top) ],
		metacell_ann_leiden = mc_annot_lei_ann [ names(mc_annot_cts_top) ]
	)
	# order metacells: by most abundant cell type, first
	ctm$metacell_top_cell_type = factor(ctm$metacell_top_cell_type, levels = levels(ctt$cell_type))
	ctm$metacell_top_leiden = factor(ctm$metacell_top_leiden, levels = ctt_leiden)
	# order metacells: by metacell clustering, second
	SeuratObject::Idents(seu) = seu@meta.data[,"metacell"]
	seu_num_pcs_v = find_pca_elbow(seu@reductions$pca@stdev)
	seu_num_pcs = seu_num_pcs_v["First derivative"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca_integrated_harmony", dims = 1:seu_num_pcs)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	ctm_metacell_factor = factor(ctm$metacell, levels = seu_phylo$tip.label)
	# apply dual order
	ctm = ctm [ order(ctm$metacell_top_cell_type, ctm$metacell_top_leiden, ctm_metacell_factor), ]
	# save
	write.table(ctm, sprintf("%s/annot.%s.mcs.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)

	# reorder metacells
	seu@meta.data$metacell = factor(seu@meta.data$metacell, levels = ctm$metacell)


	## Dimensionality reduction ##
	
	message(sprintf("metacell | %s | recolor umap...", spi))

	# first, calculate medioids
	# cell type centroids
	ct_centroids_d = data.frame(
		umap1 = seu@reductions$umap_integrated@cell.embeddings[,1], 
		umap2 = seu@reductions$umap_integrated@cell.embeddings[,2], 
		label = seu@meta.data$cell_type)
	ct_centroids_a = aggregate(cbind(ct_centroids_d$umap1,ct_centroids_d$umap2) ~ label, ct_centroids_d, median)
	# leiden centroids
	le_centroids_d = data.frame(
		umap1 = seu@reductions$umap_integrated@cell.embeddings[,1], 
		umap2 = seu@reductions$umap_integrated@cell.embeddings[,2], 
		label = seu@meta.data$leiden)
	le_centroids_a = aggregate(cbind(le_centroids_d$umap1,le_centroids_d$umap2) ~ label, le_centroids_d, median)
	# metacell centroids
	mc_centroids_d = data.frame(
		umap1 = seu@reductions$umap_integrated@cell.embeddings[,1], 
		umap2 = seu@reductions$umap_integrated@cell.embeddings[,2], 
		label = seu@meta.data$metacell)
	mc_centroids_a = aggregate(cbind(mc_centroids_d$umap1,mc_centroids_d$umap2) ~ label, mc_centroids_d, median)

	# open umap
	pdf(sprintf("%s/seu.%s.umap.pdf", out_fn, spi), height = 12, width = 12)
	layout(matrix(1:4, nrow = 2, byrow = TRUE))

	# plot integrated, by cell type
	plot(
		seu@reductions$umap_integrated@cell.embeddings[,1], seu@reductions$umap_integrated@cell.embeddings[,2], 
		col = "seashell3", pch = 19, cex = 0.5, main = "integrated, by cell type", xlab = "UMAP1", ylab = "UMAP2")
	points(mc_centroids_a[,2], mc_centroids_a[,3], font = 2, col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a[,1])], 0.4), bg = ctm_mcs_col_v[as.character(mc_centroids_a[,1])], pch = 21)
	
	# legends
	plot(0,0, xaxt = "n", xlab = "", ylab = "", yaxt = "n", col = NA, axes = FALSE)
	legend("topleft", names(ctt_cts_col_v), col = ctt_cts_col_v, pch = 19, cex = 0.5, title = "cell type", bty = "n", ncol = 2)
	legend("bottomright", names(seu_lou_col_v), col = seu_lou_col_v, pch = 19, cex = 0.5, title = "cell type", bty = "n", ncol = 4)
	
	layout(1)
	# plot integrated, by metacell and painted by leiden
	plot(seu@reductions$umap_integrated@cell.embeddings[,1], seu@reductions$umap_integrated@cell.embeddings[,2],
		col = seu@meta.data$leiden_color, pch = 19, cex = 0.5, main = "integrated, by leiden", xlab = "UMAP1", ylab = "UMAP2")
	text(le_centroids_a[,2], le_centroids_a[,3], le_centroids_a[,1], font = 1, col = colorspace::darken(ctm_mcs_col_v[as.character(le_centroids_a[,1])], 0.4), cex = 0.4)
	text(le_centroids_a[,2], le_centroids_a[,3], le_centroids_a[,1], font = 2, col = scales::alpha("gray5", 0.7), cex = 0.7)
	title(sub = sprintf("n=%i cells into %i metacells", ncol(seu), nrow(ctm)))

	layout(1)
	# plot integrated, by metacell and painted by cell type
	plot(
		seu@reductions$umap_integrated@cell.embeddings[,1], seu@reductions$umap_integrated@cell.embeddings[,2],
		col = seu@meta.data$color, pch = 19, cex = 0.5, main = "integrated, by cell type", xlab = "UMAP1", ylab = "UMAP2")
	text(mc_centroids_a[,2], mc_centroids_a[,3], mc_centroids_a[,1], font = 1, col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a[,1])], 0.4), cex = 0.4)
	text(ct_centroids_a[,2], ct_centroids_a[,3], ct_centroids_a[,1], font = 2, col = scales::alpha("gray5", 0.7), cex = 0.7)
	title(sub = sprintf("n=%i cells into %i metacells", ncol(seu), nrow(ctm)))

	dev.off()
	
	
	## Expression tables ##
	
	message(sprintf("metacell | %s | get non-orphan...", spi))
	keep_genes = rownames(seu@assays$RNA$counts)
	keep_genes = keep_genes [ !grepl("orphan", keep_genes) ]
	
	message(sprintf("metacell | %s | footprints, metacell...", spi))
	mc_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "metacell", data = "RNA", keep_genes = keep_genes)
	mc_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$metacell, sum) }))
	mc_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$metacell, function(vv) { 100*sum(vv > 0)/length(vv) } ) } ))
	mc_umifrac = sca_mc_gene_umifrac_noobj(mc_counts, multiplying_factor = 10000)

	message(sprintf("metacell | %s | footprints, cell type...", spi))
	ct_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "cell_type", data = "RNA", keep_genes = keep_genes)
	ct_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$cell_type, sum) }))
	ct_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$cell_type, function(vv) { 100*sum(vv > 0)/length(vv) } ) } ))
	ct_umifrac = sca_mc_gene_umifrac_noobj(ct_counts, multiplying_factor = 10000)
	
	message(sprintf("metacell | %s | footprints, leiden...", spi))
	le_fp = sca_seurat_cell_type_fp(seu_object = seu, cluster = "leiden", data = "RNA", keep_genes = keep_genes)
	le_counts  = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$leiden, sum) }))
	le_fracell = t(apply(as.matrix(seu@assays$RNA$counts[keep_genes,]), 1, function(v) { tapply(v, seu$leiden, function(vv) { 100*sum(vv > 0)/length(vv) } ) } ))
	le_umifrac = sca_mc_gene_umifrac_noobj(le_counts, multiplying_factor = 10000)
	
	# save expression
	# footprints
	saveRDS(ct_fp, sprintf("%s/dat.%s.expression.cts_fp.rds", out_fn, spi))
	saveRDS(mc_fp, sprintf("%s/dat.%s.expression.mcs_fp.rds", out_fn, spi))
	saveRDS(le_fp, sprintf("%s/dat.%s.expression.lei_fp.rds", out_fn, spi))

	# umifrac
	saveRDS(ct_umifrac, sprintf("%s/dat.%s.expression.cts_umifrac.rds", out_fn, spi))
	saveRDS(mc_umifrac, sprintf("%s/dat.%s.expression.mcs_umifrac.rds", out_fn, spi))
	saveRDS(le_umifrac, sprintf("%s/dat.%s.expression.lei_umifrac.rds", out_fn, spi))
		
	# umi counts
	saveRDS(ct_counts, sprintf("%s/dat.%s.expression.cts_umicount.rds", out_fn, spi))
	saveRDS(mc_counts, sprintf("%s/dat.%s.expression.mcs_umicount.rds", out_fn, spi))
	saveRDS(le_counts, sprintf("%s/dat.%s.expression.lei_umicount.rds", out_fn, spi))
		
	# cellfrac
	saveRDS(ct_fracell, sprintf("%s/dat.%s.expression.cts_fracell.rds", out_fn, spi))
	saveRDS(mc_fracell, sprintf("%s/dat.%s.expression.mcs_fracell.rds", out_fn, spi))
	saveRDS(le_fracell, sprintf("%s/dat.%s.expression.lei_fracell.rds", out_fn, spi))
	
	
	## Find cluster markers ##

	# find markers for cell types
	SeuratObject::Idents(seu) = seu@meta.data[,"cell_type"]
	mks_cts = data.frame()
	for (loi in levels(seu@meta.data[,"cell_type"])) {
		message(sprintf("seurat | %s | markers cell_type %s...", spi, loi))
		mki = Seurat::FindMarkers(seu, ident.1 = loi, layer = "RNA", min.pct = 0.05, test.use = "wilcox")
		mki = mki [ !grepl("orphan", rownames(mki)), ]
		ct_umifrac_i_v = dic_from_vecs(gsub("_", "-", rownames(ct_umifrac)), ct_umifrac[,loi])
		mki$umifrac = sprintf("%.2f", ct_umifrac_i_v [ rownames(mki) ])
		mki$focus_node = loi
		mki$gene = rownames(mki)
		rownames(mki) = NULL
		mks_cts = rbind(mks_cts, mki)
	}
	# dictionary to reconstruct gene names
	mks_cts$gene [ !grepl("orphan", mks_cts$gene) ] = gg_v [ mks_cts$gene ] [ !grepl("orphan", mks_cts$gene) ]
	mks_cts$gene_name = ogm_gv [ mks_cts$gene ]
	mks_cts$is_tf = mks_cts$gene %in% gene_annot_tfs$gene
	mks_cts$pfam = pfa_v [ mks_cts$gene ]
	mks_cts$gene_name_Spis = gnspis_v [ mks_cts$gene ]
	write.table(mks_cts, sprintf("%s/seu.%s.markers.cts.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
	# find markers for metacells
	SeuratObject::Idents(seu) = seu@meta.data[,"metacell"]
	mks_mcs = data.frame()
	for (loi in sort(unique(seu@meta.data[,"metacell"]))) {
		message(sprintf("seurat | %s | markers metacell %s...", spi, loi))
		mki = Seurat::FindMarkers(seu, ident.1 = loi, assay = "RNA", slot = "counts", fc.slot = "counts", min.pct = 0.05, test.use = "wilcox")
		mki = mki [ !grepl("orphan", rownames(mki)), ]
		mc_umifrac_i_v = dic_from_vecs(gsub("_", "-", rownames(mc_umifrac)), mc_umifrac[,loi])
		mki$umifrac = sprintf("%.2f", mc_umifrac_i_v [ rownames(mki) ])
		mki$focus_node = loi
		mki$gene = rownames(mki)
		rownames(mki) = NULL
		mks_mcs = rbind(mks_mcs, head(mki, 50))
	}
	# dictionary to reconstruct gene names
	mks_mcs$gene [ !grepl("orphan", mks_mcs$gene) ] = gg_v [ mks_mcs$gene ] [ !grepl("orphan", mks_mcs$gene) ]
	mks_mcs$gene_name = ogm_gv [ mks_mcs$gene ]
	mks_mcs$is_tf = mks_mcs$gene %in% gene_annot_tfs$gene
	mks_mcs$pfam = pfa_v [ mks_mcs$gene ]
	mks_mcs$gene_name_Spis = gnspis_v [ mks_mcs$gene ]
	write.table(mks_mcs, sprintf("%s/seu.%s.markers.mcs.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	

	## Build cell type tree ##
	
	# tree of cell types
	SeuratObject::Idents(seu) = seu@meta.data[,"cell_type"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca_integrated_harmony", dims = 1:10)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	tip_color_dictionary = dic_from_vecs(seu$cell_type, seu$color)
	mks_bipartition_cts = sca_seurat_tree_markers(
		seu = seu, seu_phylo = seu_phylo,
		output_file = sprintf("%s/seu.%s.tree.markers.cts.pdf", out_fn, spi),
		tip_color_dictionary = tip_color_dictionary,
		assay = "RNA",
		test = "wilcox",
		pval_thr = 1e-3)
	mks_bipartition_cts$gene = gg_v [ mks_bipartition_cts$gene ]
	mks_bipartition_cts$gene_name = ogm_gv [mks_bipartition_cts$gene]
	mks_bipartition_cts$is_tf = mks_bipartition_cts$gene %in% gene_annot_tfs$gene
	write.table(mks_bipartition_cts, sprintf("%s/seu.%s.tree.markers.cts.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	# tree of leiden clusters
	SeuratObject::Idents(seu) = seu@meta.data[,"leiden"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca_integrated_harmony", dims = 1:10)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	tip_color_dictionary = dic_from_vecs(ctt$cluster, ctt$color)
	mks_bipartition_lei = sca_seurat_tree_markers(
		seu = seu, seu_phylo = seu_phylo,
		output_file = sprintf("%s/seu.%s.tree.markers.lei.pdf", out_fn, spi),
		tip_color_dictionary = tip_color_dictionary,
		assay = "RNA",
		test = "wilcox",
		pval_thr = 1e-3)
	mks_bipartition_lei$gene = gg_v [ mks_bipartition_lei$gene ]
	mks_bipartition_lei$gene_name = ogm_gv [mks_bipartition_lei$gene]
	mks_bipartition_lei$is_tf = mks_bipartition_lei$gene %in% gene_annot_tfs$gene
	write.table(mks_bipartition_lei, sprintf("%s/seu.%s.tree.markers.lei.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

	# tree of metacells
	SeuratObject::Idents(seu) = seu@meta.data[,"metacell"]
	seu = Seurat::BuildClusterTree(seu, assay = "RNA", reduction = "pca_integrated_harmony", dims = 1:10)
	seu_phylo = Seurat::Tool(seu, slot = "Seurat::BuildClusterTree")
	xtc = table(seu$metacell, seu$color)
	tip_color_dictionary = apply(xtc, 1, function(v) colnames(xtc)[which.max(v)])
	mks_bipartition_mcs = sca_seurat_tree_markers(
		seu = seu, seu_phylo = seu_phylo,
		output_file = sprintf("%s/seu.%s.tree.markers.mcs.pdf", out_fn, spi),
		tip_color_dictionary = tip_color_dictionary,
		assay = "RNA",
		test = "wilcox",
		pval_thr = 1e-2,
		width = 20,
		height = 120)
	mks_bipartition_mcs$gene = gg_v [ mks_bipartition_mcs$gene ]
	mks_bipartition_mcs$gene_name = ogm_gv [mks_bipartition_mcs$gene]
	mks_bipartition_mcs$is_tf = mks_bipartition_mcs$gene %in% gene_annot_tfs$gene
	write.table(mks_bipartition_mcs, sprintf("%s/seu.%s.tree.markers.mcs.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
	
	
	
	## Plot top markers per cluster (heatmaps) ##
	
	# metacells
	message(sprintf("metacell | %s | expression maps, heatmaps...", spi))
	top_markers = select_top_markers(mc_fp, matrix_thr = 1.5, n_top_markers = 20, n_markers_rollmean = 4)
	# top_markers = unique(as.character(apply(mc_fp, 2, function(v) { names(head(sort(v, decreasing = TRUE), 20)) } )))
	mc_fp_f = mc_fp[top_markers,ctm$metacell]
	order_rows = rev(names(sort(apply(mc_fp_f,1,function(x) which.max(x)))))
	mc_fp_f = mc_fp_f[order_rows,]
	top_cl_per_gene = colnames(mc_fp_f) [ apply(mc_fp_f, 1, which.max) ]
	rownames(mc_fp_f) = paste(as.character(gg_v [ rownames(mc_fp_f) ]) , ogm_gv [ as.character(gg_v [ rownames(mc_fp_f) ]) ], sep = " | ")
	rownames(mc_fp_f) = stringr::str_trunc(rownames(mc_fp_f), 60)
	pp1 = plot_complex_heatmap(
		mc_fp_f,
		name = "Footprint",
		cluster_row = FALSE,
		cluster_col = FALSE,
		color_min = 1,
		color_max = 8,
		color_mat = c("gray95","orange","orangered2","#520c52"),
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cl_per_gene,
		colors_row = ctm_mcs_col_v,
		categories_col = colnames(mc_fp_f),
		colors_col = ctm_mcs_col_v,
		fontsize = 5
	)
	pdf(sprintf("%s/seu.%s.markers.top.mcs.pdf", out_fn, spi), width = 8 + (ncol(mc_fp_f) / 30), height = 10 + (nrow(mc_fp_f) / 20))
	print(pp1)
	dev.off()
	
	# leiden
	top_markers = select_top_markers(le_fp, matrix_thr = 1.5, n_top_markers = 20, n_markers_rollmean = 2)
	# top_markers = unique(as.character(apply(le_fp, 2, function(v) { names(head(sort(v, decreasing = TRUE), 20)) } )))
	le_fp_f = le_fp[top_markers,ctt$cluster[ctt$cluster_class=="leiden"]]
	order_rows = rev(names(sort(apply(le_fp_f,1,function(x) which.max(x)))))
	le_fp_f = le_fp_f[order_rows,]
	top_cl_per_gene = colnames(le_fp_f) [ apply(le_fp_f, 1, which.max) ]
	rownames(le_fp_f) = paste(as.character(gg_v [ rownames(le_fp_f) ]) , ogm_gv [ as.character(gg_v [ rownames(le_fp_f) ]) ], sep = " | ")
	rownames(le_fp_f) = stringr::str_trunc(rownames(le_fp_f), 60)
	pp1 = plot_complex_heatmap(
		le_fp_f,
		name = "Footprint",
		cluster_row = FALSE,
		cluster_col = FALSE,
		color_min = 1,
		color_max = 8,
		color_mat = c("gray95","orange","orangered2","#520c52"),
		cell_border = gpar(col = NA, lwd = 1, lty = 1),
		heatmap_border = gpar(col = "black", lwd = 1, lty = 1),
		categories_row = top_cl_per_gene,
		colors_row = ctt_lou_col_v,
		categories_col = colnames(le_fp_f),
		colors_col = ctt_lou_col_v,
		fontsize = 5
	)
	pdf(sprintf("%s/seu.%s.markers.top.lei.pdf", out_fn, spi), width = 8 + (ncol(le_fp_f) / 30), height = 10 + (nrow(le_fp_f) / 20))
	print(pp1)
	dev.off()
	

	## Ion channel v GPCR signal ##
	
	message(sprintf("metacell | %s | ion v gpcr signal...", spi))
	
	# ion channels		
	gset_ion = read.table(sprintf("../results_annotation/results_gene_lists/%s_list_genes.ionchannel.txt", spi_w), sep = " ", header = FALSE, col.names = c("gene","family"))
	gset_ion_v = gsub("_","-", gset_ion$gene)
	gset_ion_v = intersect(gset_ion_v, rownames(mc_umifrac))
	
	# GPCR
	gset_gpc = read.table(sprintf("../results_annotation/results_gene_lists/%s_list_genes.ionchannel.txt", spi_w), sep = " ", header = FALSE, col.names = c("gene","family"))
	gset_gpc_v = gsub("_","-", gset_gpc$gene)
	gset_gpc_v = intersect(gset_gpc_v, rownames(mc_umifrac))
	
	# secreted (with signal peptide)		
	gset_sec = read.table(sprintf("../results_annotation/results_gene_lists/%s_list_genes.secreted.txt", spi_w), sep = " ", header = FALSE, col.names = c("gene"))
	gset_sec_v = gsub("_","-", gset_sec$gene)
	gset_sec_v = intersect(gset_sec_v, rownames(mc_umifrac))

	pdf(sprintf("%s/seu.%s.markers.neural_v_gland.pdf", out_fn, spi), height = 16, width = 12)
	layout(matrix(c(1,1,1,2,2,2,3,3,3,4,5,6), nrow = 4, byrow = TRUE))
	# umifrac aggregated plot, leiden
	umif_ion = colSums(le_umifrac[gset_ion_v,])
	umif_gpc = colSums(le_umifrac[gset_gpc_v,])
	umif_sec = colSums(le_umifrac[gset_sec_v,])
	# umifrac aggregated plot, leiden barplots
	barplot(umif_ion[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "Ion UMIfrac", border = NA, space = 0)
	barplot(umif_gpc[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "GPCR UMIfrac", border = NA, space = 0)
	barplot(umif_sec[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "Secreted UMIfrac", border = NA, space = 0)
	# umifrac aggregated plot, leiden
	plot(umif_ion[ctt$cluster[ctt$cluster_class=="leiden"]], umif_sec[ctt$cluster[ctt$cluster_class=="leiden"]], pch = 19, col = ctt$color[ctt$cluster_class=="leiden"], xlab = "Ion UMIfrac", ylab = "secreted UMIfrac", log = "xy", main = "Neural v gland signal, leiden")
	text(umif_ion[ctt$cluster[ctt$cluster_class=="leiden"]], umif_sec[ctt$cluster[ctt$cluster_class=="leiden"]], ctt$cluster[ctt$cluster_class=="leiden"], cex = 0.5, col = colorspace::darken(ctt$color, 0.6)[ctt$cluster_class=="leiden"])
	# umifrac aggregated plot, metacell
	umif_ion = colSums(mc_umifrac[gset_ion_v,])
	umif_gpc = colSums(mc_umifrac[gset_gpc_v,])
	umif_sec = colSums(mc_umifrac[gset_sec_v,])
	plot(umif_ion[ctm$metacell], umif_sec[ctm$metacell], pch = 19, col = ctm$color, xlab = "Ion UMIfrac", ylab = "secreted UMIfrac", log = "xy", main = "Neural v gland signal, metacell")
	text(umif_ion[ctm$metacell], umif_sec[ctm$metacell], ctm$metacell, cex = 0.5, col = colorspace::darken(ctm$color, 0.6))
	# umifrac aggregated plot, cell type
	umif_ion = colSums(ct_umifrac[gset_ion_v,])
	umif_hpc = colSums(ct_umifrac[gset_gpc_v,])
	umif_sec = colSums(ct_umifrac[gset_sec_v,])
	plot(umif_ion[levels(ctt$cell_type)], umif_sec[levels(ctt$cell_type)], pch = 19, col = ctt_cts_col_v[levels(ctt$cell_type)], xlab = "Ion UMIfrac", ylab = "secreted UMIfrac", log = "xy", main = "Neural v gland signal, cell type")
	text(umif_ion[levels(ctt$cell_type)], umif_sec[levels(ctt$cell_type)], levels(ctt$cell_type), cex = 0.5, col = colorspace::darken(ctt_cts_col_v[levels(ctt$cell_type)], 0.6))

	# num specific genes
	spef_ion = apply(le_fp[gset_ion_v,], 2, function(v){ length(which(v>1.5)) })
	spef_gpc = apply(le_fp[gset_gpc_v,], 2, function(v){ length(which(v>1.5)) })
	spef_sec = apply(le_fp[gset_sec_v,], 2, function(v){ length(which(v>1.5)) })
	# umifrac aggregated plot, leiden barplots
	barplot(spef_ion[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "Ion num specific", border = NA, space = 0)
	barplot(spef_gpc[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "GPCR num specific", border = NA, space = 0)
	barplot(spef_sec[ctt$cluster[ctt$cluster_class=="leiden"]], las = 2, col = ctt$color[ctt$cluster_class=="leiden"], main = "Secreted num specific", border = NA, space = 0)
	# umifrac aggregated plot, leiden
	plot(spef_ion[ctt$cluster[ctt$cluster_class=="leiden"]], spef_sec[ctt$cluster[ctt$cluster_class=="leiden"]], pch = 19, col = ctt$color[ctt$cluster_class=="leiden"], xlab = "Ion num specific", ylab = "secreted num specific", main = "Neural v gland signal, leiden")
	text(spef_ion[ctt$cluster[ctt$cluster_class=="leiden"]], spef_sec[ctt$cluster[ctt$cluster_class=="leiden"]], ctt$cluster[ctt$cluster_class=="leiden"], cex = 0.5, col = colorspace::darken(ctt$color, 0.6)[ctt$cluster_class=="leiden"])
	# spefrac aggregated plot, metacell
	spef_ion = apply(mc_fp[gset_ion_v,], 2, function(v){ length(which(v>1.5)) })
	spef_gpc = apply(mc_fp[gset_gpc_v,], 2, function(v){ length(which(v>1.5)) })
	spef_sec = apply(mc_fp[gset_sec_v,], 2, function(v){ length(which(v>1.5)) })
	plot(spef_ion[ctm$metacell], spef_sec[ctm$metacell], pch = 19, col = ctm$color, xlab = "Ion num specific", ylab = "secreted num specific", main = "Neural v gland signal, metacell")
	text(spef_ion[ctm$metacell], spef_sec[ctm$metacell], ctm$metacell, cex = 0.5, col = colorspace::darken(ctm$color, 0.6))
	# spefrac aggregated plot, cell type
	spef_ion = apply(ct_fp[gset_ion_v,], 2, function(v){ length(which(v>1.5)) })
	spef_gpc = apply(ct_fp[gset_gpc_v,], 2, function(v){ length(which(v>1.5)) })
	spef_sec = apply(ct_fp[gset_sec_v,], 2, function(v){ length(which(v>1.5)) })
	plot(spef_ion[levels(ctt$cell_type)], spef_sec[levels(ctt$cell_type)], pch = 19, col = ctt_cts_col_v[levels(ctt$cell_type)], xlab = "Ion num specific", ylab = "secreted num specific", main = "Neural v gland signal, cell type")
	text(spef_ion[levels(ctt$cell_type)], spef_sec[levels(ctt$cell_type)], levels(ctt$cell_type), cex = 0.5, col = colorspace::darken(ctt_cts_col_v[levels(ctt$cell_type)], 0.6))
	
	dev.off()
	

	## Differentiation status (Cytotrace) ##
	
	cytt = CytoTRACE::CytoTRACE(as.matrix(seu@assays$RNA$counts), ncores = 24, subsamplesize = 1000)
	seu@meta.data$cytotrace = cytt$CytoTRACE [ rownames(seu@meta.data) ]
	seu@meta.data$cytotrace_rank = cytt$CytoTRACErank [ rownames(seu@meta.data) ]
	seu@meta.data$cytotrace_num_genes = cytt$Counts [ rownames(seu@meta.data) ]
	seu@meta.data$cytotrace_gcs = cytt$GCS [ rownames(seu@meta.data) ]
	
	pdf(sprintf("%s/seu.%s.cytotrace.pdf", out_fn, spi), height = 8, width = 20)
	boxplot(cytotrace ~ cell_type, seu@meta.data, las = 2, col = ctt$color[!duplicated(ctt$cell_type)])
	boxplot(cytotrace ~ leiden, seu@meta.data, las = 2, col = ctt$color)
	boxplot(cytotrace ~ metacell, seu@meta.data, las = 2, col = ctm$color)
	cyt_mcs_v = aggregate(cytotrace ~ metacell, seu@meta.data, median)
	cyt_mcs_n = aggregate(nFeature_RNA ~ metacell, seu@meta.data, median)
	col_mcs_v = dic_from_vecs(ctm$metacell, ctm$color)
	plot(cyt_mcs_n[,2], cyt_mcs_v[,2], col = col_mcs_v[cyt_mcs_v[,1]], pch = 19, log = "x", cex = 0.5, ylab = "cytotrace score, median", xlab = "num transcripts, median", main = "Cytotrace per metacell", ylim = c(0,1))
	text(cyt_mcs_n[,2], cyt_mcs_v[,2], cyt_mcs_v[,1], cex = 0.5, col = scales::alpha(colorspace::darken(col_mcs_v[cyt_mcs_v[,1]]), 0.5))
	dev.off()

	# cytotrace
	saveRDS(cytt, sprintf("%s/dat.%s.cytotrace.rds", out_fn, spi))
	
	
	## Functional enrichments in markers ##
	
	go_annot = gsa_topgo_load_emapper(emapper_fn = sprintf("../data/reference/%s_ensembl.GO.csv", spi_w), index_col_GOs = 2)
	pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi_w))
	names(go_annot) = oga_gtv [ names(go_annot)  ]
	names(pf_annot) = oga_gtv [ names(pf_annot)  ]
	
	dir.create(sprintf("%s/tests/", out_fn), showWarnings = FALSE)
	for (loi in levels(seu@meta.data[,"cell_type"])) {
		
		# get markers
		mki = mks_cts [ mks_cts$focus_node == loi & mks_cts$p_val_adj < 0.01 & mks_cts$avg_log2FC > log2(1.5), ]
		mki_genes = mki$gene
		mki_genes = mki_genes [ !is.na(mki_genes) ]
		bgi_genes = gg_v [rownames(seu)]
		bgi_genes = bgi_genes [ !is.na(bgi_genes) ]
		
		message(sprintf("functional enrichments | %s, n=%i genes", loi, length(mki_genes)))
		if (length(mki_genes)>0) {
			# test pfam
			enrichment_table = gsa_enrichment_hypergeometric(
				annotation = pf_annot,
				genes_fg = mki_genes,
				genes_bg = bgi_genes,
				output_prefix = sprintf("%s/tests/enrichment.%s", out_fn, spi),
				name_fg = make.names(loi)
			)
			# test GOs
			go_enrichment_table = gsa_topgo_enrichment(
				annotation = go_annot,
				genes_fg = mki_genes,
				genes_bg = bgi_genes,
				output_prefix = sprintf("%s/tests/enrichment.%s", out_fn, spi),
				name_fg = make.names(loi),
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)

			# list GOs
			if (nrow(go_enrichment_table) > 0) {
				mki_genes_v = factor(as.numeric(mki_genes %in% mki_genes), levels = c(0,1))
				names(mki_genes_v) = mki_genes
				mki_genes_gos_d = data.frame()
				for (oti in c("BP","MF","CC")) {
					go_annot_object = suppressMessages(new("topGOdata", ontology = oti, allGenes = mki_genes_v, annot = topGO::annFUN.gene2GO, gene2GO = go_annot))
					mki_genes_gos = topGO::genesInTerm(go_annot_object)
					mki_genes_gos_d = rbind(
						mki_genes_gos_d,
						data.frame(
							GO = unlist(as.list(sapply(1:length(mki_genes_gos), function(vv) { rep(names(mki_genes_gos)[vv], lengths(mki_genes_gos)[vv])  } ))),
							gene = as.character(unlist(mki_genes_gos)),
							ontology = oti
						)
					)
				}
				mki_genes_gos_d$gene_name = ogm_gv [mki_genes_gos_d$gene]
				mki_genes_gos_d = mki_genes_gos_d [ order(mki_genes_gos_d$ontology, mki_genes_gos_d$GO), ]
				write.table(mki_genes_gos_d, sprintf("%s/tests/enrichment.%s.%s.topgo.genes.csv", out_fn, spi, make.names(loi)), sep = "\t", quote = FALSE, row.names = FALSE)
			}
			
		}
	}
	
	
	## Save outputs ##
	
	# end
	message(sprintf("metacell | %s | save object...", spi))
	saveRDS(seu, sprintf("%s/dat.%s.seurat_final.rds", out_fn, spi))
	write.table(seu@meta.data, sprintf("%s/dat.%s.cell_metadata.csv", out_fn, spi), sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
		
}

message("All done!")
