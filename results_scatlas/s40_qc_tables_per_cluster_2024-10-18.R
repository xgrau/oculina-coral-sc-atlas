# libraries
suppressMessages(library("Seurat"))
suppressMessages(library("SeuratObject"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
graphics.off()


# data index
sps_list = c("Ocupat","Spin","Spis","Amil")

# load orthology
ogm = read.table("../data/orthology_Metazoa/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)
oga_v = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
gg_v = dic_from_vecs(gsub("_","-",oga$gene), oga$gene)

# load gene names
gna_fn = "../data/orthology_Metazoa/orthogroup_conservation.csv"
gna = read.table(gna_fn, sep = "\t", header = TRUE)

# dictionaries transcript to gene
gna_gtv = dic_from_vecs(names = gna$transcript, terms = gna$gene)
gna_tgv = dic_from_vecs(names = gna$gene, terms = gna$transcript)



# output
out_fn = sprintf("results_qc_tables/")
dir.create(out_fn, showWarnings = FALSE)


# loop over species
for (spi in sps_list) { 
	
	# set working species
	if (spi == "Spio" | spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}

	## Load annotations ##

	# input
	inp_fn = sprintf("results_metacell_%s_filt/", spi)

	# load gene annotations for this species
	message(sprintf("qc tables | %s load", spi))
	gna_i = gna [ gna$species == spi, ]
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
	list_tfs = gsub("_","-", names(gene_annot_tfs_v))
	gene_annot_v = dic_from_vecs(gsub("_","-",gene_annot$gene), gene_annot$name)
	
	# cell types
	message(sprintf("qc tables | %s ct info", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", inp_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$cluster = factor(ctt$cluster, levels = unique(ctt$cluster))
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
	ctt_cts_bct_v = dic_from_vecs(ctt$cell_type, ctt$broad_cell_type)
	
	# metacell
	message(sprintf("qc tables | %s mc info", spi))
	ctm = read.table(sprintf("%s/annot.%s.mcs.csv", inp_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctm$metacell = factor(ctm$metacell, levels = unique(ctm$metacell))
	ctm$broad_cell_type = ctt_cts_bct_v [ ctm$metacell_top_cell_type ]
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
	ctm_mcs_cts_v = dic_from_vecs(ctm$metacell, ctm$metacell_top_cell_type)
	
	# seurat object with batch info
	message(sprintf("qc tables | %s load Seurat...", spi))
	seu = readRDS(sprintf("%s/dat.%s.seurat_final.rds", inp_fn, spi))
	seu$batch_method = factor(seu$batch_method)

	# expression
	ct_fp = readRDS(sprintf("%s/dat.%s.expression.cts_fp.rds", inp_fn, spi))
	mc_fp = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", inp_fn, spi))
	
	# num tfs specific to each cluster
	ct_fp_f = ct_fp [ rownames(ct_fp) %in% list_tfs, ]
	mc_fp_f = mc_fp [ rownames(mc_fp) %in% list_tfs, ]
	ct_fp_num_specific_tfs = apply(ct_fp_f, 2, function(vv) { sum(vv >= 1.5) })
	mc_fp_num_specific_tfs = apply(mc_fp_f, 2, function(vv) { sum(vv >= 1.5) })
	
	# doublet scores	
	dbs_scs = read.table(sprintf("results_metacell_%s_prefilt/seu.%s.doublet_report.cell.csv", spi, spi), header = TRUE)
	dbs_scs$metacell = factor(dbs_scs$metacell, levels = levels(seu$metacell))
	dbs_scs$cell_type = factor(ctm_mcs_cts_v[dbs_scs$metacell], levels = levels(seu$cell_type))
	dbs_scs$class = factor(dbs_scs$class, levels = c("singlet","doublet"))

	# doublet info from clicktags
	dbs_ct_fl = list.files("mapping/scdb_clicktag", pattern = "*.tsv", full.names = TRUE)
	dbs_ct_fl = dbs_ct_fl [ grepl(spi_w, dbs_ct_fl) ]
	dbs_ct_d = data.frame()
	for (ffi in dbs_ct_fl) {
		dbs_ct_d = rbind(dbs_ct_d, read.table(ffi, header = TRUE))
	}
	if (nrow(dbs_ct_d) == 0) {
		dbs_ct_d = data.frame(matrix(nrow = nrow(seu@meta.data), ncol = 2), row.names = rownames(seu@meta.data))
		colnames(dbs_ct_d) = c("ct_relative_size_fn","ct_total_counts")
		dbs_ct_d[,1] = dbs_ct_d[,2] = 0
	}
	seu_ct_v = dic_from_vecs(rownames(seu@meta.data), seu@meta.data$cell_type)
	seu_mc_v = dic_from_vecs(rownames(seu@meta.data), seu@meta.data$metacell)
	dbs_ct_d$metacell = seu_mc_v [ rownames(dbs_ct_d) ]
	dbs_ct_d$cell_type = seu_ct_v [ rownames(dbs_ct_d) ]


	# qc summary table per metacell
	qcmc_d = data.frame(
		metacell        = ctm$metacell,
		cell_type       = ctm$metacell_top_cell_type,
		broad_cell_type = ctm$broad_cell_type,
		num_cells        = as.numeric(table(seu$metacell)[ctm$metacell]),
		UMI_per_cell_mean = aggregate(nCount_RNA ~ metacell, data = seu@meta.data, mean)[ctm$metacell,2],
		UMI_per_cell_Q1   = aggregate(nCount_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.25))[ctm$metacell,2],
		UMI_per_cell_Q2   = aggregate(nCount_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.50))[ctm$metacell,2],
		UMI_per_cell_Q3   = aggregate(nCount_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.75))[ctm$metacell,2],
		gene_per_cell_mean = aggregate(nFeature_RNA ~ metacell, data = seu@meta.data, mean)[ctm$metacell,2],
		gene_per_cell_Q1   = aggregate(nFeature_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.25))[ctm$metacell,2],
		gene_per_cell_Q2   = aggregate(nFeature_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.50))[ctm$metacell,2],
		gene_per_cell_Q3   = aggregate(nFeature_RNA ~ metacell, data = seu@meta.data, function(vv) quantile(vv, 0.75))[ctm$metacell,2],
		num_specific_TFs   = mc_fp_num_specific_tfs [ ctm$metacell ],
		doublet_score_median = aggregate(score ~ metacell, data = dbs_scs, drop = FALSE, function(vv) quantile(vv, 0.5))[ctm$metacell,2],
		clicktag_foldchange_median   = aggregate(ct_relative_size_fn ~ metacell, data = dbs_ct_d, function(vv) quantile(vv, 0.5))[ctm$metacell,2],
		clicktag_UMI_per_cell_median = aggregate(ct_total_counts ~ metacell, data = dbs_ct_d, function(vv) quantile(vv, 0.5))[ctm$metacell,2],
		fraction_ribo_median = aggregate(percent_ribo ~ metacell, data = seu@meta.data, function(vv) quantile(vv/100, 0.5))[ctm$metacell,2],
		fraction_cells_per_batch = NA
	)
	
	# string of batch composition
	tt = table(seu@meta.data$metacell, seu@meta.data$batch)
	tt = tt / rowSums(tt)[ctm$metacell]
	batch_string_per_cluster = sapply(1:nrow(tt), function(vv) { paste(sprintf("%s:%.2f", colnames(tt), tt[vv,]), collapse = ";") })
	qcmc_d$fraction_cells_per_batch = batch_string_per_cluster
	
	# save
	write.table(qcmc_d, sprintf("%s/seu.%s.qc_table.mcs.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)
	
	# qc summary table per cell_type
	qcct_d = data.frame(
		cell_type       = levels(ctt$cell_type),
		broad_cell_type = ctt_cts_bct_v[levels(ctt$cell_type)],
		num_cells        = as.numeric(table(seu$cell_type)),
		UMI_per_cell_mean = aggregate(nCount_RNA ~ cell_type, data = seu@meta.data, mean)[,2],
		UMI_per_cell_Q1   = aggregate(nCount_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.25))[,2],
		UMI_per_cell_Q2   = aggregate(nCount_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.50))[,2],
		UMI_per_cell_Q3   = aggregate(nCount_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.75))[,2],
		gene_per_cell_mean = aggregate(nFeature_RNA ~ cell_type, data = seu@meta.data, mean)[,2],
		gene_per_cell_Q1   = aggregate(nFeature_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.25))[,2],
		gene_per_cell_Q2   = aggregate(nFeature_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.50))[,2],
		gene_per_cell_Q3   = aggregate(nFeature_RNA ~ cell_type, data = seu@meta.data, function(vv) quantile(vv, 0.75))[,2],
		num_specific_TFs   = ct_fp_num_specific_tfs,
		doublet_score_median = aggregate(score ~ cell_type, data = dbs_scs, drop = FALSE, function(vv) quantile(vv, 0.5))[,2],
		clicktag_foldchange_median   = aggregate(ct_relative_size_fn ~ cell_type, data = dbs_ct_d, function(vv) quantile(vv, 0.5))[,2],
		clicktag_UMI_per_cell_median = aggregate(ct_total_counts ~ cell_type, data = dbs_ct_d, function(vv) quantile(vv, 0.5))[,2],
		fraction_ribo_median = aggregate(percent_ribo ~ cell_type, data = seu@meta.data, function(vv) quantile(vv/100, 0.5))[,2],
		fraction_cells_per_batch = NA
	)
	
	# string of batch composition
	tt = table(seu@meta.data$cell_type, seu@meta.data$batch)
	tt = tt / rowSums(tt)
	batch_string_per_cluster = sapply(1:nrow(tt), function(vv) { paste(sprintf("%s:%.2f", colnames(tt), tt[vv,]), collapse = ";") })
	qcct_d$fraction_cells_per_batch = batch_string_per_cluster
	
	# save
	write.table(qcct_d, sprintf("%s/seu.%s.qc_table.cts.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)

}