# libraries
suppressMessages(library("Seurat"))
suppressMessages(library("SeuratObject"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
graphics.off()

# data index
sps_list = c("Ocupat","Spin","Amil")

# load orthology
ogm = read.table("../data/orthology_Metazoa/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)

# load gene names
gna_fn = "../data/orthology_Metazoa/orthogroup_conservation.csv"
gna = read.table(gna_fn, sep = "\t", header = TRUE)

# dictionaries transcript to gene
gna_gtv = dic_from_vecs(names = gna$transcript, terms = gna$gene)
gna_tgv = dic_from_vecs(names = gna$gene, terms = gna$transcript)

# assignment threshold
ass_thr = 0.9

# loop over species
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


	## Load seurat ##

	# seurat object with batch info
	message(sprintf("check host | %s | load Seurat...", spi))
	seu = readRDS(sprintf("%s/dat.%s.seurat_final.rds", out_fn, spi))
	seu$batch_method = factor(seu$batch_method)

	## Load annotations ##

	# load gene annotations for this species
	message(sprintf("check host | load %s", spi))
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
	
	# seurat gene dict
	gg_v = dic_from_vecs(gsub("_","-", gna_p$gene), gna_p$gene)
	
	# cell types
	message(sprintf("check host | load %s ct info", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt$cluster = factor(ctt$cluster, levels = unique(ctt$cluster))
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
	ctt_cts_bct_v = dic_from_vecs(ctt$cell_type, ctt$broad_cell_type)
	ctt_leiden    = ctt$cluster [ ctt$cluster_class == "leiden" ]

	# metacell
	message(sprintf("check host | load %s mc info", spi))
	ctm = read.table(sprintf("%s/annot.%s.mcs.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctm$metacell = factor(ctm$metacell, levels = unique(ctm$metacell))
	ctm$broad_cell_type = ctt_cts_bct_v [ ctm$metacell_top_cell_type ]
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)


	## Load host cells ##
	
	message(sprintf("check host | load %s mars", spi))

	list_mars_mats = list.files(sprintf("mapping/map_%s_mars/umi.tab/", spi_w), full.names = TRUE)
	for (ffi in list_mars_mats) {  
		message(sprintf("check host | load %s mars, %s", spi, ffi))
		mdi = as.data.frame(data.table::fread(ffi))
		mdi_genes = as.character(mdi$V1)
		mdi$V1 = NULL
		mmi = as.matrix(mdi)
		rownames(mmi) = mdi_genes
		if (ffi == list_mars_mats[1]) {
			mmt = mmi
		} else {
			mmt = cbind(mmt, mmi)
		}
	}
	mah = mmt [ grepl(sprintf("^%s_", spi_w), rownames(mmt)), ]
	maa = mmt [ !grepl(sprintf("^%s_", spi_w), rownames(mmt)), ]
	rownames(maa) [ grepl("^Smic[^_]", rownames(maa)) ] = paste("Symmic", rownames(maa) [ grepl("^Smic[^_]", rownames(maa)) ], sep = "_")
	maa = maa [ !grepl("^orphan", rownames(maa)) , ]

	# create seurat object with batch info
	seu_h = Seurat::CreateSeuratObject(mah, project = "mars_host")
	seu_a = Seurat::CreateSeuratObject(maa, project = "mars_alga")
	seu_h@meta.data$batch_method = "mars"

	# thresholds
	if (spi_w == "Ocupat") {
		umi_host_thr = 250
		umi_host_thr_max = 5e3
	} else if (spi_w == "Amil") { 
		umi_host_thr = 200
		umi_host_thr_max = 5e3
	} else if (spi_w == "Spis") { 
		umi_host_thr = 100
		umi_host_thr_max = 5e3
	}

	# plot
	message(sprintf("check host | plot qcs %s", spi))
	pdf(sprintf("%s/seu.%s.hostcells.mars_qc.pdf", out_fn, spi), height = 6, width = 6)
	hist(log10(seu_h@meta.data$nCount_RNA), breaks = 80, main = "cell sizes to host")
	abline(v=log10(umi_host_thr), lty = 2, col = "red")
	hist(log10(seu_a@meta.data$nCount_RNA), breaks = 80, main = "cell sizes to alga")
	plot(
		seu_h@meta.data$nCount_RNA,
		seu_a@meta.data$nCount_RNA,
		pch = 19,
		col = ifelse(seu_h@meta.data$nCount_RNA >= umi_host_thr, "blue", "gray"),
		cex = 0.6,
		log ="xy",
		xlab = "cell sizes to host", ylab = "cell sizes to alga",
		sub = sprintf("%i out of %i cells remain at umi >= %i", sum(seu_h@meta.data$nCount_RNA >= umi_host_thr), dim(seu_h)[2], umi_host_thr)
	)
	abline(v=umi_host_thr, lty = 2, col = "red")
	dev.off()
	
	# drop bad cells
	message(sprintf("check host | filter bad cells %s", spi))
	seu_h = subset(seu_h, subset = nCount_RNA >= umi_host_thr & nCount_RNA < umi_host_thr_max)
	
	# add algal signal per cell
	seu_h@meta.data$algal_umi = seu_a$nCount_RNA [ rownames(seu_h@meta.data) ]

	# map query to reference
	seu_r = seu
	Seurat::DefaultAssay(seu_r) = "RNA"
	Seurat::DefaultAssay(seu_h) = "RNA"

	# subset to common features
	shared_feat = intersect(rownames(seu_r), rownames(seu_h))
	Seurat::VariableFeatures(seu_r) = seu_r@assays$SCT@var.features [ seu_r@assays$SCT@var.features %in% shared_feat ]
	seu_r$SCT = NULL
	seu_r = seu_r [ shared_feat , ]
	seu_h = seu_h [ shared_feat , ]

	# normalize
	message(sprintf("check host | seurat normalise %s", spi))
	seu_h = Seurat::NormalizeData(seu_h)
	seu_r = Seurat::NormalizeData(seu_r)

	# find anchors
	message(sprintf("check host | seurat transfer anchors %s", spi))
	seu_num_pcs_v = find_pca_elbow(seu_r@reductions$pca@stdev)
	seu_num_pcs = seu_num_pcs_v["First derivative"]
	seu_anchors = Seurat::FindTransferAnchors(
		reference = seu_r,
		query = seu_h,
		dims = 1:seu_num_pcs,
		reduction = "pcaproject",
		normalization.method = "LogNormalize",
		reference.reduction = "pca")

	# transfer data
	message(sprintf("check host | seurat transfer data %s", spi))
	seu_h_transfer_le = Seurat::TransferData(anchorset = seu_anchors, refdata = seu_r$leiden, dims = 1:seu_num_pcs, n.trees = 500)
	seu_h_transfer_mc = Seurat::TransferData(anchorset = seu_anchors, refdata = seu_r$metacell, dims = 1:seu_num_pcs, n.trees = 500)
	seu_h_transfer_ct = Seurat::TransferData(anchorset = seu_anchors, refdata = seu_r$cell_type, dims = 1:seu_num_pcs, n.trees = 500)
	seu_h@meta.data$cell_type = factor(seu_h_transfer_ct$predicted.id, levels = levels(seu$cell_type))
	seu_h@meta.data$leiden = factor(seu_h_transfer_le$predicted.id, levels = levels(seu$leiden))
	seu_h@meta.data$metacell = factor(seu_h_transfer_mc$predicted.id, levels = levels(seu$metacell))
	seu_h@meta.data$cell_type_score = seu_h_transfer_ct$prediction.score.max
	seu_h@meta.data$leiden_score = seu_h_transfer_le$prediction.score.max         
	seu_h@meta.data$metacell_score = seu_h_transfer_mc$prediction.score.max

	# filter
	seu_h_d = seu_h@meta.data
	seu_h_d_f = seu_h_d [ seu_h_d$cell_type_score >= ass_thr,  ]

	# save transfer scores
	message(sprintf("check host | write transfer scores %s", spi))
	seu_h_transfer_mc_f = seu_h_transfer_mc
	seu_h_transfer_le_f = seu_h_transfer_le
	seu_h_transfer_ct_f = seu_h_transfer_ct
	colnames(seu_h_transfer_mc_f) = gsub("^prediction.score.","", colnames(seu_h_transfer_mc_f))
	colnames(seu_h_transfer_le_f) = gsub("^prediction.score.","", colnames(seu_h_transfer_le_f))
	colnames(seu_h_transfer_ct_f) = gsub("^prediction.score.","", colnames(seu_h_transfer_ct_f))
	write.table(seu_h_transfer_mc_f, sprintf("%s/seu.%s.hostcells.scores.mcs.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(seu_h_transfer_le_f, sprintf("%s/seu.%s.hostcells.scores.lei.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)
	write.table(seu_h_transfer_ct_f, sprintf("%s/seu.%s.hostcells.scores.cts.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)

	# plot
	message(sprintf("check host | plot signal per cluster %s", spi))

	# summarise algal signal per metacell
	pdf(sprintf("%s/seu.%s.hostcells.scores.mcs.pdf", out_fn, spi), height = 4, width = 24)
	barplot(colSums(seu_h_transfer_mc_f[,levels(ctm$metacell)]), col = ctm_mcs_col_v[levels(ctm$metacell)], space = 0, las = 2, ylab = "sum transfer score", cex.names = 0.3)
	title(main = sprintf("transfer to %s metacells", spi))
	barplot(apply(seu_h_transfer_mc_f[,levels(ctm$metacell)], 2, mean), col = ctm_mcs_col_v[levels(ctm$metacell)], space = 0, las = 2, ylab = "mean transfer score", cex.names = 0.4)
	barplot(apply(seu_h_transfer_mc_f[,levels(ctm$metacell)], 2, function(v) { sum(v >= 0.5) }), col = ctm_mcs_col_v[levels(ctm$metacell)], space = 0, las = 2, ylab = "num positive transfer score", cex.names = 0.4)
	dev.off()

	# summarise algal signal per cell type
	pdf(sprintf("%s/seu.%s.hostcells.scores.cts.pdf", out_fn, spi), height = 4, width = 12)
	barplot(colSums(seu_h_transfer_ct_f[,make.names(levels(ctt$cell_type))]), col = ctt_cts_col_v[levels(ctt$cell_type)], space = 0, las = 2, ylab = "sum transfer score", cex.names = 0.5)
	title(main = sprintf("transfer to %s cell types", spi))
	barplot(apply(seu_h_transfer_ct_f[,make.names(levels(ctt$cell_type))], 2, mean), col = ctt_cts_col_v[levels(ctt$cell_type)], space = 0, las = 2, ylab = "mean transfer score", cex.names = 0.5)
	barplot(apply(seu_h_transfer_ct_f[,make.names(levels(ctt$cell_type))], 2, function(v) { sum(v >= ass_thr) }), col = ctt_cts_col_v[levels(ctt$cell_type)], space = 0, las = 2, ylab = "num positive transfer score", cex.names = 0.5)
	dev.off()

	# summarise algal signal per leiden
	pdf(sprintf("%s/seu.%s.hostcells.scores.lei.pdf", out_fn, spi), height = 4, width = 12)
	barplot(colSums(seu_h_transfer_le_f[,as.character(ctt_leiden)]), col = ctt_lou_col_v[as.character(ctt_leiden)], space = 0, las = 2, ylab = "sum transfer score", cex.names = 0.5)
	title(main = sprintf("transfer to %s leiden", spi))
	barplot(apply(seu_h_transfer_le_f[,as.character(ctt_leiden)], 2, mean), col = ctt_lou_col_v[as.character(ctt_leiden)], space = 0, las = 2, ylab = "mean transfer score", cex.names = 0.5)
	barplot(apply(seu_h_transfer_le_f[,as.character(ctt_leiden)], 2, function(v) { sum(v >= ass_thr	) }), col = ctt_lou_col_v[as.character(ctt_leiden)], space = 0, las = 2, ylab = "num positive transfer score", cex.names = 0.5)
	dev.off()
	
	# algal signal per cell: gastrodermis v alga
	gastro_types = levels(ctt$cell_type) [ grepl("gastrodermis",levels(ctt$cell_type)) ]
	gasalg_types = gastro_types [ grepl("alga_hosting", gastro_types) ]
	gastro_types = gastro_types [ !gastro_types %in% gasalg_types  ]

	pdf(sprintf("%s/seu.%s.hostcells.pairwise_scores.cts.pdf", out_fn, spi), height = 6, width = 6)
	for (gast_i in gastro_types) {
		for (gaal_i in gasalg_types) {
			plot(seu_h_transfer_ct_f[,gaal_i], seu_h_transfer_ct_f[,gast_i], col = scales::alpha("blue", 0.5), pch = 19, cex = 1, xlim = c(0,1), ylim = c(0,1), xlab = gaal_i, ylab = gast_i)
			title(main = sprintf("%s: cell transfer score %s v %s", spi, gast_i, gaal_i))
			abline(v = ass_thr, h = ass_thr, lty = 2, col = "red")
		}
	}
	dev.off()

	# summarise algal signal per cell type
	message(sprintf("check host | write signal per cluster %s", spi))
	asig = data.frame(
		cell_type = levels(seu_h_d_f$cell_type),
		num_cells = as.numeric(table(seu_h_d_f$cell_type)),
		num_umis_alga = aggregate(algal_umi ~ cell_type, seu_h_d_f, drop = FALSE, sum)[,2],
		num_umis_host = aggregate(nCount_RNA ~ cell_type, seu_h_d_f, drop = FALSE, sum)[,2]
	)
	asig$num_umis_alga [ is.na(asig$num_umis_alga) ] = 0
	asig$num_umis_host [ is.na(asig$num_umis_host) ] = 0
	# asig$umis_alga_per_cell = asig$num_umis_alga / asig$num_cells
	asig$fraction_cells_alga = as.numeric(sprintf("%.5f",asig$num_cells / sum(asig$num_cells)))
	asig$fraction_umis_alga = as.numeric(sprintf("%.5f",asig$num_umis_alga / sum(asig$num_umis_alga)))
	asig$fraction_umis_host = as.numeric(sprintf("%.5f",asig$num_umis_host / sum(asig$num_umis_host)))
	# asig$ratio_umi_alga_host = asig$num_umis_alga / asig$num_umis_host
	
	# save
	write.table(asig, sprintf("%s/seu.%s.hostcells.assignment.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)

	# plots of ct transfer
	message(sprintf("check host | transfer plots %s", spi))
	pdf(sprintf("%s/seu.%s.hostcells.assignment.pdf", out_fn, spi), height = 6, width = 6)
	
	# pie
	pie_label = sprintf("%s\nn=%i(%.1fpc)", asig$cell_type, asig$num_cells, asig$fraction_cells_alga*100)
	pie_label [ grepl("n=0", pie_label) ] = ""
	pie(asig$num_cells, col = ctt_cts_col_v [asig$cell_type], labels = pie_label)
	title(main = spi)
	title(sub = sprintf("n=%i alga-hosting cells, mapped", sum(asig$num_cells)))
	
	# cell size
	message(sprintf("check host | cell size plots %s", spi))
	plot(
		seu_h_d_f$nCount_RNA,
		seu_a@meta.data[rownames(seu_h_d_f),]$nCount_RNA,
		pch = 19,
		col = ctt_cts_col_v [ seu_h_d_f$cell_type ],
		cex = 0.6,
		log ="xy",
		xlab = "cell sizes to host", ylab = "cell sizes to alga", main = "cell size"
	)
	ppl_label = sprintf("%s | n=%i", asig$cell_type, asig$num_cells)
	legend("bottomright", col = ctt_cts_col_v[asig$cell_type[asig$num_cells>0]], legend = ppl_label[asig$num_cells>0], pch = 19, bty = "n", cex = 0.5)
	dev.off()


	# summarise signal per alga	
	message(sprintf("check host | per-alga signal plots %s", spi))
	if (spi != "Spis") {
		
		alg_sps = unique(gsub("-.*", "", rownames(seu_a@assays$RNA$counts)))
		alg_sps = alg_sps [ alg_sps != "ERCC" ]
		alg_sps = alg_sps [ !grepl("^Spis", alg_sps) ]
		alg_sps = alg_sps [ !alg_sps %in% c("Trnag","Trnap") ]
		alg_col = viridisLite::viridis(length(alg_sps), begin = 0.1, end = 0.8)
		
		for (syi in alg_sps) {
			mat_syi = seu_a@assays$RNA$counts [ grepl(sprintf("^%s-", syi), rownames(seu_a@assays$RNA$counts)), ]
			mat_syi_cz = colSums(mat_syi)
			seu_h_d_f[,sprintf("algal_umi_%s", syi)] = mat_syi_cz [ rownames(seu_h_d_f) ]
		}
		
		# open
		pdf(sprintf("%s/seu.%s.hostcells.alga_per_ct.pdf", out_fn, spi), height = 12, width = 12)
		layout(1:2)

		# aggregate UMIs per species
		signal_per_cell = t(as.matrix(seu_h_d_f[ , sprintf("algal_umi_%s", alg_sps)]))
		signal_per_cell_per_ct = apply(signal_per_cell, 1, function(vv) tapply(vv, seu_h_d_f$cell_type, function(vvv) sum(vvv, na.rm = TRUE)))
		signal_per_cell_per_ct [ is.na(signal_per_cell_per_ct) ] = 0
		barplot(t(signal_per_cell_per_ct), las = 2, col  = alg_col)
		title(main = sprintf("num alga UMIs per %s cts", spi))
		legend("topright", legend = alg_sps, fill = alg_col)
		plot(0, bty="n", frame.plot=FALSE, xaxt = "n", yaxt = "n")
		
		# cells with a min num of UMIs per species
		signal2_per_cell_per_ct = apply(signal_per_cell, 1, function(vv) tapply(vv, seu_h_d_f$cell_type, function(vvv) sum(sum(vvv > 10), na.rm = TRUE)))
		signal2_per_cell_per_ct [ is.na(signal2_per_cell_per_ct) ] = 0
		barplot(t(signal2_per_cell_per_ct), las = 2, col  = alg_col)
		title(main = sprintf("num cells with >10 UMI per %s cts and alga sps", spi))
		legend("topright", legend = alg_sps, fill = alg_col)
		plot(0, bty="n", frame.plot=FALSE, xaxt = "n", yaxt = "n")
		
		# loop per algal species
		for (syi in alg_sps) {
			boxplot(seu_h_d_f[,sprintf("algal_umi_%s", syi)] ~ seu_h_d_f$cell_type, las = 2, cex.names = 0.5, col = ctt_cts_col_v[levels(seu_h_d_f$cell_type)], ylab = "UMI")
			title(main = sprintf("%s signal in %s cts", syi, spi))
			seu_h_d_a = aggregate(formula(sprintf("algal_umi_%s ~ cell_type", syi)), data = seu_h_d_f, function(vv) { sum(vv > 10) }, drop = FALSE)
			seu_h_d_a[,2] [ is.na(seu_h_d_a[,2]) ] = 0
			barplot(seu_h_d_a[,2], las = 2, cex.names = 0.5, col = ctt_cts_col_v[levels(seu_h_d_f$cell_type)], ylab = sprintf("num cells >10 %s UMI", syi))
		}
		
		# top alga per cell
		signal_per_cell_f = signal_per_cell [ , colSums(signal_per_cell) >= 50 ]
		xtt = table(rownames(signal_per_cell_f)[apply(signal_per_cell_f, 2, which.max)])
		xtl = sprintf("%s\nn=%i(%.1fpc)", names(xtt), xtt, 100*xtt/sum(xtt))
		pie(xtt, xtl, main = "top alga per cell (only those with >50UMIs)")
		
		dev.off()
		
		# UMI signal alga v alga
		pdf(sprintf("%s/seu.%s.hostcells.alga_v_alga.pdf", out_fn, spi), height = 6, width = 6)
		for (syi in alg_sps) {
			for	 (syj in alg_sps) {
				maxy = log10(max(c(signal_per_cell[sprintf("algal_umi_%s", syi),], signal_per_cell[sprintf("algal_umi_%s", syj),])))
				plot(log10(signal_per_cell[sprintf("algal_umi_%s", syi),]), log10(signal_per_cell[sprintf("algal_umi_%s", syj),]), las = 2, cex = 0.6, col = ctt_cts_col_v[as.character(seu_h_d_f$cell_type)], xlab = sprintf("UMI %s", syi), ylab = sprintf("UMI %s", syj), pch = 19, xlim = c(0, maxy), ylim = c(0, maxy))
				abline(a = 0, b = 1, col = "black", lty = 2)
			}
		}
		dev.off()

		
	}
	
	# write
	write.table(seu_h_d, sprintf("%s/seu.%s.hostcells.cell_metadata.csv", out_fn, spi), sep = "\t", row.names = FALSE, quote = FALSE)
	saveRDS(seu_h,       sprintf("%s/seu.%s.hostcells.seurat_object.rds", out_fn, spi))
	saveRDS(seu_anchors, sprintf("%s/seu.%s.hostcells.seurat_anchors.rds", out_fn, spi))
	
	
	
	## Dimensionality reduction ##

	message(sprintf("metacell | %s | recolor umap...", spi))

	# first, calculate medioids
	# cell type centroids
	ct_centroids_d = data.frame(
		umap1 = seu@reductions$umap_integrated@cell.embeddings[,1], 
		umap2 = seu@reductions$umap_integrated@cell.embeddings[,2], 
		label = seu@meta.data$cell_type)
	ct_centroids_a = aggregate(cbind(ct_centroids_d$umap1,ct_centroids_d$umap2) ~ label, ct_centroids_d, median)
	# metacell centroids
	mc_centroids_d = data.frame(
		umap1 = seu@reductions$umap_integrated@cell.embeddings[,1], 
		umap2 = seu@reductions$umap_integrated@cell.embeddings[,2], 
		label = seu@meta.data$metacell)
	mc_centroids_a = aggregate(cbind(mc_centroids_d$umap1,mc_centroids_d$umap2) ~ label, mc_centroids_d, median)
	mc_centroids_a_x_v = dic_from_vecs(mc_centroids_a[,1], mc_centroids_a[,2])
	mc_centroids_a_y_v = dic_from_vecs(mc_centroids_a[,1], mc_centroids_a[,3])

	# positions of query metacells (good only)
	set.seed(1)
	x_jitter_amount = (max(range(seu@reductions$umap_integrated@cell.embeddings[,1])) - min(range(seu@reductions$umap_integrated@cell.embeddings[,1]))) / 60
	y_jitter_amount = (max(range(seu@reductions$umap_integrated@cell.embeddings[,2])) - min(range(seu@reductions$umap_integrated@cell.embeddings[,2]))) / 60
	query_mc_d = data.frame(
		assigned_metacell = seu_h_d_f$metacell,
		assigned_x = jitter(mc_centroids_a_x_v[ seu_h_d_f$metacell ], amount = x_jitter_amount),
		assigned_y = jitter(mc_centroids_a_y_v[ seu_h_d_f$metacell ], amount = y_jitter_amount))

	# open umap
	pdf(sprintf("%s/seu.%s.hostcells.umap.pdf", out_fn, spi), height = 12, width = 12)

	# plot reference cells
	layout(matrix(1:4, nrow = 2, byrow = TRUE))
	plot(
		seu@reductions$umap_integrated@cell.embeddings[,1], seu@reductions$umap_integrated@cell.embeddings[,2], 
		col = "seashell3", pch = 19, cex = 0.5, main = "integrated, by cell type", xlab = "UMAP1", ylab = "UMAP2")
	
	# plot reference metacells
	points(mc_centroids_a[,2], mc_centroids_a[,3], font = 2, col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a[,1])], 0.4), bg = ctm_mcs_col_v[as.character(mc_centroids_a[,1])], pch = 21)
	
	# plot query cells around assigned metacells
	points(query_mc_d[,2], query_mc_d[,3], pch = 19, col = scales::alpha("red", 0.4), cex = 0.25)
	
	# another plot with cell-wise colors
	layout(1)
	plot(
		seu@reductions$umap_integrated@cell.embeddings[,1], seu@reductions$umap_integrated@cell.embeddings[,2],
		col = seu@meta.data$color, pch = 19, cex = 0.5, main = "integrated, by cell type", xlab = "UMAP1", ylab = "UMAP2")
	# text(mc_centroids_a[,2], mc_centroids_a[,3], mc_centroids_a[,1], font = 1, col = colorspace::darken(ctm_mcs_col_v[as.character(mc_centroids_a[,1])], 0.4), cex = 0.4)
	text(ct_centroids_a[,2], ct_centroids_a[,3], ct_centroids_a[,1], font = 2, col = scales::alpha("gray5", 0.7), cex = 0.7)
	title(sub = sprintf("n=%i cells into %i metacells", ncol(seu), nrow(ctm)))
	points(query_mc_d[,2], query_mc_d[,3], pch = 19, col = scales::alpha("red", 0.4), cex = 0.4)
	
	dev.off()
	

}

message("all done!")
