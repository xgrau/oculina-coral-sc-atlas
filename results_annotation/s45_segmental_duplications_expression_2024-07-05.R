# load libs
suppressMessages(require("ape"))
suppressMessages(require("igraph"))
suppressMessages(require("rtracklayer"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
graphics.off()

## Scope of analysis ##

prob_thr = 0.95
distance_thr = 1e5

list_tests = list(
	"allsclerplusplus" = list(sps_list = c("Ocupat","Ocuarb","Gasp","Fspp","Spis","Pocdam","Amil","Adig","Gfas"), anc_list = c("Scleractinia","Hexacorallia","Anthozoa"))
)

sps_focus = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp")


## Define input ##

# input
daa_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"

# output
out_fn = "results_segmental_duplications_plus/"
dir.create(out_fn, showWarnings = FALSE)


# read gene classification, taxonomy, etc
message("posteriors | read orthology")
daa = read.table(daa_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
oga_gtv = dic_from_vecs(daa$transcript, daa$gene)

## Identify ancestral single-copy genes ##

for (nn in 1:length(list_tests)) {
	
	tid = names(list_tests)[nn]
	anc_list = list_tests[[nn]]$anc_list
	sps_list = list_tests[[nn]]$sps_list
	sps_focus_f = sps_focus [ sps_focus %in% sps_list ]
	message(sprintf("posteriors | %s | expression analysis in %i species in %i ancestors", tid, length(sps_focus_f), length(anc_list)))

	# reload
	segments_d = read.table(sprintf("%s/summary_segmental_duplications.%s.csv", out_fn, tid), sep = "\t", header = TRUE)
	segments_d = segments_d [ segments_d$species %in% sps_focus_f, ]
	segments_d$species = factor(segments_d$species, levels = sps_focus_f)
	message(sprintf("posteriors | %s | expression analysis of %i segment genes across all species (%i families)", tid, nrow(segments_d), length(unique(segments_d$orthogroup))))
	
	# create output folder for cds alignments
	dir.create(sprintf("%s/alignments_%s/", out_fn, tid), showWarnings = FALSE)
	
	for (spi in sps_focus_f) { 
		
		if (spi == "Spis") {
			spi_dd = "Spin"
		} else {
			spi_dd = spi
		}
		
		# log
		message(sprintf("posteriors | %s | %s: load data...", tid, spi))

		# load ct annot
		ctt = read.table(sprintf("../results_scatlas/results_metacell_%s_filt/annot.%s.leiden.csv", spi_dd, spi_dd), sep = "\t", header = TRUE, comment.char = "")
		ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
		ctt_lou_cts_v = dic_from_vecs(ctt$cluster, ctt$cell_type)
		ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
		ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	
		# metacell
		ctm = read.table(sprintf("../results_scatlas/results_metacell_%s_filt/annot.%s.mcs.csv", spi_dd, spi_dd), sep = "\t", header = TRUE, comment.char = "")
		ctm$metacell = factor(ctm$metacell, levels = unique(ctm$metacell))
		ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
		ctm_mcs_cts_v = dic_from_vecs(ctm$metacell, ctm$metacell_annotation)
	
		# load expression
		ctfp = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi_dd, spi_dd))
		ctfp_v = apply(ctfp, 1, function(vv) { zz = which.max(vv); sprintf("%s(%.2f)", names(zz), max(vv)) })
		ctfp_vn = apply(ctfp, 1, function(vv) { zz = which.max(vv); sprintf("%s", names(zz)) })
		# mcfp = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.mcs_fp.rds", spi_dd, spi_dd))
		# mcfp_v = apply(mcfp, 1, function(vv) { zz = which.max(vv); sprintf("%s(%.2f)", ctm_mcs_cts_v[names(zz)], max(vv)) })

		# load pairwise alignments
		message(sprintf("posteriors | %s | %s: sequence distance of paralogs, load...", tid, spi))
		bls_d = as.data.frame(data.table::fread(sprintf("../results_csps/results_blast/samap.%s-%s.blast.csv", spi, spi)))
		bls_d[,1] = gsub("_","-",bls_d[,1])
		bls_d[,2] = gsub("_","-",bls_d[,2])
		# bls_m = reshape(bls_d[,c(1:3)], idvar = "V1", timevar = "V2", direction = "wide")
		
		# # load CDS fasta
		# cds = seqinr::read.fasta(sprintf("../data/reference/%s_long.cds.fasta", spi))
		# names(cds) = gsub("_", "-", as.character(oga_gtv [ names(cds) ]))
		
		# get segments spi
		segments_i = segments_d [ segments_d$species == spi, ]
		message(sprintf("posteriors | %s | %s: expression analysis of %i segment genes (%i families)", tid, spi, nrow(segments_i), length(unique(segments_i$orthogroup))))
		
		# set gene-level table
		dupgenes_i = segments_i [ c("gene","orthogroup","duplicate_set","gene_name") ]
		dupgenes_i$has_redundant_copy = NA
		dupgenes_i$most_similar_copy  = NA
		dupgenes_i$corr_similar_copy  = NA
		dupgenes_i$pid_similar_copy   = NA
		dupgenes_i$status             = NA
		rownames(dupgenes_i) = gsub("_","-", dupgenes_i$gene)

		# validate thresholds with distribution of within-family correlations
		ggc_t_v = vector(mode = "numeric")
		for (ogi in unique(dupgenes_i$orthogroup)) {
			
			genes_ogi = dupgenes_i$gene [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = dupgenes_i$gene_name [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = paste(unique(sort(ognam_ogi)), collapse = " / ")
			genes_ogi_v = gsub("_","-", genes_ogi)
			genes_ogi_v = genes_ogi_v [ genes_ogi_v %in% rownames(ctfp) ]
			if (length(genes_ogi_v)>1) {
				gge = ctfp[genes_ogi_v,]
				ggc_t_i = cor(t(gge), method = "pearson", use = "everything")
				ggc_t_v = c(ggc_t_v, ggc_t_i [ upper.tri(ggc_t_i) | diag(ggc_t_i) ])
			}
		
		}
		ggc_t_v = ggc_t_v [ !is.na(ggc_t_v) ]
		
		# plot thresholds
		pdf(sprintf("%s/exp_segdups.%s.%s.thresholds.pdf", out_fn, tid, spi), height = 8, width = 8)
		hist(ggc_t_v, breaks = 80, xlim = c(-1,1))
		abline(v=0.6, col = "blue", lty = 2)
		abline(v=0.4, col = "darkblue", lty = 2)
		dev.off()

		# family-wise heatmaps
		pdf(sprintf("%s/exp_segdups.%s.%s.heatmaps.pdf", out_fn, tid, spi), height = 8, width = 16)
		for (ogi in unique(dupgenes_i$orthogroup)) {
			
			genes_ogi = dupgenes_i$gene [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = dupgenes_i$gene_name [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = paste(unique(sort(ognam_ogi)), collapse = " / ")
			genes_ogi_v = gsub("_","-", genes_ogi)
			
			# get expression matrix
			if (all(genes_ogi_v %in% rownames(ctfp))) {
				gge = ctfp[genes_ogi_v,]
			} else {
				genes_ogi_vv = genes_ogi_v [ genes_ogi_v %in% rownames(ctfp) ]
				genes_ogi_vn = setdiff(genes_ogi_v, rownames(ctfp))
				ctfp_pp = ctfp
				for (ggin in genes_ogi_vn) {
					ctfp_pp = rbind(ctfp, rep(NA, ncol(ctfp)))
				}
				rownames(ctfp_pp) = c(rownames(ctfp), genes_ogi_vn)
				gge = ctfp_pp[genes_ogi_v,]
			}
			
			# gene-gene correlation
			ggc_i = cor(t(gge), method = "pearson", use = "everything")
			ggc_i [ is.na(ggc_i) ] = 0
			
			# # gene-gene distance using ka/ks (omit, slow...)
			# # compute codon alignments clustalw
			# cds_i = cds [ genes_ogi_v ]
			# seqinr::write.fasta(cds_i, names = names(cds_i), file.out = sprintf("%s/alignments_%s/%s.%s.cds.fasta", out_fn, tid, ogi, spi))
			# invisible(seqinr::reverse.align(
			# 	nucl.file =  sprintf("%s/alignments_%s/%s.%s.cds.fasta", out_fn, tid, ogi, spi),
			# 	out.file =   sprintf("%s/alignments_%s/%s.%s.cds.l.fasta", out_fn, tid, ogi, spi),
			# 	align.prot = TRUE,
			# 	clustal.path = "/home/xavi/mambaforge/bin/clustalw"
			# ))
			# # read alignment
			# ali_i = suppressMessages(seqinr::read.alignment(sprintf("%s/alignments_%s/%s.%s.cds.l.fasta", out_fn, tid, ogi, spi), format = "fasta"))
			# ali_i_kaks = suppressMessages(seqinr::kaks(ali_i))
			# if (length(ali_i_kaks) > 1) {
			# 	ali_i_m_ka = as.matrix(ali_i_kaks$ka)
			# 	ali_i_m_ks = as.matrix(ali_i_kaks$ks)
			# 	ali_i_m_kaks = ali_i_m_ka / ali_i_m_ks
			# 	write.table(ali_i_m_ka, sprintf("%s/alignments_%s/%s.%s.kaks.ka.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
			# 	write.table(ali_i_m_ks, sprintf("%s/alignments_%s/%s.%s.kaks.ks.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
			# 	write.table(ali_i_m_kaks, sprintf("%s/alignments_%s/%s.%s.kaks.kaks.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
				
			# 	# similarity
			# 	ali_i_m_ks_sim = 1 - ali_i_m_ks / 10
			# 	rownames(ali_i_m_ks_sim) = rownames(ali_i_m_ks)
			# 	colnames(ali_i_m_ks_sim) = colnames(ali_i_m_ks)
			# } else {
			# 	ali_i_m_ks_sim = ggc_i
			# 	ali_i_m_ks_sim [ ali_i_m_ks_sim != 0 ] = 0
			# 	ali_i_m_ks = ali_i_m_ks_sim
			# }
			
			# # instead, set to zero...
			# ali_i_m_ks = matrix(nrow = length(genes_ogi_v), ncol = length(genes_ogi_v))
			# rownames(ali_i_m_ks) = genes_ogi_v
			# colnames(ali_i_m_ks) = genes_ogi_v
			# ali_i_m_ks_sim = ali_i_m_ks
			
			# instead, use pident...
			bls_d_ff = bls_d [ bls_d[,1] %in% genes_ogi_v & bls_d[,2] %in% genes_ogi_v, 1:3 ] 
			bls_d_ff = bls_d_ff [ bls_d_ff[,1] != bls_d_ff[,2],  ]
			bls_d_ff = bls_d_ff [ !duplicated(bls_d_ff[,1]),  ]
			bls_d_ff[,1] = factor(bls_d_ff[,1], levels = genes_ogi_v)
			bls_d_ff[,2] = factor(bls_d_ff[,2], levels = genes_ogi_v)
			bls_m_fm = tidyr::spread(bls_d_ff, key = V1, value = V3, drop = FALSE)
			rownames(bls_m_fm) = bls_m_fm [,1]
			bls_m_fm [,1] = NULL
			ali_i_m_ks = as.matrix(bls_m_fm)
			bls_m_fm [ is.na(bls_m_fm) ] = 0
			diag(ali_i_m_ks) = 100
			ali_i_m_ks_sim = ali_i_m_ks
			
			# # with test?
			# gge_t = t(gge)
			# ggp_i = sapply(1:ncol(gge_t), function(nn1) {
			# 	sapply(1:ncol(gge_t), function(nn2) {
			# 		tt = cor.test(gge_t[,nn1], gge_t[,nn2], method = "pearson", use = "everything")
			# 		tt$p.value
		  	# 	} )
			# } )
			# ggp_i_bool = ggp_i < 1e-3
			
			# tag copies according to similarity of expression
			ggc_i_f = ggc_i
			ali_i_m_ks_f = ali_i_m_ks
			diag(ggc_i_f) = NA
			diag(ali_i_m_ks_f) = NA
			ggc_has_divergent_copy = apply(ggc_i_f, 1, function(vv) min(vv, na.rm = TRUE)) <  0.4
			ggc_has_redundant_copy = apply(ggc_i_f, 1, function(vv) max(vv, na.rm = TRUE)) >= 0.6
			ggc_most_similar_copy  = apply(ggc_i_f, 1, function(vv) names(which.max(vv)) )
			ggc_corr_similar_copy  = apply(ggc_i_f, 1, function(vv) max(vv, na.rm = TRUE) )
			pid_similar_copy       = apply(ali_i_m_ks[rownames(ggc_i_f),colnames(ggc_i_f)], 1, function(vv) max(vv, na.rm = TRUE) )
			
			ggc_label = ifelse(any(ggc_has_redundant_copy), "redundant", "")
			ggc_label = paste(ggc_label, ifelse(any(ggc_has_divergent_copy), "divergent", ""), sep = "/")
			ggc_label = paste(ggc_label, ifelse(!all(genes_ogi_v %in% rownames(ctfp)), "nonexp", ""), sep = "/")
			ggc_label = gsub("^/", "", ggc_label)
			ggc_label = gsub("/$", "", ggc_label)
			ggc_label = gsub("/$", "", ggc_label)
			ggc_label [ ggc_label == "" ] = "undetermined"
			
			# save info per gene
			dupgenes_i[ rownames(ggc_i), "has_redundant_copy"] = ggc_has_redundant_copy
			dupgenes_i[ rownames(ggc_i), "most_similar_copy"]  = ggc_most_similar_copy
			dupgenes_i[ rownames(ggc_i), "corr_similar_copy"]  = ggc_corr_similar_copy
			dupgenes_i[ rownames(ggc_i), "pid_similar_copy"]    = pid_similar_copy
			dupgenes_i[ rownames(ggc_i), "status"] = ggc_label
			dupgenes_i[ rownames(ggc_i), "top_fc"] = ctfp_v [ rownames(ggc_i) ]
			# dupgenes_i[ rownames(ggc_i), "top_umifrac"] = mcuf_v [ rownames(ggc_i) ]

			# order genes		
			ggc_i_dis = as.dist(1 - cor(ggc_i))	
			ggc_i_dis [ is.na(ggc_i_dis) ] = 1
			ggc_i_clu = hclust(ggc_i_dis, method = "ward.D2")
			ggc_i_labels = ggc_i_clu$labels [ ggc_i_clu$order ]
			
			# plot gene-gene correlation
			pp1 = plot_complex_heatmap(
				ggc_i[ggc_i_labels,ggc_i_labels],
				name = "cor",
				color_max = 1,
				color_min = 0.5,
				cluster_row = FALSE,
				cluster_col = FALSE,
				use_raster = FALSE,
				color_mat = c("gray95", "deepskyblue","dodgerblue3","midnightblue"),
				cell_border = gpar(col = "white", lwd = 1, lty = 1),
				heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
			)
			pp1@matrix_param$width  = unit(3, "cm")
			pp1@matrix_param$height = unit(3, "cm")
			pp3 = plot_complex_heatmap(
				ali_i_m_ks_sim[ggc_i_labels,ggc_i_labels],
				name = "pident",
				color_max = 100,
				color_min = 50,
				cluster_row = FALSE,
				cluster_col = FALSE,
				use_raster = FALSE,
				color_mat = c("gray95", "slategray4"),
				cell_border = gpar(col = "white", lwd = 1, lty = 1),
				heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
			)
			pp3@matrix_param$width  = unit(3, "cm")
			pp3@matrix_param$height = unit(3, "cm")
			pp2 = plot_complex_heatmap(
				gge[ggc_i_labels,],
				name = "fp",
				color_mat = c("gray95","orange","orangered2","#520c52"),
				color_min = 1,
				cluster_row = FALSE,
				cluster_col = FALSE,
				use_raster = FALSE,
				categories_col = colnames(gge),
				colors_col = ctt_cts_col_v,
				title_col = sprintf("%s (%s)\n%s", ogi, ggc_label, ognam_ogi),
				color_max = 1.5,
				cell_border = gpar(col = "white", lwd = 1, lty = 1),
				heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
			)
			print(pp1 + pp3 + pp2)
			
		}
		dev.off()
		
		# aggregate status at OG level
		message(sprintf("posteriors | %s | %s: aggregate results per orthogroup...", tid, spi))
		dupgenes_i$gene = gsub("_", "-", dupgenes_i$gene)
		dupgenes_a = unique(dupgenes_i [ ,c("orthogroup","status") ])
		dupgenes_a$status = factor(dupgenes_a$status, levels = c("divergent","divergent/nonexp","redundant","redundant/nonexp","redundant/divergent","redundant/divergent/nonexp","undetermined"))
		
		# get gene ages per OG
		gea_d = read.table(sprintf("results_gene_family_evolution_anthozoa_plus/geneages.%s.posteriors.csv", spi), header = TRUE, sep = "\t")
		gea_v = dic_from_vecs(gea_d$orthogroup, gea_d$gain_node)
		dupgenes_a$gain_node_orthogroup = gea_v [ dupgenes_a$orthogroup ]
		write.table(dupgenes_a, file = sprintf("%s/exp_segdups.%s.%s.ogs.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = FALSE)

		# plot pies
		pdf(sprintf("%s/exp_segdups.%s.%s.pie.pdf", out_fn, tid, spi), height = 4, width = 4)
		
		# status by OG
		xtt = table(dupgenes_a$status)
		xtt_labels = sprintf("%s | n=%i", names(xtt), xtt)
		pie(xtt, col = viridisLite::turbo(nlevels(dupgenes_a$status), begin = 0.1, end = 0.8), label = xtt_labels, cex = 0.6)
		title(main = sprintf("%s families: %s", tid, spi), sub = sprintf("n=%i families", sum(xtt)))
		
		# status by best gene pair
		dupgenes_i_gene_pair_unique = data.frame(
			gene_pair = apply(dupgenes_i[,c("gene","most_similar_copy")], 1, function(v) { paste(sort(v), collapse = ",") }),
			status = "undetermined"
		)
		dupgenes_i_gene_pair_unique$status [ dupgenes_i$has_redundant_copy ] = "has_redundant_copy"
		dupgenes_i_gene_pair_unique$status [ !dupgenes_i$has_redundant_copy ] = "no_redundant_copy"
		dupgenes_i_gene_pair_unique$status = factor(dupgenes_i_gene_pair_unique$status, levels = c("no_redundant_copy","has_redundant_copy","undetermined"))
		dupgenes_i_gene_pair_unique = dupgenes_i_gene_pair_unique [ order(dupgenes_i_gene_pair_unique$gene_pair, dupgenes_i_gene_pair_unique$status), ]
		dupgenes_i_gene_pair_unique = dupgenes_i_gene_pair_unique [ !duplicated(dupgenes_i_gene_pair_unique$gene_pair), ]
		
		xtt = table(dupgenes_i_gene_pair_unique$status)
		xtt_labels = sprintf("%s | n=%i", names(xtt), xtt)
		pie(xtt, col = viridisLite::turbo(nlevels(dupgenes_i_gene_pair_unique$status), begin = 0.1, end = 0.8), label = xtt_labels, cex = 0.6)
		title(main = sprintf("%s genes: %s", tid, spi), sub = sprintf("n=%i unique gene pairs", sum(xtt)))
		plot(dupgenes_i$corr_similar_copy, dupgenes_i$pident, pch = 19, col= "blue", cex = 0.4)
		
		# status by top cell type, dup genes
		xtt = table(factor(gsub("\\(.*", "", dupgenes_i$top_fc), levels = levels(ctt$cell_type)))
		xtt_labels = sprintf("%s | n=%i", names(xtt), xtt)
		pie(xtt, col = ctt_cts_col_v[names(xtt)], label = xtt_labels, cex = 0.6)
		title(main = sprintf("%s top cell type per TD gene: %s", tid, spi), sub = sprintf("n=%i genes", sum(xtt)))
		
		# status by top cell type, all genes
		xtt = table(factor(ctfp_vn, levels = levels(ctt$cell_type)))
		xtt_labels = sprintf("%s | n=%i", names(xtt), xtt)
		pie(xtt, col = ctt_cts_col_v[names(xtt)], label = xtt_labels, cex = 0.6)
		title(main = sprintf("%s top cell type per gene (all): %s", tid, spi), sub = sprintf("n=%i genes", sum(xtt)))

		# with tests
		xtb = table(factor(ctfp_vn, levels = levels(ctt$cell_type)), factor(names(ctfp_vn) %in% dupgenes_i$gene, levels = c("TRUE","FALSE")))
		pvs = apply(xtb, 1, function(v) { cq = fisher.test(matrix(c(v, colSums(xtb) - v), nrow = 2), alternative = "two.side") })
		pva = p.adjust(sapply(pvs, function(vv) { vv$p.value }), method = "fdr")
		pvo = sapply(pvs, function(vv) { vv$estimate })
		pvob = sapply(pvs, function(vv) { sprintf("%.2f-%.2f", vv$conf.int[1], vv$conf.int[2]) })
		b = barplot(log2(pvo), col = "lightblue", las = 2, ylab = "log2(OR)", space = 0, ylim = c(-3,3), cex.names = 0.3)
		points(x = b, y = rep(-3, length(b)), pch = 22, bg = ctt_cts_col_v[rownames(xtb)], col = ctt_cts_col_v[rownames(xtb)])
		title(sub = "Fisher test of symbiotic sample enrichment, FDR-adjusted", adj = 0)
		text(x = b, y = -2, sprintf("p=%.1E | OR=%.2f(%s)", pva,pvo,pvob), cex = 0.3, col = "darkblue", pos = 3, srt = 90, font = 1)

		# with tests, aggregated
		ctfp_vn_agg = gsub("_.*", "", ctfp_vn)
		ctfp_vn_agg [ grepl("alga", ctfp_vn) ] = "gastrodermis_alga_hosting"
		ctfp_vn_agg [ grepl("muscle", ctfp_vn) ] = "gastrodermis_muscle"
		xtb = table(ctfp_vn_agg, factor(names(ctfp_vn) %in% dupgenes_i$gene, levels = c("TRUE","FALSE")))
		pvs = apply(xtb, 1, function(v) { cq = fisher.test(matrix(c(v, colSums(xtb) - v), nrow = 2), alternative = "two.side") })
		pva = p.adjust(sapply(pvs, function(vv) { vv$p.value }), method = "fdr")
		pvo = sapply(pvs, function(vv) { vv$estimate })
		pvob = sapply(pvs, function(vv) { sprintf("%.2f-%.2f", vv$conf.int[1], vv$conf.int[2]) })
		b = barplot(log2(pvo), col = "lightblue", las = 2, ylab = "log2(OR)", space = 0, ylim = c(-3,3), cex.names = 0.3)
		points(x = b, y = rep(-3, length(b)), pch = 22, bg = ctt_cts_col_v[rownames(xtb)], col = ctt_cts_col_v[rownames(xtb)])
		title(sub = "Fisher test of symbiotic sample enrichment, FDR-adjusted", adj = 0)
		text(x = b, y = -2, sprintf("p=%.1E | OR=%.2f(%s)", pva,pvo,pvob), cex = 0.3, col = "darkblue", pos = 3, srt = 90, font = 1)
		
		# close pie
		dev.off()
	
		# get alignment distance of paralogs
		message(sprintf("posteriors | %s | %s: sequence distance of paralogs, filter...", tid, spi))
		bls_d_f = bls_d [ bls_d[,1] %in%  c(dupgenes_i$gene, dupgenes_i$most_similar_copy) & bls_d[,2] %in%  c(dupgenes_i$gene, dupgenes_i$most_similar_copy), ]
		
		# pident of best paralog pairs (based on expression)
		message(sprintf("posteriors | %s | %s: sequence distance of paralogs, fetch...", tid, spi))
		dupgenes_i_pident = sapply(1:nrow(dupgenes_i), function(vv) { 
			if (dupgenes_i$gene[vv] %in% bls_d_f[,1] & dupgenes_i$most_similar_copy[vv] %in% bls_d_f[,2]){
				pid = bls_d_f[bls_d_f[,1] == dupgenes_i$gene[vv] & bls_d_f[,2] == dupgenes_i$most_similar_copy[vv], 3][1]
			} else {
				pid = NA
			}
		})
		dupgenes_i$pident = dupgenes_i_pident
		dupgenes_i_f = dupgenes_i [ !is.na(dupgenes_i$pident), ]

		# relationship between sequence distance and expression similarity		
		message(sprintf("posteriors | %s | %s: sequence distance of paralogs, plot...", tid, spi))
		pdf(sprintf("%s/exp_segdups.%s.%s.pidentcor.pdf", out_fn, tid, spi), height = 6, width = 6)

		# # plot, based on ks distance (alignments)
		# ks_bool = dupgenes_i_f$ks_similar_copy < 6
		# dupgenes_i_f$ks_category = NA
		# dupgenes_i_f$ks_category[ks_bool] = cut(dupgenes_i_f$ks_similar_copy[ks_bool], 5)
		# plot(dupgenes_i_f$ks_similar_copy[ks_bool], dupgenes_i_f$corr_similar_copy[ks_bool], pch = 19, cex = 0.6, col = "lightblue3", xlab = "ks distance", ylab = "expr corr most similar paralog", ylim = c(-1,1), las = 2)
		# boxplot(dupgenes_i_f$corr_similar_copy[ks_bool] ~ dupgenes_i_f$ks_category[ks_bool], pch = 19, cex = 0.6, col = "lightblue3", las = 2, ylab = "expr corr most similar paralog", xlab = "", ylim = c(-1,1))
		# tt = cor.test(dupgenes_i_f$ks_similar_copy[ks_bool], dupgenes_i_f$corr_similar_copy[ks_bool], method = "spearman")
		# title(sub = sprintf("spearman rho = %.3f, p=%.2E, n=%i", tt$estimate, tt$p.value, sum(ks_bool)), main = spi)

		# plot, based on percentage identity (blast)
		dupgenes_i_f$pident_category = cut(dupgenes_i_f$pident, 5)
		plot(dupgenes_i_f$pident, dupgenes_i_f$corr_similar_copy, pch = 19, cex = 0.6, col = "lightblue3", xlab = "% identity most similar paralog", ylab = "expr corr most similar paralog", ylim = c(-1,1), las = 2)
		boxplot(dupgenes_i_f$corr_similar_copy ~ dupgenes_i_f$pident_category, pch = 19, cex = 0.6, col = "lightblue3", las = 2, ylab = "expr corr most similar paralog", xlab = "", ylim = c(-1,1))
		tt = cor.test(dupgenes_i_f$pident, dupgenes_i_f$corr_similar_copy, method = "spearman")
		title(sub = sprintf("spearman rho = %.3f, p=%.2E, n=%i", tt$estimate, tt$p.value, nrow(dupgenes_i_f)), main = spi)

		# close plot
		dev.off()
		
		# output
		message(sprintf("posteriors | %s | %s: save...", tid, spi))
		write.table(dupgenes_i, file = sprintf("%s/exp_segdups.%s.%s.genes.csv", out_fn, tid, spi), sep = "\t", quote = FALSE, row.names = FALSE)

		# # get ks distance of paralogs
		# kst_d = read.table(sprintf("results_wgd/paralog_distributions/wgd_%s/%s.ks.tsv", spi, spi), header = TRUE, sep = "\t")
		# kst_d_f = kst_d [ ,c("Paralog1","Paralog2","Ka","Ks","Distance","AlignmentIdentity") ]
		# kst_d_f$Paralog1 = gsub("_", "-", oga_gtv [ kst_d_f$Paralog1 ])
		# kst_d_f$Paralog2 = gsub("_", "-", oga_gtv [ kst_d_f$Paralog2 ])
		# kst_d_f = kst_d_f [ kst_d_f$Paralog1 %in% c(dupgenes_i$gene, dupgenes_i$most_similar_copy) | kst_d_f$Paralog2 %in% c(dupgenes_i$gene, dupgenes_i$most_similar_copy), ]
		# # pairwise table...
		# kst_m_f = xtabs(kst_d_f$AlignmentIdentity ~ kst_d_f$Paralog1 + kst_d_f$Paralog2)
		# kst_m_f = table_to_matrix(kst_m_f)
		
		# dupgenes_i$ks_distance = sapply(1:nrow(dupgenes_i), function(vv) { 
		# 	if (dupgenes_i$gene[vv] %in% rownames(kst_m_f) & dupgenes_i$most_similar_copy[vv] %in% colnames(kst_m_f)){
		# 		ks = kst_m_f[ dupgenes_i$gene[vv],dupgenes_i$most_similar_copy[vv] ]
		# 	} else {
		# 		ks = NA
		# 	}
		# })
		# dupgenes_i_f = dupgenes_i [ !is.na(dupgenes_i$ks_distance), ]
		# pdf("hola.pdf", w=7, h=7)
		# plot(dupgenes_i_f$ks_distance, dupgenes_i_f$corr_similar_copy, pch = 19, cex = 0.6, col = "blue")
		# dev.off()
		# cor.test(dupgenes_i_f$ks_distance, dupgenes_i_f$corr_similar_copy, method = "spearman")
		
	}
	
	
}