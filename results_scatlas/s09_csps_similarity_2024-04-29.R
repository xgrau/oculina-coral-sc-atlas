# libraries
suppressMessages(source("../scripts/Seurat_functions.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/Cross_species_functions.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/helper.R"))
graphics.off()
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))


# data index
sps_list = c("Ocupat","Ocuarb","Spin","Amil","Nvec","Xesp")
sps_refl = "Ocupat"

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)


# loop
for (spi in sps_refl) {

	# output
	out_fn = sprintf("results_metacell_%s_filt/csps/", spi)
	dir.create(out_fn, showWarnings = FALSE)
	
	# set working species
	if (spi == "Spio" | spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}


	# compare to loop
	for (spj in sps_list [ sps_list != spi ]) {
		
		# set working species
		if (spj == "Spio" | spj == "Spin") {
			spj_w = "Spis"
		} else {
			spj_w = spj
		}
		
		message(sprintf("csps | %s | compare to %s, load ICC...", spi, spj))
		icc_ecv = read.table(sprintf("%s/dat.icc.%s-%s.ec_scores.csv", out_fn, spi, spj), header = TRUE, sep = "\t")
		icc_ecv_v = dic_from_vecs(icc_ecv$sp1, icc_ecv$ec_value)
		icc_obj = readRDS(sprintf("%s/dat.icc.%s-%s.obj.rds", out_fn, spi, spj))

		message(sprintf("csps | %s | compare to %s, load footprints...", spi, spj))
		inp_fn_i = sprintf("results_metacell_%s_filt/", spi)
		inp_fn_j = sprintf("results_metacell_%s_filt/", spj)
		mc_fp_i = readRDS(sprintf("%s/dat.%s.expression.cts_fp.rds", inp_fn_i, spi))
		mc_fp_j = readRDS(sprintf("%s/dat.%s.expression.cts_fp.rds", inp_fn_j, spj))
		
		message(sprintf("csps | %s | compare to %s, load cell type annotations...", spi, spj))
		ctt_i = read.table(sprintf("%s/annot.%s.leiden.csv", inp_fn_i, spi), sep = "\t", header = TRUE, comment.char = "")
		ctt_j = read.table(sprintf("%s/annot.%s.leiden.csv", inp_fn_j, spj), sep = "\t", header = TRUE, comment.char = "")
		ctt_i$cell_type_species = sprintf("%s|%s", spi, ctt_i$cell_type)
		ctt_j$cell_type_species = sprintf("%s|%s", spj, ctt_j$cell_type)
		ctt_i_cts_col_v = dic_from_vecs(ctt_i$cell_type_species, ctt_i$color)
		ctt_j_cts_col_v = dic_from_vecs(ctt_j$cell_type_species, ctt_j$color)
		ctt_i_cts_col_vv = dic_from_vecs(ctt_i$cell_type, ctt_i$color)
		ctt_j_cts_col_vv = dic_from_vecs(ctt_j$cell_type, ctt_j$color)

		message(sprintf("csps | %s | compare to %s, conform matrices...", spi, spj))
		mc_fp_i_f = mc_fp_i [ icc_ecv$sp1, ]
		mc_fp_j_f = mc_fp_j [ icc_ecv$sp2, ]
		mc_fp_m_f = cbind(mc_fp_i_f, mc_fp_j_f)

		message(sprintf("csps | %s | compare to %s, get covariable genes...", spi, spj))
		genes_covariable = csps_select_covariable_genes(
			sp1_fp = mc_fp_i_f,
			sp2_fp = mc_fp_j_f, 
			merged = mc_fp_m_f,
			cross_fc_thrs = 2,
			cross_n = 1,
			method = "min_fc"	
		)
		
		# compute similarity based on weighted pearson
		message(sprintf("csps | %s | compare to %s, wpearson...", spi, spj))
		icc_out = csps_correlation_matrix_noobj(
			mm = mc_fp_m_f,
			m1 = mc_fp_i_f,
			m2 = mc_fp_j_f,
			prefix_sp1 = spi,
			prefix_sp2 = spj,
			quantile_discretisation = FALSE,
			use_var_genes = genes_covariable[[1]],
			gene_weights = icc_ecv_v [ genes_covariable[[1]] ],
			cor_method = "wpearson"
		)
		
		# plot
		pp1 = plot_complex_heatmap(
			icc_out$cor_matrix,
			cluster_row = FALSE,
			cluster_col = FALSE,
			color_min = 0,
			color_max = 0.7,
			fontsize = 10,
			categories_row = names(ctt_i_cts_col_v), 
			categories_col = names(ctt_j_cts_col_v), 
			colors_row = ctt_i_cts_col_v, 
			colors_col = ctt_j_cts_col_v, 
			color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
			cell_border = gpar(col = "white", lwd = 1, lty = 1),
			heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
		)
		
		message(sprintf("csps | %s | compare to %s, save...", spi, spj))
		saveRDS(icc_out, sprintf("%s/csps.%s-%s.cts.rds", out_fn, spi, spj))
		pdf(sprintf("%s/csps.%s-%s.cts.pdf", out_fn, spi, spj),width = 8+dim(icc_out$cor_matrix)[2]/20, height = 8+dim(icc_out$cor_matrix)[1]/20)
		print(pp1)
		dev.off()
		
		
		# plot expression of a duplicate set		
		message(sprintf("csps | %s | compare to %s, EC duplicate plots...", spi, spj))
		icc_dup_ecv = read.table(sprintf("%s/dat.icc.%s-%s.ec_scores.duplicates.csv", out_fn, spi, spj), header = TRUE, sep = "\t")
		cc_xt = table(icc_dup_ecv$cluster)
		cc_xt_f = names(cc_xt) [ cc_xt < 4 ]
		icc_dup_ecv_f = icc_dup_ecv [ icc_dup_ecv$cluster %in% cc_xt_f, ]
		
		# pdf
		pdf(sprintf("%s/csps.%s-%s.cts.ec_scores.duplicates.pdf", out_fn, spi, spj), width = 8, height = 20)
		for (clu in cc_xt_f) {
			
			icc_dup_ecv_ff = icc_dup_ecv_f [ icc_dup_ecv_f$cluster == clu, ]
			icc_dup_ecv_ff$ec_value [ is.na(icc_dup_ecv_ff$ec_value) ] = 0
			icc_dup_ecv_ff$ec_value [ is.na(icc_dup_ecv_ff$ec_value) ] = 0
			icc_dup_ecv_ff = icc_dup_ecv_ff [ order(icc_dup_ecv_ff$ec_value, decreasing = TRUE), ]
			layout(matrix(1:20, ncol = 2, byrow = TRUE))
			
			for (nn in 1:nrow(icc_dup_ecv_ff)) {
				
				ggi = icc_dup_ecv_ff[nn, "sp1"]
				ggj = icc_dup_ecv_ff[nn, "sp2"]
				ecv = icc_dup_ecv_ff[nn, "ec_value"]
				max_y = max(mc_fp_i[ggi,], mc_fp_j[ggj,], na.rm = TRUE)
				barplot(mc_fp_i[ggi,], las = 2, ylab = "Footprint", col = ctt_i_cts_col_vv[colnames(mc_fp_i)], main = sprintf("%s\n%s\nEC=%.2f", ggi, ogm_gv[ggi], ecv), cex.names = 0.6, cex.main = 0.8)
				abline(h = 1, lty = 2)
				barplot(mc_fp_j[ggj,], las = 2, ylab = "Footprint", col = ctt_j_cts_col_vv[colnames(mc_fp_j)], main = sprintf("%s\n%s\n", ggj, ogm_gv[ggj]), cex.names = 0.6, cex.main = 0.8)
				abline(h = 1, lty = 2)
				
			}
			
		}
		dev.off()

	}
	
}

message("All done!")
