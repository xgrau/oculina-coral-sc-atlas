# libraries
require("topGO")
require("VennDiagram")
require("grid")

# Load pfam annotations as list, for functional enrichments
#' 
#' @param pfam_architecture_file tab-separated file where first column is gene name and second column is a string of domain architectures for that gene (by default, separated by whitespace)
#' @param do_remove_empty remove genes without annotation (default: TRUE)
#' @param do_unique remove redundant annotations for each gene (default: TRUE)
#' @param architecture_sep character separating the annotations in each row (by default, whitespace) (default: whitespace)
#' 
gsa_enrichment_load_pfam_list = function(pfam_architecture_file, do_unique = TRUE, do_remove_empty = TRUE, architecture_sep = " ") {
	
	# load pfam species i
	pfm_i = read.table(pfam_architecture_file, col.names = c("gene","architecture"), sep = "\t")
	if (do_unique) {
		pfm_i_l = lapply(pfm_i$architecture, function(v) unique(unlist(strsplit(v, split = architecture_sep))) )
		names(pfm_i_l) = pfm_i$gene
	} else {
		pfm_i_l = lapply(pfm_i$architecture, function(v) unlist(strsplit(v, split = architecture_sep)) )
		names(pfm_i_l) = pfm_i$gene
	}
	
	if (do_remove_empty) {
		pfm_i_l = pfm_i_l [ lengths(pfm_i_l) > 0 ]
	}
	
	# output
	return(pfm_i_l)
	
}

# Functional enrichment of gene lists with a hypergeometric test
#' 
#' @param annotation named list, where each entry is a gene, and each entry contains a list of annotations mapped to this gene (e.g. a unique list of pfam domains in the gene)
#' @param genes_fg vector of genes in the foreground (list of genes of interest)
#' @param genes_bg vector of genes in the background (this defaults to all genes in the annotation object, but you can specify a narrower background if needed)
#' @param do_print_output,name_fg,top_markers control output plot
#' @param p_adj_method pvalue adjustment method for p.adj function (default is BH)
#' 
#' @return data.frame with annotation enrichments
#'
gsa_enrichment_hypergeometric = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		do_print_output = TRUE,
		output_prefix = "enrichment_output",
		name_fg = "fg",
		top_markers = 30,
		p_adj_method = "BH") {
	
	if (is.null(genes_bg)) {
		genes_bg = names(annotation)
	}
	
	# subset list of annotations to genes in foreground and background
	annotation = annotation [ names(annotation) %in% c(genes_fg, genes_bg) ]
	
	# count distributions in background and foreground list
	dist_fg = as.data.frame(table(unlist(annotation [ genes_fg ])))
	dist_bg = as.data.frame(table(unlist(annotation)))
	
	if (nrow(dist_fg) > 0) {
		
		# prepare table
		tab_in_all  = merge(dist_fg, dist_bg, by.x = "Var1", by.y = "Var1", all.x = TRUE, all.y = FALSE)
		tab_in_all  = tab_in_all[apply(tab_in_all!=0, 1, all),]
		colnames(tab_in_all) = c("annot","freq_in_fg","freq_in_bg")
		tab_in_all$total_annot_in_fg = length(unlist(annotation [ genes_fg ]))
		tab_in_all$total_annot_in_bg = length(unlist(annotation))
		
		# hypergeometric test
		tab_in_all$pval = phyper(
			tab_in_all$freq_in_fg - 1,
			tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_bg - tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_fg,
			lower.tail = FALSE)
		
		# pvalue adjustment
		tab_in_all$pval_adj = p.adjust(tab_in_all$pval, method = p_adj_method)
		
		# sort result by pvalue
		tab_in_all = tab_in_all [ order(tab_in_all$pval_adj) , ]
		
		# Output 
		if(do_print_output) {
			
			pdf(paste(output_prefix,".",name_fg,".hypergeo.pdf",sep=""), height = 4 + length(top_markers / 10), width = 8)
			# par(mar=c(5,10,5,2))
			layout(matrix(1:3, ncol = 3))
			
			# padding plot			
			plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
			
			# pval plot
			ploti = barplot(
				height = -rev(head(log10(tab_in_all$pval_adj), top_markers)),
				names.arg = rev(head(tab_in_all$annot, top_markers)),
				xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = name_fg,
				sub = sprintf("n=%i/%i input genes with annotations", nrow(dist_fg), length(genes_fg)),
				xlab = "-log10(p)"
			)
			abline(v=log10(0.01),lty=2,lwd=0.5,col="indianred2")
			abline(v=log10(0.05),lty=2,lwd=0.5,col="indianred2")
			text(
				x= 0, ploti,
				labels = sprintf("p=%.1E | n=%i", rev(head(tab_in_all$pval_adj, top_markers)) , rev(head(tab_in_all$freq_in_fg, top_markers))),
				col = "gray20", pos = 4)
			
			# fg/bg fraction
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)),
				horiz = TRUE, las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction"
			)
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers)),
				horiz = TRUE, las = 1,col = "azure4", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction",
				add = TRUE, xaxt = "n"
			)
			text(
				x= 0, ploti,
				labels = sprintf("fg=%.2f | bg=%.2f", rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) , rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers))),
				col = "gray20", pos = 4)
			# ploti = barplot(
			# 	height = rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) / rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers)),
			# 	names.arg = rev(head(tab_in_all$annot, top_markers)),
			# 	horiz = TRUE, las = 1, col = "azure3", border = NA, ylim = c(0, top_markers + 10), log= "x",
			# 	xlim = c(1, 100),
			# 	main = "FC",
			# 	xlab = "FC",
			# )
			# text(
			# 	x=1, ploti,
			# 	labels = sprintf("fc=%.2f", rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) / rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers))),
			# 	col = "gray20", pos = 4)
			
			
			dev.off()
			
			# table
			write.table(
				tab_in_all,
				file = sprintf("%s.%s.hypergeo.csv", output_prefix, name_fg),
				row.names = FALSE, sep="\t", quote = FALSE)
		}
		
	} else {
		print("skip, no annotations in interest list!")
		tab_in_all = NULL
	}
	return(tab_in_all)
	
}

# Functional enrichment of gene lists with a hypergeometric test
#' 
#' @param annotation named list, where each entry is a gene, and each entry contains a list of annotations mapped to this gene (e.g. a unique list of pfam domains in the gene)
#' @param genes_fg vector of genes in the foreground (list of genes of interest)
#' @param genes_bg vector of genes in the background (this defaults to all genes in the annotation object, but you can specify a narrower background if needed)
#' @param do_print_output,name_fg,top_markers control output plot
#' @param p_adj_method pvalue adjustment method for p.adj function (default is BH)
#' 
#' @return data.frame with annotation enrichments
#'
gsa_enrichment_hypergeometric_v2 = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		do_print_output = TRUE,
		out_fn = "enrichment_output",
		name_fg = "fg",
		top_markers = 30,
		p_adj_method = "BH") {
	
	if (is.null(genes_bg)) {
		genes_bg = names(annotation)
	}
	
	# subset list of annotations to genes in foreground and background
	annotation = annotation [ names(annotation) %in% c(genes_fg, genes_bg) ]
	
	# count distributions in background and foreground list
	dist_fg = as.data.frame(table(unlist(annotation [ genes_fg ])))
	dist_bg = as.data.frame(table(unlist(annotation)))
	
	if (nrow(dist_fg) > 0) {
		
		# prepare table
		tab_in_all  = merge(dist_fg, dist_bg, by.x = "Var1", by.y = "Var1", all.x = TRUE, all.y = FALSE)
		tab_in_all  = tab_in_all[apply(tab_in_all!=0, 1, all),]
		colnames(tab_in_all) = c("annot","freq_in_fg","freq_in_bg")
		tab_in_all$total_annot_in_fg = length(unlist(annotation [ genes_fg ]))
		tab_in_all$total_annot_in_bg = length(unlist(annotation))
		tab_in_all$fraction_in_fg = tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg
		tab_in_all$fraction_in_bg = tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg
		
		# hypergeometric test
		tab_in_all$pval = phyper(
			tab_in_all$freq_in_fg - 1,
			tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_bg - tab_in_all$freq_in_bg,
			tab_in_all$total_annot_in_fg,
			lower.tail = FALSE)
		
		# pvalue adjustment
		tab_in_all$pval_adj = p.adjust(tab_in_all$pval, method = p_adj_method)
		
		# sort result by pvalue
		tab_in_all = tab_in_all [ order(tab_in_all$pval_adj) , ]
		
		# Output 
		if(do_print_output) {
			
			pdf(sprintf("%s.pdf", gsub(".csv$", "", out_fn)), height = 4 + length(top_markers / 10), width = 8)
			# par(mar=c(5,10,5,2))
			layout(matrix(1:3, ncol = 3))
			
			# padding plot			
			plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
			
			# pval plot
			ploti = barplot(
				height = -rev(head(log10(tab_in_all$pval_adj), top_markers)),
				names.arg = rev(head(tab_in_all$annot, top_markers)),
				xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = name_fg,
				sub = sprintf("n=%i/%i input genes with annotations", nrow(dist_fg), length(genes_fg)),
				xlab = "-log10(p)"
			)
			abline(v=log10(0.01),lty=2,lwd=0.5,col="indianred2")
			abline(v=log10(0.05),lty=2,lwd=0.5,col="indianred2")
			text(
				x= 0, ploti,
				labels = sprintf("p=%.1E | n=%i", rev(head(tab_in_all$pval_adj, top_markers)) , rev(head(tab_in_all$freq_in_fg, top_markers))),
				col = "gray20", pos = 4)
			
			# fg/bg fraction
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)),
				horiz = TRUE, las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction"
			)
			ploti = barplot(
				height = rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers)),
				horiz = TRUE, las = 1,col = "azure4", border = NA, ylim = c(0, top_markers + 10),
				main = "fraction genes in fg and bg",
				xlab = "fraction",
				add = TRUE, xaxt = "n"
			)
			text(
				x= 0, ploti,
				labels = sprintf("fg=%.2f | bg=%.2f", rev(head(tab_in_all$freq_in_fg / tab_in_all$total_annot_in_fg, top_markers)) , rev(head(tab_in_all$freq_in_bg / tab_in_all$total_annot_in_bg, top_markers))),
				col = "gray20", pos = 4)
			
			dev.off()
			
			# table
			write.table(
				tab_in_all,
				file = out_fn,
				row.names = FALSE, sep="\t", quote = FALSE)
		}
		
	} else {
		message("skip, no annotations in fg list!")
		tab_in_all = NULL
	}
	
	# mappings table
	annotation_f = annotation [ genes_fg [ genes_fg %in% names(annotation)] ]
	if (length(annotation_f) > 0) {
		annotation_d = data.frame(
			gene = names(annotation_f),
			annot = unlist(sapply(1:length(annotation_f), function(nn) { paste(annotation_f[[nn]], collapse = " ") } ))
		)
	} else {
		annotation_d = data.frame(matrix(nrow = 0, ncol = 2))
		colnames(annotation_d) = c("gene","annot")
	}
	
	# return
	return(list(enrichment_table = tab_in_all, gene_mapping = annotation_d))
	
}

### Topgo ###

gsa_topgo_load_emapper = function(
	emapper_fn,
	sep_col = "\t",
	sep_gos = ",",
	index_col_GOs = 10,
	index_col_gen = 1,
	make_unique = TRUE
) {
	
	gos_i = read.table(emapper_fn, sep = sep_col)
	gos_i = data.frame("gene" = gos_i[,index_col_gen], "GO" = gos_i[,index_col_GOs])
	if (make_unique) {
		gos_i_l = lapply(gos_i$GO, function(v) unique(unlist(strsplit(v, split = sep_gos))) )
	} else {
		gos_i_l = lapply(gos_i$GO, function(v) unlist(strsplit(v, split = sep_gos)) )
	}
	names(gos_i_l) = gos_i$gene
	return(gos_i_l)
	
}


gsa_topgo_enrichment  = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		do_print_output = TRUE,
		output_prefix = "enrichment_output",
		name_fg = "fg",
		ontologyset = c("BP","MF","CC"),
		tg_test = "fisher",
		tg_algorithm = "elim",
		top_markers = 30,
		nodesize = 10,
		printfile = TRUE,
		p_adj_method="BH") {
	
	# Input 
	genes_fg = unique(genes_fg)
	
	
	if (!is.null(genes_bg)) {
		genes_bg = genes_bg [ genes_bg %in% names(annotation) ]
	} else  {
		genes_bg = names(annotation)
	}
	genes_fg_ix = factor(as.integer(genes_bg %in% genes_fg))
	names(genes_fg_ix) = genes_bg
	
	# shortened go mappings without empty transcripts
	gomap_nonempty = annotation[lapply(annotation,length)>0]
	
	if (do_print_output) {
		
		pdf(paste(output_prefix,".", name_fg,".topgo.pdf", sep=""), height = (4 + length(top_markers / 10)) * 3, width = 8)
		# par(mar=c(5,10,5,2))
		layout(matrix(1:9, ncol = 3, byrow = TRUE))
		
	}
	
	# test each go domain
	topgo_tau_tot = data.frame()
	
	if (length(genes_fg[genes_fg %in% names(gomap_nonempty)])>1) {
		
		for (ontologyseti in ontologyset) {
			
			# topgo setup
			GOdata = suppressMessages(new("topGOdata", ontology = ontologyseti, allGenes = genes_fg_ix, annot = annFUN.gene2GO, gene2GO = annotation))
			num_interest_feasible = sum(GOdata@feasible & genes_bg %in% genes_fg)
			
			# topGO test 
			topgo_res = suppressMessages(topGO::runTest(GOdata, algorithm = tg_algorithm, statistic = tg_test))
			topgo_tau = suppressMessages(topGO::GenTable(GOdata, pval_test = topgo_res, orderBy = "pval_test",  topNodes = length(usedGO(object = GOdata))))
			topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
			topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
			topgo_tau$ontology = ontologyseti
			topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
			
			if (do_print_output) {
				
				# padding plot
				plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
				
				# pval plot
				ploti = barplot(
					height = -rev(head(log(topgo_tau$pval_test,10), top_markers)),
					names.arg = rev(head(paste(topgo_tau$Term, topgo_tau$GO.ID), top_markers)),
					xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = sprintf("GO:%s\n%s", ontologyseti, name_fg),
					sub = sprintf("n=%i/%i input genes with annotations", num_interest_feasible, length(genes_fg)),
					xlab="-log(p)")
				abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
				abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
				text(
					x= 0, ploti,
					labels = sprintf("p=%.1E | n=%i", rev(head(topgo_tau$pval_test, top_markers)), rev(head(topgo_tau$Significant, top_markers))),
					col = "gray20", pos = 4)
				
				# fraction significant
				ploti = barplot(
					height = rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)),
					horiz = TRUE, las = 1, col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = "fraction genes in fg and expected value",
					xlab = "fraction"
				)
				ploti = barplot(
					height = rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)),
					horiz = TRUE, las = 1, col = "azure4", border = NA, ylim = c(0, top_markers + 10),
					add = TRUE, xaxt = "n"
				)
				text(
					x= 0, ploti,
					labels = sprintf("fg=%.2f | bg=%.2f", rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)), rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)) ),
					col = "gray20", pos = 4)
				
				
			}
			
		}
		
	} else {
		print("skip, no annotations in interest list!")
	}
	
	if (do_print_output) {
		# table
		write.table(
			topgo_tau_tot [ topgo_tau_tot$Significant > 0, ],
			file = sprintf("%s.%s.topgo.csv", output_prefix, name_fg),
			row.names = FALSE, sep="\t", quote = FALSE)
		dev.off()
	}
	
	# return
	return(topgo_tau_tot)
	
}

# # topgo, prefetch go_annot_object_list for faster enrichment analysis
# gsa_topgo_prefetch_object_v2 = function(
# 	annotation,
# 	genes_annot = NULL,
# 	ontologyset = c("BP","MF","CC")) {
	
# 	if (!is.null(genes_annot)) {
# 		genes_annot = genes_annot [ genes_annot %in% names(annotation) ]
# 	} else  {
# 		genes_annot = names(annotation)
# 	}
# 	genes_fg_ix = factor(as.integer(genes_bg %in% genes_fg), levels = c(0,1))
# 	names(genes_fg_ix) = genes_bg


# 	for (ontologyseti in ontologyset) {
# 		go_annot_object = suppressMessages(new("topGOdata", ontology = ontologyseti, allGenes = genes_fg_ix, annot = annFUN.gene2GO, gene2GO = annotation))
# 	}
	
	
# 	return()
	
# }

# topgo enrichments, new version
gsa_topgo_enrichment_v2  = function(
		annotation,
		genes_fg,
		genes_bg = NULL,
		go_annot_object_list = NULL,
		do_print_output = TRUE,
		out_fn = "enrichment_output",
		name_fg = "fg",
		ontologyset = c("BP","MF","CC"),
		tg_test = "fisher",
		tg_algorithm = "elim",
		top_markers = 30,
		nodesize = 10,
		printfile = TRUE,
		verbose = FALSE,
		p_adj_method="BH") {
	
	# Input 
	genes_fg = unique(genes_fg)
	
	# browser()
	
	if (!is.null(genes_bg)) {
		genes_bg = genes_bg [ genes_bg %in% names(annotation) ]
	} else  {
		genes_bg = names(annotation)
	}
	genes_fg_ix = factor(as.integer(genes_bg %in% genes_fg), levels = c(0,1))
	names(genes_fg_ix) = genes_bg
	
	# shortened go mappings without empty transcripts
	gomap_nonempty = annotation[lapply(annotation,length)>0]
	
	if (do_print_output) {
		
		pdf(sprintf("%s.pdf", gsub(".csv$", "", out_fn)), height = (4 + length(top_markers / 10)) * 3, width = 8)
		layout(matrix(1:9, ncol = 3, byrow = TRUE))
		
	}
	
	# test each go domain
	topgo_tau_tot = data.frame()
	fgi_genes_gos_d = data.frame()
	for (ontologyseti in ontologyset) {
		
		# topgo setup
		if (verbose) { message(sprintf("go enrichments | %s setup...", ontologyseti)) }
		if (is.null(go_annot_object_list[[ontologyseti]])) {
			go_annot_object = suppressMessages(new("topGOdata", ontology = ontologyseti, allGenes = genes_fg_ix, annot = annFUN.gene2GO, gene2GO = annotation))
		} else {
			go_annot_object = go_annot_object_list[[ontologyseti]]
		}
		num_interest_feasible = sum(go_annot_object@feasible & genes_bg %in% genes_fg)

		# browser()
		if (length(genes_fg[genes_fg %in% names(gomap_nonempty)])>1 & num_interest_feasible > 0) {

			# store GOs mapped to each gene
			if (verbose) { message(sprintf("go enrichments | %s genes per term...", ontologyseti)) }
			fgi_genes_gos_i = topGO::genesInTerm(go_annot_object)
			fgi_genes_gos_d = rbind(
				fgi_genes_gos_d,
				data.frame(
					GO = unlist(as.list(sapply(1:length(fgi_genes_gos_i), function(vv) { rep(names(fgi_genes_gos_i)[vv], lengths(fgi_genes_gos_i)[vv])  } ))),
					gene = as.character(unlist(fgi_genes_gos_i)),
					ontology = ontologyseti
				)
			)

			# topGO test 
			if (verbose) { message(sprintf("go enrichments | %s test...", ontologyseti)) }
			topgo_res = suppressMessages(topGO::runTest(go_annot_object, algorithm = tg_algorithm, statistic = tg_test))
			topgo_tau = suppressMessages(topGO::GenTable(go_annot_object, pval_test = topgo_res, orderBy = "pval_test",  topNodes = length(usedGO(object = go_annot_object))))
			topgo_tau$pval_test = as.numeric(topgo_tau$pval_test)
			topgo_tau$pval_adj  = p.adjust(topgo_tau$pval_test, method=p_adj_method)
			topgo_tau$ontology = ontologyseti
			topgo_tau_tot = rbind(topgo_tau_tot,topgo_tau)
			
			if (do_print_output) {
				
				# padding plot
				if (verbose) { message(sprintf("go enrichments | %s plot...", ontologyseti)) }
				plot(NA,NA,xaxt="n",yaxt="n",xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",bty="n")
				
				# pval plot
				ploti = barplot(
					height = -rev(head(log(topgo_tau$pval_test,10), top_markers)),
					names.arg = rev(head(paste(topgo_tau$Term, topgo_tau$GO.ID), top_markers)),
					xlim = c(0,5), horiz = TRUE,las = 1,col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = sprintf("GO:%s\n%s", ontologyseti, name_fg),
					sub = sprintf("n=%i/%i input genes with annotations", num_interest_feasible, length(genes_fg)),
					xlab="-log(p)")
				abline(v=log(0.01,10),lty=2,lwd=0.5,col="pink")
				abline(v=log(0.05,10),lty=2,lwd=0.5,col="pink")
				text(
					x= 0, ploti,
					labels = sprintf("p=%.1E | n=%i", rev(head(topgo_tau$pval_test, top_markers)), rev(head(topgo_tau$Significant, top_markers))),
					col = "gray20", pos = 4)
				
				# fraction significant
				ploti = barplot(
					height = rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)),
					horiz = TRUE, las = 1, col = "azure3", border = NA, ylim = c(0, top_markers + 10),
					main = "fraction genes in fg and expected value",
					xlab = "fraction"
				)
				ploti = barplot(
					height = rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)),
					horiz = TRUE, las = 1, col = "azure4", border = NA, ylim = c(0, top_markers + 10),
					add = TRUE, xaxt = "n"
				)
				text(
					x= 0, ploti,
					labels = sprintf("fg=%.2f | bg=%.2f", rev(head(topgo_tau$Significant / num_interest_feasible, top_markers)), rev(head(topgo_tau$Expected / num_interest_feasible, top_markers)) ),
					col = "gray20", pos = 4)
				
				
			}
			
		}
		
	}
	# ignore genes not in fg for mapping
	fgi_genes_gos_d_f = fgi_genes_gos_d [ fgi_genes_gos_d$gene %in% genes_fg, ]
		
	if (do_print_output) {
		# table
		write.table(
			topgo_tau_tot [ topgo_tau_tot$Significant > 0, ],
			file = out_fn,
			row.names = FALSE, sep="\t", quote = FALSE)
		dev.off()
	}
	
	# return
	return(list(enrichment_table = topgo_tau_tot, gene_mapping = fgi_genes_gos_d_f))
	
}

### Merge enrichment files ###

gsa_merge_tables_v2 = function(
	main_enrichment,
	main_mapping,
	pfam_enrichment,
	pfam_mapping,
	output_fn = "summarised_enrichments.tsv",
	main_gene = "gene",
	main_annot = "GO",
	main_enrichment_fields = c("GO.ID","Term","pval_test","Significant"),
	main_mapping_fields = c("GO","gene","ontology"),
	main_sort = "pval_test",
	pfam_gene = "gene",
	pfam_annot = "annot",
	pfam_freq  = "freq_in_fg",
	pfam_sort = "pval_adj",
	extra_annot_dict = NULL
	) {
	
	if (nrow(main_enrichment) == 0 & ncol(main_enrichment) == 0) {
		main_enrichment = data.frame(matrix(ncol = length(main_enrichment_fields), nrow = 0))
		colnames(main_enrichment) = main_enrichment_fields
	}
	
	if (nrow(main_mapping) == 0 & ncol(main_mapping) == 0) {
		main_mapping = data.frame(matrix(ncol = 3, nrow = 0))
		colnames(main_mapping) = main_mapping_fields
	}
	
	# in mappings, only keep annotations that are present in the enrichment tables
	main_enrichment_f = main_enrichment [ main_enrichment$Significant > 0, main_enrichment_fields ]
	main_mapping_f = main_mapping [ main_mapping[,main_annot] %in% main_enrichment_f[,main_enrichment_fields[1]], main_mapping_fields ]
	pfam_enrichment_f = pfam_enrichment [ pfam_enrichment[,pfam_freq] > 0, ]
	pfam_mapping_f = pfam_mapping #[ pfam_mapping[,pfam_annot] %in% as.character(pfam_enrichment_f[,pfam_annot]), ]

	# rename columns in first enrichment
	colnames(main_enrichment_f) = c("GO","GO_term","GO_pval","GO_fg_freq")
	
	# vector of pvals per secondary annotation
	pfam_enrichment_p_v = pfam_enrichment[,pfam_sort]
	names(pfam_enrichment_p_v) = pfam_enrichment[,pfam_annot]

	# vector of num per secondary annotation
	pfam_enrichment_n_v = pfam_enrichment[,pfam_freq]
	names(pfam_enrichment_n_v) = pfam_enrichment[,pfam_annot]
	
	# dictionary of pfam architectures
	pfam_mapping_f_v = stringr::str_split(pfam_mapping_f[,pfam_annot], pattern = " ")
	names(pfam_mapping_f_v) = pfam_mapping_f[,pfam_gene]
	
	# add counts and pval of best domain in each architecture
	if (nrow(pfam_mapping_f) != 0) {
		pfam_mapping_f_unfolded = data.frame(
			gene = as.vector(unlist(sapply(1:nrow(pfam_mapping_f), function(nn)  rep(pfam_mapping_f[nn,pfam_gene], stringr::str_count(pfam_mapping_f[nn,pfam_annot], pattern = "\\S+")) ))),
			pfam = as.vector(unlist(sapply(pfam_mapping_f[,pfam_annot], function(vv) unique(unlist(stringr::str_split(vv, pattern = " ")))))),
			pfam_architecture = as.vector(unlist(sapply(1:nrow(pfam_mapping_f), function(nn)  rep(pfam_mapping_f[nn,pfam_annot], stringr::str_count(pfam_mapping_f[nn,pfam_annot], pattern = "\\S+")) )))
		)
		pfam_mapping_f_unfolded$pfam_top_pvalue = pfam_enrichment_p_v [ pfam_mapping_f_unfolded$pfam ]
		pfam_mapping_f_unfolded$pfam_fg_freq    = pfam_enrichment_n_v [ pfam_mapping_f_unfolded$pfam ]
		pfam_mapping_f_unfolded$list_pfam_annot = pfam_mapping_f_v [ pfam_mapping_f_unfolded$gene ]
		pfam_mapping_f_unfolded$pfam_pvalue_string = sapply(pfam_mapping_f_unfolded$list_pfam_annot, function(vv) {
			tpv = sort(pfam_enrichment_p_v[vv])
			tpv = tpv [ tpv < 0.05 ]
			tpn = names(tpv)
			tps = sprintf("%s(%.1E)", tpn, tpv)
			tps = paste(tps, collapse = ";")
			}
		)
		pfam_mapping_f_unfolded$list_pfam_annot = NULL

		# sort by pfam value
		pfam_mapping_f_unfolded = pfam_mapping_f_unfolded [ order(pfam_mapping_f_unfolded$gene, pfam_mapping_f_unfolded$pfam_top_pvalue), ]
		pfam_mapping_f_unfolded = pfam_mapping_f_unfolded [ !duplicated(pfam_mapping_f_unfolded$gene), ]
	} else {
		pfam_mapping_f_unfolded = data.frame(gene = NA, pfam = NA, pfam_architecture = NA, pfam_top_pvalue = NA, pfam_fg_freq = NA, pfam_pvalue_string = NA)
	}
	
		
	# merge mappings of first and secondary enrichments
	gme = merge(main_mapping_f, pfam_mapping_f_unfolded, by.x = main_gene, by.y = pfam_gene, all.x = TRUE, all.y = TRUE)

	# add test info for first enrichment
	gme = merge(gme, main_enrichment_f, by.x = main_annot, by.y = main_annot, all.x = TRUE, all.y = TRUE)
	# gme = head(gme, 1e4)
	
	# fix column names
	gmo = gme
	colnames(gmo) = c("GO","gene","GO_ontology","pfam_top_domain","pfam_architecture","pfam_top_pvalue","pfam_fg_freq","pfam_enrichment_string","GO_term","GO_pval","GO_fg_freq")
	gmo = gmo [ , c("gene","GO","GO_term","GO_ontology","GO_pval","GO_fg_freq","pfam_top_domain","pfam_top_pvalue","pfam_fg_freq","pfam_architecture","pfam_enrichment_string") ]


	# sort
	gmo$GO_pval = sprintf("%.1E", gmo$GO_pval)
	gmo$pfam_top_pvalue = sprintf("%.1E", gmo$pfam_top_pvalue)
	gmo = gmo [ order(as.numeric(gmo$GO_pval), gmo$GO, as.numeric(gmo$pfam_top_pvalue)), ]

	# ignore rows without gene (they really shouldn't be here but just in case)
	gmo_f = gmo
	gmo_f = gmo_f [ !is.na(gmo_f$gene), ]
	
	# if available, add extra annotations to summary
	if (!is.null(extra_annot_dict)) {
		if (is.list(extra_annot_dict)) {
			for (ndi in 1:length(extra_annot_dict)) {
				gmo_f[,sprintf("extra_annotation_%i", ndi)] = extra_annot_dict[[ndi]] [ gmo_f$gene ]
			}
		} else {
			gmo_f$extra_annotation_1 = extra_annot_dict [ gmo_f$gene ]	
		}
		
	}

	# save
	write.table(gmo_f, output_fn, row.names = FALSE, sep = "\t", quote = FALSE)
	
}

#### VENN DIAGRAMS ####
# plot venn diagrams from lists

venn.two = function(
		list1,
		list2,
		catname1,
		catname2,
		main="Venn",
		col1="green3",
		col2="magenta3",
		eulerbool=TRUE,
		print.mode=c("raw","percent")) {
	
	# compute set overlaps, intersections, etc
	list_uniq      = base::unique(c(list1, list2))
	list_intersect = base::intersect(list1, list2)
	list_diff_i    = base::setdiff(list1, list2)
	list_diff_j    = base::setdiff(list2, list1)
	
	# draw venn
	venn = draw.pairwise.venn(
		area1 = length(list1),
		area2 = length(list2),
		cross.area = length(list_intersect),
		category=c(catname1,catname2),
		fontfamily="Helvetica",
		cat.fontfamily = "Helvetica",
		col=c(col1,col2),ext.text = F,
		cat.col=c(col1,col2),
		ind = F, print.mode = print.mode,
		euler.d = eulerbool,
		scaled = eulerbool)
	grid::grid.newpage()
	grid::grid.draw(venn)
	grid::grid.text(main,x=0.5,y=0.95,
					gp = gpar(fontsize = 8, fontface = "bold"))
	
	# output
	output = list(
		"catname_i"      = catname1,
		"catname_j"      = catname2,
		"title"          = main,
		"list_uniq"      = list_uniq,
		"list_intersect" = list_intersect,
		"list_diff_i"    = list_diff_i,
		"list_diff_j"    = list_diff_j,
		"venn" = venn
	)
	
	return(output)
	
}



venn.three = function(
		list1,
		list2,
		list3,
		catname1,
		catname2,
		catname3,
		main="Venn",
		col1="green3",
		col2="magenta3",
		col3="orange",
		eulerbool=TRUE,
		print.mode=c("raw","percent")) {
	
	# compute set overlaps, intersections, etc
	intersect_n12  = base::intersect(list1, list2)
	intersect_n13  = base::intersect(list1, list3)
	intersect_n23  = base::intersect(list2, list3)
	intersect_n123 = base::intersect(intersect_n12, intersect_n13)
	
	# draw venn
	venn = draw.triple.venn(
		area1 = length(list1),
		area2 = length(list2),
		area3 = length(list3),
		n12 = length(intersect_n12),
		n13 = length(intersect_n13),
		n23 = length(intersect_n23),
		n123 = length(intersect_n123),
		category=c(catname1,catname2,catname3),
		fontfamily="Helvetica",
		cat.fontfamily = "Helvetica",
		col=c(col1,col2,col3),ext.text = F,
		cat.col=c(col1,col2,col3),
		ind = F, print.mode = print.mode,
		euler.d = eulerbool, scaled = T)
	grid::grid.newpage()
	grid::grid.draw(venn)
	grid::grid.text(main,x=0.5,y=0.95,
					gp = gpar(fontsize = 8, fontface = "bold"))
	
	# output
	output = list(
		"venn" = venn,
		"catname_1"      = catname1,
		"catname_2"      = catname2,
		"catname_3"      = catname3,
		"intersect_12" = intersect_n12,
		"intersect_13" = intersect_n13,
		"intersect_23" = intersect_n23,
		"intersect_123" = intersect_n123
	)
	
	return(output)
	
}
