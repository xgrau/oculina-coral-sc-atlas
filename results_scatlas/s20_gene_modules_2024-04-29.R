# libraries
suppressMessages(library("scales"))
suppressMessages(source("../scripts/Gene_module_functions.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/helper.R"))

# data index
sps_list = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp")

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)

# load gene names
gna_fn = "../data/orthology_Metazoa_plus/orthogroup_conservation.csv"
gna = read.table(gna_fn, sep = "\t", header = TRUE)

# dictionaries transcript to gene
gna_gtv = dic_from_vecs(names = gna$transcript, terms = gna$gene)
gna_tgv = dic_from_vecs(names = gna$gene, terms = gna$transcript)
oga_v = dic_from_vecs(names = gna$gene, terms = gna$orthogroup_name)

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
	
	message(sprintf("wgcna | %s | load annotations...", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)

	# define variable genes
	message(sprintf("wgcna | %s | define variable genes...", spi))
	fc_thr = 1.25
	var_genes = names(which(apply(mc_fp, 1, max) > fc_thr))
	var_genes = var_genes [ !grepl("orphan", var_genes) ]
	write.table(var_genes, sprintf("%s/gmod.%s.wgcna.var_genes.txt", out_fn, spi), quote = FALSE, row.names = FALSE, col.names = FALSE)
	
	# determine soft power
	message(sprintf("wgcna | %s | soft power...", spi))
	gmod_determineSoftPowerWGCNA(
		data = mc_fp[var_genes,], 
		output_file = sprintf("%s/gmod.%s.wgcna.softpower.pdf", out_fn, spi)
	)
	
	# run wgcna with soft power value determined above (first value above line)
	message(sprintf("wgcna | %s | wgcna...", spi))
	# mc_wgcna = readRDS(sprintf("%s/gmod.%s.wgcna.wgcna.rds", out_fn, spi))
	mc_wgcna = gmod_runWGCNA(
		data = mc_fp[var_genes,],
		propGenes = 1,
		softPower = 5,
		cor_method = "pearson",
		signedNetwork = TRUE)
	
	# plot dendrogram
	message(sprintf("wgcna | %s | dendrogram...", spi))
	gmod_plotModulesCut(
		mc_wgcna, 
		output_file = sprintf("%s/gmod.%s.wgcna.dendrogram.pdf", out_fn, spi)
	)
	
	# calculate module eigengenes
	message(sprintf("wgcna | %s | eigenvalues...", spi))
	mc_wgcna_me = gmod_calculate_module_eigengenes_v2(
		mc_wgcna,
		split = 2, 
		minClusterSize = 10, 
		cutHeight = 0.99)
	message(sprintf("wgcna | %s | eigenvalues, n=%i modules across %i clusters...", spi, ncol(mc_wgcna_me), nrow(mc_wgcna_me)))
	
	# overlapping memberships
	message(sprintf("wgcna | %s | gene memberships...", spi))
	mc_wgcna_gmods = gmod_moduleOverlapingMembership(
		mc_wgcna, 
		mc_wgcna_me, 
		kME_threshold = 0.5)
	mc_wgcna_gmods = mc_wgcna_gmods [ lengths(mc_wgcna_gmods) > 0 ]
	
	# make into a table...
	mc_wgcna_gmods_d = data.frame()
	for (nn in 1:length(mc_wgcna_gmods)) {
		mc_wgcna_gmods_d = rbind(mc_wgcna_gmods_d, data.frame(gene = mc_wgcna_gmods[[nn]], module = names(mc_wgcna_gmods)[nn]))
	}
	mc_wgcna_gmods_d$gene_clean = gg_v [ mc_wgcna_gmods_d$gene ]
	mc_wgcna_gmods_d$gene_clean [ is.na(mc_wgcna_gmods_d$gene_clean) ] = gsub("-","_", mc_wgcna_gmods_d$gene [ is.na(mc_wgcna_gmods_d$gene_clean) ])
	mc_wgcna_gmods_d$gene = mc_wgcna_gmods_d$gene_clean
	mc_wgcna_gmods_d$gene_clean = NULL
	mc_wgcna_gmods_d = merge(mc_wgcna_gmods_d, gene_annot, by = "gene", all.x = TRUE, all.y = FALSE)
	mc_wgcna_gmods_d$is_tf [ is.na(mc_wgcna_gmods_d$is_tf) ] = FALSE
	mc_wgcna_gmods_d$is_tf = factor(mc_wgcna_gmods_d$is_tf, levels = c(TRUE,FALSE))
	mc_wgcna_gmods_d = mc_wgcna_gmods_d [ order(mc_wgcna_gmods_d$module, mc_wgcna_gmods_d$is_tf, mc_wgcna_gmods_d$name) ,]
	write.table(mc_wgcna_gmods_d, sprintf("%s/gmod.%s.wgcna.memberships.csv", out_fn, spi), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
	
	# save objects for later use
	saveRDS(mc_wgcna,       sprintf("%s/gmod.%s.wgcna.wgcna.rds", out_fn, spi))
	saveRDS(mc_wgcna_gmods, sprintf("%s/gmod.%s.wgcna.memberships.rds", out_fn, spi))
	saveRDS(mc_wgcna_me,    sprintf("%s/gmod.%s.wgcna.ME.rds", out_fn, spi))
	
}

