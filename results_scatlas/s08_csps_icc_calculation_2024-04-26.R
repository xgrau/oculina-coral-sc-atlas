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
sps_refl = c("Ocupat")

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

# loop
for (spi in sps_refl) {

	# set working species
	if (spi == "Spio" | spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}

	# output
	out_fn = sprintf("results_metacell_%s_filt/csps/", spi)
	dir.create(out_fn, showWarnings = FALSE)

	## Load annotations ##

	# load gene annotations for this species
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
	
	# seurat gene dict
	gg_v = dic_from_vecs(gsub("_","-", gna_p$gene), gna_p$gene)
	
	
	for (spj in sps_list [ sps_list != spi ]) {
		
		# set working species
		if (spj == "Spio" | spj == "Spin") {
			spj_w = "Spis"
		} else {
			spj_w = spj
		}
		message(sprintf("csps | %s | compare to %s, load footprints...", spi, spj))
		inp_fn_i = sprintf("results_metacell_%s_filt/", spi)
		inp_fn_j = sprintf("results_metacell_%s_filt/", spj)
		mc_fp_i = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", inp_fn_i, spi))
		mc_fp_j = readRDS(sprintf("%s/dat.%s.expression.mcs_fp.rds", inp_fn_j, spj))

		message(sprintf("csps | %s | compare to %s, load counts...", spi, spj))
		mc_umi_i = readRDS(sprintf("%s/dat.%s.expression.mcs_umicount.rds", inp_fn_i, spi))
		mc_umi_j = readRDS(sprintf("%s/dat.%s.expression.mcs_umicount.rds", inp_fn_j, spj))

		message(sprintf("csps | %s | compare to %s, load ortholog pairs...", spi, spj))
		og_pairs_ij = clean_og_pairs(
			og_pairs_fn = "../data/orthology_Anthozoa_plus/orthogroup_pairs.genes.csv", 
			sp1 = spi_w, 
			sp2 = spj_w, 
			t2g = FALSE,
			header = TRUE
		)
		og_pairs_ij[,1] = gsub("_","-", og_pairs_ij[,1])
		og_pairs_ij[,2] = gsub("_","-", og_pairs_ij[,2])

		message(sprintf("csps | %s | compare to %s, drop genes with too few counts...", spi, spj))
		counts_per_gene_i = rowSums(mc_umi_i)
		counts_per_gene_j = rowSums(mc_umi_j)
		good_genes_i = names(counts_per_gene_i) [ counts_per_gene_i >= 20 ]
		good_genes_j = names(counts_per_gene_j) [ counts_per_gene_j >= 20 ]
		good_genes_i = good_genes_i [ good_genes_i %in% rownames(mc_fp_i) ]
		good_genes_j = good_genes_j [ good_genes_j %in% rownames(mc_fp_j) ]

		message(sprintf("csps | %s | compare to %s, run ICC...", spi, spj))
		icc = csps_markers_icc_noobj(
			mat_sp1 = mc_fp_i[good_genes_i,],
			mat_sp2 = mc_fp_j[good_genes_j,],
			og_pairs = og_pairs_ij,
			do_duplicates = TRUE,
			nthreads_icc = 32,
			nthreads_dup = 1
		)
		write.table(icc$ec_markers, sprintf("%s/dat.icc.%s-%s.ec_scores.csv", out_fn, spi, spj), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		write.table(icc$ec_duplicates, sprintf("%s/dat.icc.%s-%s.ec_scores.duplicates.csv", out_fn, spi, spj), row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
		saveRDS(icc, sprintf("%s/dat.icc.%s-%s.obj.rds", out_fn, spi, spj))
		
	}
	
}

message("All done!")
