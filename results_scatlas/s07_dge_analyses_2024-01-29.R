# libraries
suppressMessages(source("../scripts/Seurat_functions.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))
graphics.off()

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)

sps_list = c("Ocupat","Ocuarb")

# list of comparisons to perform for each species (g1 and g2 = groups of cells 1 and 2, class = metadata field used to define groups of cells,
# subset = groups of cells to subset in this analysis (e.g. a cell type), subset_class = metadata field to define subset, i.e. "cell_type")
list_dge = list(
	"Ocupat" = list(
		"ct_hostcells.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gastrodermis_alga_hosting"), subset_class = c("cell_type")),
		"ct_gastrodermis.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gastrodermis"), subset_class = c("cell_type")),
		"ct_immune.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("immune_macrophage_like","immune_like_1","immune_like_2"), subset_class = c("cell_type")),
		"ct_immune_macrophage_like.sam_ubl_v_bl" = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("immune_macrophage_like"), subset_class = c("cell_type")),
		"ct_immune_like1.sam_ubl_v_bl" = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("immune_like_1"), subset_class = c("cell_type")),
		"ct_immune_like2.sam_ubl_v_bl" = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("immune_like_2"), subset_class = c("cell_type")),
		"ct_neuron.sam_ubl_v_bl"      = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Pou4_1","neuron_Pou4_2","neuron_Isl","neuron_Pou4_Gsx","neuron_Pou4_Otp_1","neuron_Pou4_Otp_2"), subset_class = c("cell_type")),
		"ct_neuron_Isl.sam_ubl_v_bl"  = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Isl"), subset_class = c("cell_type")),
		"ct_neuron_Pou.sam_ubl_v_bl"  = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Pou4_1","neuron_Pou4_2","neuron_Pou4_Gsx","neuron_Pou4_Otp_1","neuron_Pou4_Otp_2"), subset_class = c("cell_type")),
		"ct_neuron_Pou12.sam_ubl_v_bl"  = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Pou4_1","neuron_Pou4_2"), subset_class = c("cell_type")),
		"ct_neuron_PouGsx.sam_ubl_v_bl"  = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Pou4_Gsx"), subset_class = c("cell_type")),
		"ct_neuron_PouOtp.sam_ubl_v_bl"  = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neuron_Pou4_Otp_1","neuron_Pou4_Otp_2"), subset_class = c("cell_type")),
		"ct_gland1.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_1_Xbp"), subset_class = c("cell_type")),
		"ct_gland2.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_2"), subset_class = c("cell_type")),
		"ct_gland3.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_3"), subset_class = c("cell_type")),
		"ct_gland4.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_4"), subset_class = c("cell_type")),
		"ct_gland5.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_5"), subset_class = c("cell_type")),
		"ct_gland6.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_6"), subset_class = c("cell_type")),
		"ct_gland7.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_7"), subset_class = c("cell_type")),
		"ct_gland8.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_8"), subset_class = c("cell_type")),
		"ct_gland.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gland_1_Xbp","gland_2","gland_3","gland_4","gland_5","gland_6","gland_7","gland_8","gland_9","gland_10","gland_11"), subset_class = c("cell_type")),
		"ct_epidermis.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("epidermis","epidermis_like_1","epidermis_like_2"), subset_class = c("cell_type")),
		"ct_digestive_filaments.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("digestive_filaments"), subset_class = c("cell_type")),
		"ct_calicoblast.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("calicoblast"), subset_class = c("cell_type")),
		"ct_neurosecretory_progenitors.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("neurosecretory_progenitors"), subset_class = c("cell_type")),
		"ct_cnidocyte.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("cnidocyte_1","cnidocyte_2"), subset_class = c("cell_type")),
		"ct_germline_oocytes.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("germline_oocytes"), subset_class = c("cell_type")),
		"ct_gastrodermis_like_1.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gastrodermis_like_1"), subset_class = c("cell_type")),
		"ct_gastrodermis_like_2.sam_ubl_v_bl"   = list(g1 = c("Ocupat02_OPUB","Ocupat04_Opat02U"), g2 = c("Ocupat01_OPB","Ocupat03_Opat02B"), class = "batch_method", subset = c("gastrodermis_like_2"), subset_class = c("cell_type"))
	),
	"Ocuarb" = list(
		"ct_hostcells_v_gast.sam_all" = list(g1 = c("gastrodermis_alga_hosting"), g2 = c("gastrodermis"), class = "cell_type", subset = NULL, subset_class = NULL),
		"ct_hostcells_v_gast.sam_sym" = list(g1 = c("gastrodermis_alga_hosting"), g2 = c("gastrodermis"), class = "cell_type", subset = "Ocuarb01_sym_SRR29367137", subset_class = "batch_method"),
		"ct_hostcells_v_gast.sam_apo" = list(g1 = c("gastrodermis_alga_hosting"), g2 = c("gastrodermis"), class = "cell_type", subset = "Ocuarb02_apo_SRR29367138", subset_class = "batch_method"),
		"ct_hostcells_v_gast.sam_all" = list(g1 = c("gastrodermis_alga_hosting"), g2 = c("gastrodermis_1","gastrodermis_2"), class = "cell_type", subset = NULL, subset_class = NULL),
		"ct_hostcells.sam_ubl_v_bl" = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = "gastrodermis_alga_hosting", subset_class = "cell_type"),
		"ct_gastrodermis.sam_ubl_v_bl" = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = "gastrodermis", subset_class = "cell_type"),
		"ct_gastrodermis.sam_ubl_v_bl" = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gastrodermis_1","gastrodermis_2"), subset_class = "cell_type"),
		"ct_immune.sam_ubl_v_bl" = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = "immune", subset_class = "cell_type"),
		"ct_gland1.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_1_Xbp"), subset_class = c("cell_type")),
		"ct_gland2.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_2"), subset_class = c("cell_type")),
		"ct_gland3.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_3"), subset_class = c("cell_type")),
		"ct_gland4.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_4"), subset_class = c("cell_type")),
		"ct_gland5.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_5"), subset_class = c("cell_type")),
		"ct_gland6.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_6"), subset_class = c("cell_type")),
		"ct_gland8.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_8"), subset_class = c("cell_type")),
		"ct_epidermis.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("epidermis"), subset_class = c("cell_type")),
		"ct_gland.sam_ubl_v_bl"   = list(g1 = c("Ocuarb01_sym_SRR29367137"), g2 = c("Ocuarb02_apo_SRR29367138"), class = "batch_method", subset = c("gland_1_Xbp","gland_2","gland_3","gland_4","gland_5","gland_6","gland_8"), subset_class = c("cell_type"))
	)
)

# significance thresholds
sig_thr_p = 1e-6 # or 1e-3 for nonstrict
sig_thr_log2fc_for_plot = 1
sig_thr_fc_for_dge = 1.1


for (spi in sps_list) {

	## Create dir ##
	
	# set working species
	if (spi == "Spio" | spi == "Spin") {
		spi_w = "Spis"
	} else {
		spi_w = spi
	}

	# output
	out_fn = sprintf("results_metacell_%s_filt/", spi)
	dir.create(out_fn, showWarnings = FALSE)
	dir.create(sprintf("%s/tests_dge_nov24/", out_fn), showWarnings = FALSE)

	## Load annotations ##

	# load gene annotations for this species
	gna_i = ogm [ ogm$species == spi, ]
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

	# functional enrichments
	go_annot = gsa_topgo_load_emapper(emapper_fn = sprintf("../data/reference/%s_ensembl.GO.csv", spi_w), index_col_GOs = 2)
	pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi_w))
	names(go_annot) = oga_gtv [ names(go_annot)  ]
	names(pf_annot) = oga_gtv [ names(pf_annot)  ]

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
	
	# gene names from spis, for reference
	if (spi_w != "Spis") {
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


	## Load Seurat ##

	# seurat object with batch info
	message(sprintf("dge | %s | load Seurat...", spi))
	seu = readRDS(sprintf("%s/dat.%s.seurat_final.rds", out_fn, spi))
	seu$batch_method = factor(seu$batch_method)
	seu_lou_col_v = dic_from_vecs(seu$leiden, seu$leiden_color)
	seu_lou_col_v = seu_lou_col_v[levels(seu$leiden)]
	
	message(sprintf("dge | %s | load annotated clusters...", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", out_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
	ctt_lou_cts_v = dic_from_vecs(ctt$cluster, ctt$cell_type)
	ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	

	## Loop through dge list ##

	for (dei in names(list_dge[[spi]])) {
		
		# log
		d_g1 = list_dge[[spi]][[dei]][["g1"]]
		d_g2 = list_dge[[spi]][[dei]][["g2"]]
		d_cl = list_dge[[spi]][[dei]][["class"]]
		d_su = list_dge[[spi]][[dei]][["subset"]]
		d_cs = list_dge[[spi]][[dei]][["subset_class"]]
		
		# get cells list
		if (!is.null(d_cs)) {
			message(sprintf("dge | %s | markers %s, %s ~ %s (%s)...", spi, dei, paste(d_g1, collapse = ","), paste(d_g2, collapse = ","), paste(d_su, collapse = ",")))
			cells1 = rownames(seu@meta.data) [ seu@meta.data[,d_cl] %in% d_g1 & seu@meta.data[,d_cs] %in% d_su ]
			cells2 = rownames(seu@meta.data) [ seu@meta.data[,d_cl] %in% d_g2 & seu@meta.data[,d_cs] %in% d_su ]
		} else {
			message(sprintf("dge | %s | markers %s, %s ~ %s...", spi, dei, paste(d_g1, collapse = ","), paste(d_g2, collapse = ",")))
			cells1 = rownames(seu@meta.data) [ seu@meta.data[,d_cl] %in% d_g1 ]
			cells2 = rownames(seu@meta.data) [ seu@meta.data[,d_cl] %in% d_g2 ]
		}
		
		# dge per se
		mki = Seurat::FindMarkers(seu@assays$RNA, slot = "counts", fc.slot = "counts", cells.1 = cells1, cells.2 = cells2, min.pct = 0.01, test.use = "wilcox")

		# counts in cells of interest
		counts_per_gene = rowSums(seu@assays$RNA$counts[,c(cells1, cells2)])
		names(counts_per_gene) = gg_v [ names(counts_per_gene) ]
		counts_per_gene = counts_per_gene [ !is.na(names(counts_per_gene)) ]
		expressed_genes_in_focus = names(counts_per_gene) [ counts_per_gene > 0 ]
		
		# save dge info
		mki$focus = dei
		mki$group1 = paste(d_g1, collapse = ",")
		mki$group2 = paste(d_g2, collapse = ",")
		mki$class = d_cl
		mki$subset = paste(d_su, collapse = ",")
		mki$subset_class = d_cs
	
		# dictionary to reconstruct gene names
		mki$gene = rownames(mki)
		rownames(mki) = NULL
		mki$gene = gg_v [ mki$gene ]
		mki$gene_name = ogm_gv [mki$gene]
		mki$pfam = pfa_v [ mki$gene ]
		mki$gene_name_Spis = gnspis_v [ mki$gene ]
		mki$is_tf = mki$gene %in% gene_annot_tfs$gene
		mki$p_val_adj = as.numeric(sprintf("%.2E", mki$p_val_adj))
		mki$p_val = as.numeric(sprintf("%.2E", mki$p_val))
		mki$avg_log2FC = as.numeric(sprintf("%.2E", mki$avg_log2FC))
		mki$counts = counts_per_gene [ mki$gene ]
		
		# sort
		mki = mki [ order(-mki$avg_log2FC, -mki$p_val_adj), ]

		# save
		mki = mki [ !is.na(mki$gene), ]
		write.table(mki, sprintf("%s/tests_dge_nov24/dge.%s.%s.csv", out_fn, spi, dei), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
		xlsx::write.xlsx2(mki, sprintf("%s/tests_dge_nov24/dge.%s.%s.xlsx", out_fn, spi, dei), row.names = FALSE)

		# prepare plot
		# mki$pct_log2FC = log2((mki$pct.1 + 1e-4) / (mki$pct.2 + 1e-4))
		mki$log_p = -log10(mki$p_val_adj)
		pct_num_col = 4
		pct_log_lim = ceiling(max(abs(mki$log_p)))
		color_palette = colorRampPalette(interpolate="l",c("gray90", "#9ae445","#44aa2c","#1a622d"))
		mki$log_p_cat = cut(mki$log_p, breaks = c(-1,-log10(sig_thr_p), -log10(sig_thr_p) * 2,  -log10(sig_thr_p) * 3, Inf))
		pct_col = color_palette(pct_num_col)
		names(pct_col) = levels(mki$log_p_cat)
		bool_sig  = mki$log_p >= -log10(sig_thr_p)
		bool_sig_1 = bool_sig & mki$avg_log2FC >= log2(1.1) 
		bool_sig_2 = bool_sig & mki$avg_log2FC <= log2(1.1) 
		
		# plot
		pdf(sprintf("%s/tests_dge_nov24/dge.%s.%s.volcano.pdf", out_fn, spi, dei), height = 8, width = 8)
		plot(mki$avg_log2FC[!bool_sig], mki$counts[!bool_sig], pch = 19, cex=0.7, col = pct_col [ mki$log_p_cat[!bool_sig] ], ylab = "counts", xlab = "", log = "y", las = 1, xlim = c(-max(abs(mki$avg_log2FC)), max(abs(mki$avg_log2FC))), ylim = c(1,max(counts_per_gene)))
		points(mki$avg_log2FC[bool_sig], mki$counts[bool_sig], pch = 19, cex=0.7, col = pct_col [ mki$log_p_cat[ bool_sig] ])
		abline(v=0, lty=2)
		ixx = mki$p_val_adj <= sig_thr_p & abs(mki$avg_log2FC) >= sig_thr_log2fc_for_plot & grepl(":", mki$gene_name)
		ixg = mki$gene_name[ixx]
		ixg = stringr::str_trunc(ixg, 50)
		if (length(ixg) > 0) {
			text(mki$avg_log2FC[ixx], mki$counts[ixx], ixg, col = scales::alpha("darkblue", 0.3), cex = 0.6, pos = ifelse(mki$avg_log2FC[ixx]>0, 4, 2))
		}
		title(sub = sprintf("n=%i (+%i / -%i) DGE genes (p<%.1E)\n-%s / +%s", sum(bool_sig), sum(bool_sig_1), sum(bool_sig_2), sig_thr_p, paste(d_g2, collapse = ","), paste(d_g1, collapse = ",")), main = sprintf("%s\nlog2(fc)", dei))
		legend("bottomleft", names(pct_col), col = pct_col, cex = 0.7, pch = 19, title = "-log10p")
		dev.off()
		
		# functional enrichments, first group
		mki_genes = mki [ mki$p_val_adj <= sig_thr_p & mki$avg_log2FC >= log2(sig_thr_fc_for_dge) , "gene" ]
		mki_genes = mki_genes [ !is.na(mki_genes) ]
		# bgi_genes = gg_v [rownames(seu)]
		# bgi_genes = bgi_genes [ !is.na(bgi_genes) ]
		bgi_genes = expressed_genes_in_focus
		
		message(sprintf("functional enrichments | %s G1, n=%i genes", dei, length(mki_genes)))
		if (length(mki_genes)>0) {
			
			# # KEGG enrichment
			# kegg_enrichment = gsa_enrichment_hypergeometric_v2(
			# 	annotation = keg_lll,
			# 	genes_fg = mki_genes,
			# 	genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(keg_lll) ],
			# 	out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G1.kegg.csv", spi, spi, make.names(dei)),
			# 	name_fg = "G1"
			# )
			# keg_ld_i = keg_ld [ keg_ld$gene %in% mki_genes, ]
			# keg_ld_i$kegg_name = factor(keg_ld_i$kegg_name, levels = kegg_enrichment[[1]]$annot)
			# keg_ld_i = keg_ld_i [ order(keg_ld_i$kegg_name, keg_ld_i$ko, keg_ld_i$gene), ]
			
			# # KEGG enrichment pvalue
			# kegg_enrichment_pv_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$pval_adj)			
			# kegg_enrichment_fg_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$freq_in_fg)	
			# keg_ld_i$kegg_pvalue = kegg_enrichment_pv_v[keg_ld_i$kegg_name]		
			# keg_ld_i$kegg_n_in_fg = kegg_enrichment_fg_v[keg_ld_i$kegg_name]		
			# write.table(keg_ld_i, sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G1.keggsum.csv", spi, spi, make.names(dei)), sep = "\t", quote = FALSE, row.names = FALSE)
			
			# test domain enrichment
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = mki_genes,
				genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(pf_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G1.pfam.csv", spi, spi, make.names(dei)),
				name_fg = "G1"
			)

			# test GOs
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = mki_genes,
				genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(go_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G1.topgo.csv", spi, spi, make.names(dei)),
				name_fg = "G1",
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)
			
			# integrate GOs and pfam in a single table
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G1.summary.csv", spi, spi, make.names(dei)),
				extra_annot_dict = ogm_gv
			) 

		}

		# functional enrichments, second group
		mki_genes = mki [ mki$p_val_adj <= sig_thr_p & mki$avg_log2FC <= -log2(sig_thr_fc_for_dge) , "gene" ]
		mki_genes = mki_genes [ !is.na(mki_genes) ]
		# bgi_genes = gg_v [rownames(seu)]
		# bgi_genes = bgi_genes [ !is.na(bgi_genes) ]
		bgi_genes = expressed_genes_in_focus
		
		message(sprintf("functional enrichments | %s G2, n=%i genes", dei, length(mki_genes)))
		if (length(mki_genes)>0) {
			
			# # KEGG enrichment
			# kegg_enrichment = gsa_enrichment_hypergeometric_v2(
			# 	annotation = keg_lll,
			# 	genes_fg = mki_genes,
			# 	genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(keg_lll) ],
			# 	out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G2.kegg.csv", spi, spi, make.names(dei)),
			# 	name_fg = "G2"
			# )
			# keg_ld_i = keg_ld [ keg_ld$gene %in% mki_genes, ]
			# keg_ld_i$kegg_name = factor(keg_ld_i$kegg_name, levels = kegg_enrichment[[1]]$annot)
			# keg_ld_i = keg_ld_i [ order(keg_ld_i$kegg_name, keg_ld_i$ko, keg_ld_i$gene), ]
			
			# # KEGG enrichment pvalue
			# kegg_enrichment_pv_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$pval_adj)			
			# kegg_enrichment_fg_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$freq_in_fg)	
			# keg_ld_i$kegg_pvalue = kegg_enrichment_pv_v[keg_ld_i$kegg_name]		
			# keg_ld_i$kegg_n_in_fg = kegg_enrichment_fg_v[keg_ld_i$kegg_name]		
			# write.table(keg_ld_i, sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G2.keggsum.csv", spi, spi, make.names(dei)), sep = "\t", quote = FALSE, row.names = FALSE)
			

			# test domain enrichment
			pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
				annotation = pf_annot,
				genes_fg = mki_genes,
				genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(pf_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G2.pfam.csv", spi, spi, make.names(dei)),
				name_fg = "G2"
			)

			# test GOs
			go_enrichment_l = gsa_topgo_enrichment_v2(
				annotation = go_annot,
				genes_fg = mki_genes,
				genes_bg = gg_v[rownames(seu)] [ gg_v[rownames(seu)] %in% names(go_annot) ],
				out_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G2.topgo.csv", spi, spi, make.names(dei)),
				name_fg = "G2",
				ontologyset = c("BP","MF","CC"),
				tg_test = "fisher",
				tg_algorithm = "elim",
				top_markers = 30,
				nodesize = 10
			)
			
			# integrate GOs and pfam in a single table
			gsa_merge_tables_v2(
				main_enrichment = go_enrichment_l[[1]],
				main_mapping = go_enrichment_l[[2]],
				pfam_enrichment = pf_enrichment_l[[1]],
				pfam_mapping = pf_enrichment_l[[2]],
				output_fn = sprintf("results_metacell_%s_filt/tests_dge_nov24/dge.%s.%s.G2.summary.csv", spi, spi, make.names(dei)),
				extra_annot_dict = ogm_gv
			) 

		}
		
	}
	
}

message("All done!")
