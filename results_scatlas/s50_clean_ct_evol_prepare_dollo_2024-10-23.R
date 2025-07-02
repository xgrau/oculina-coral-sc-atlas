# libraries
suppressMessages(library("Seurat"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
graphics.off()
suppressMessages(library("ape"))

# output
out_fn = "results_ct_evolution/"
dir.create(out_fn, showWarnings = FALSE)

# list of comparisons to perform: for each cell type, select which clusters in each species are meant to be included in the reconstruction
comp_list = list(
	"gastrodermis_alga_hosting_symbio" = list(Ocupat = "gastrodermis_alga_hosting", Ocuarb = "gastrodermis_alga_hosting", Spin = "gastrodermis_alga_hosting", Amil = "gastrodermis_alga_hosting"),
	"gastrodermis" = list(Ocupat = "gastrodermis", Ocuarb = "gastrodermis", Spin = "gastrodermis", Amil = c("gastrodermis_1","gastrodermis_2"), Nvec = "gastrodermis", Xesp = "gastrodermis"),
	"glanddigestive" = list(Ocupat = "gland_3", Ocuarb = "gland_3", Spin = "gland_2", Amil = "gland_4", Nvec = "gland_17"),
	"glandmucin" = list(Ocupat = "gland_8", Ocuarb = "gland_8", Spin = "gland_3", Amil = "gland_2", Nvec = "gland_4"),
	"digestfil" = list(Ocupat = "digestive_filaments", Ocuarb = "digestive_filaments", Spin = "digestive_filaments", Amil = "digestive_filaments", Nvec = "digestive_filaments")
)

# thresholds
thr_fps = 1.1
thr_pv  = 1e-3

# load orthology
ogm = read.table("../data/orthology_Metazoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
ogm_gv = dic_from_vecs(names = ogm$gene, terms = ogm$orthogroup_name)
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_gov = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
ogm_tv = dic_from_vecs(names = ogm$transcript, terms = ogm$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)
gg_v = dic_from_vecs(gsub("_","-", oga$gene), oga$gene)


# based on markers fc
seu_l = list()
ufs_l = list()
umi_l = list()
# load seurat and umifrac
for (spi in c("Ocupat","Ocuarb","Spin","Amil","Nvec","Xesp")) {
		
	message(sprintf("load %s", spi))
	ufs_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umifrac.rds", spi, spi))
	umi_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umicount.rds", spi, spi))
	seu_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.seurat_final.rds", spi, spi))
	
}


# loop over marker sets
for (nn in 1:length(comp_list)) {

	cti = names(comp_list)[[nn]]
	ctv = comp_list[[nn]]
	message(sprintf("evo | %s | init...", cti))
	
	# tree specific for each case...
	sps_list = names(comp_list[[nn]])

	# get markers
	mks_l = list()
	mks_t = data.frame()
	for (spi in sps_list) {
		
		if (spi == "Spin") {
			spi_w = "Spis"
		} else {
			spi_w = spi
		}

		message(sprintf("evo | %s | calculate markers %s...", cti, spi))
		# umi and umifrac in that ct/group
		umi_i_v = rowSums(as.matrix(umi_l[[spi]][,ctv[[spi]]]))
		names(umi_i_v) = gg_v [ names(umi_i_v) ]
		ufs_i_v = umi_i_v / sum(umi_i_v) * 1e4
		
		# test
		# SeuratObject::Idents(seu_l[[spi]]) = seu_l[[spi]]@meta.data[,"cell_type"]
		# mks_i = Seurat::FindMarkers(seu_l[[spi]], ident.1 = ctv[[spi]], layer = "RNA", slot = "counts", fc.slot = "counts", min.pct = 0, test.use = "wilcox")
		# first select good cells
		if (grepl("symbio", cti) & spi %in% c("Ocupat","Ocuarb")) {
			cells_fg = rownames(seu_l[[spi]]@meta.data) [ seu_l[[spi]]@meta.data[,"cell_type"] %in% ctv[[spi]] & seu_l[[spi]]@meta.data$batch_method %in% c("Ocupat02_OPUB","Ocupat04_Opat02U","Ocuarb01_sym_SRR29367137") ]
			cells_bg = rownames(seu_l[[spi]]@meta.data) [ ! rownames(seu_l[[spi]]@meta.data) %in% cells_fg     & seu_l[[spi]]@meta.data$batch_method %in% c("Ocupat02_OPUB","Ocupat04_Opat02U","Ocuarb01_sym_SRR29367137") ]
		} else if (grepl("aposym", cti) & spi %in% c("Ocupat","Ocuarb")) {
			cells_fg = rownames(seu_l[[spi]]@meta.data) [ seu_l[[spi]]@meta.data[,"cell_type"] %in% ctv[[spi]] & !seu_l[[spi]]@meta.data$batch_method %in% c("Ocupat02_OPUB","Ocupat04_Opat02U","Ocuarb01_sym_SRR29367137") ]
			cells_bg = rownames(seu_l[[spi]]@meta.data) [ ! rownames(seu_l[[spi]]@meta.data) %in% cells_fg     & !seu_l[[spi]]@meta.data$batch_method %in% c("Ocupat02_OPUB","Ocupat04_Opat02U","Ocuarb01_sym_SRR29367137") ]
		} else {
			cells_fg = rownames(seu_l[[spi]]@meta.data) [ seu_l[[spi]]@meta.data[,"cell_type"] %in% ctv[[spi]] ]
			cells_bg = rownames(seu_l[[spi]]@meta.data) [ ! rownames(seu_l[[spi]]@meta.data) %in% cells_fg     ]
		}
		mks_i = Seurat::FindMarkers(seu_l[[spi]]@assays$RNA, slot = "counts", fc.slot = "counts", cells.1 = cells_fg, cells.2 = cells_bg, min.pct = 0, test.use = "wilcox")
		# add info to markers table
		mks_i = mks_i [ !grepl("orphan", rownames(mks_i)), ]
		mks_i$gene = gg_v [ rownames(mks_i) ]
		mks_i = mks_i [ !is.na(mks_i$gene) , ]
		mks_i$orthogroup = oga_gov [ mks_i$gene ]
		mks_i$gene_name = ogm_gv [ mks_i$gene ]
		mks_i$cell_type = ctv[[spi]][1]
		mks_i$umi = umi_i_v [ mks_i$gene ]
		mks_i$umifrac = ufs_i_v [ mks_i$gene ]
		mks_i$fc = as.numeric(2^mks_i$avg_log2FC)
		rownames(mks_i) = NULL
		# store awaky for this species
		write.table(mks_i, sprintf("%s/mks.%s.%s.all.tsv", out_fn, cti, spi), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		mks_l[[spi]] = mks_i
		
		# merge with other species
		mks_f = mks_l[[spi]] [ , c("gene","orthogroup","umi","umifrac","fc","p_val_adj","pct.1")]
		# add nonexpressed markers
		mks_f_nonexp = data.frame(gene = names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ], orthogroup = oga_gov [ names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ] ])
		mks_f_nonexp$umi = 0
		mks_f_nonexp$umifrac = 0
		mks_f_nonexp$fc = 0
		mks_f_nonexp$p_val_adj = 1
		mks_f_nonexp$pct.1 = 0
		mks_f = rbind(mks_f, mks_f_nonexp)
		# is significant?
		mks_f$significant = mks_f$p_val_adj < thr_pv & mks_f$fc >= thr_fps
		# add gene name
		mks_f$gene_name = ogm_gv [ mks_f$gene ]
		mks_f = mks_f [ !is.na(mks_f$gene_name), ]
		mks_f = mks_f [ order(mks_f$orthogroup, -as.numeric(mks_f$significant), mks_f$p_val_adj, -mks_f$fc) , ]
		# add pfam domain
		pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi_w))
		names(pf_annot) = oga_gtv[names(pf_annot)]
		pf_annot_v = lapply(pf_annot, function(vv) { paste(vv, collapse = "/") })
		mks_f$pfam_domains = as.character(pf_annot_v[mks_f$gene])
		mks_f$pfam_domains [ mks_f$pfam_domains == "NULL" ] = ""
		# add species
		mks_f$species = spi
		mks_t = rbind(mks_t, mks_f)
	
	}
	# merged table
	mks_t$species = factor(mks_t$species, levels = sps_list)

	# load set of markers for each species
	ogs_d = data.frame()
	gen_l = list()
	for (spi in sps_list) {
		message(sprintf("evo | %s | select markers %s...", cti, spi))
		# save og+gene info
		gen_i = mks_l[[spi]] [ mks_l[[spi]]$fc > thr_fps & mks_l[[spi]]$p_val_adj < thr_pv ,]
		gen_i = gen_i [ order(gen_i$orthogroup, -gen_i$fc), ]
		gen_l[[spi]] = gen_i
		write.table(gen_i, sprintf("%s/mks.%s.%s.filt.tsv", out_fn, cti, spi), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
		# store og info only
		ogs_d_i = data.frame(orthogroup = oga_gov [ gen_i$gene ], gene = gen_i$gene, species = spi)
		ogs_d_i = ogs_d_i [ !is.na(ogs_d_i$orthogroup), ]
		ogs_d = rbind(ogs_d, ogs_d_i)
	}

	# save
	write.table(ogs_d, sprintf("%s/anc.%s.mat.extant.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	# resort merged table
	mks_t_o = unique(mks_t [ , c("orthogroup","significant","species") ])
	mks_t_median_uf_per_og = aggregate(umifrac ~ orthogroup, mks_t, mean)
	mks_t_numsignif_v = table(mks_t_o$orthogroup, mks_t_o$significant)[,"TRUE"]
	mks_t_median_uf_per_og$numsignif = mks_t_numsignif_v [ mks_t_median_uf_per_og$orthogroup ]
	mks_t_median_uf_per_og = mks_t_median_uf_per_og [ order(-mks_t_median_uf_per_og$numsignif, -mks_t_median_uf_per_og$umifrac), ]
	mks_t$orthogroup = factor(mks_t$orthogroup, levels = mks_t_median_uf_per_og$orthogroup)
	# filter ogs that are in markers list, and also best gene per og
	mks_t = mks_t [ order(mks_t$orthogroup, mks_t$species), ]
	mks_t_f = mks_t [ mks_t$orthogroup %in% ogs_d$orthogroup, ] 
	mks_t_f_best = mks_t_f [ !duplicated(paste(mks_t_f$orthogroup, mks_t_f$species)), ]

	# save
	write.table(mks_t_f, sprintf("%s/mks.%s.allsps.all_genes_per_og.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	write.table(mks_t_f_best, sprintf("%s/mks.%s.allsps.top_genes_per_og.tsv", out_fn, cti), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
	
	
	# loop over pairs of species
	for (spi in sps_list) {
		
		pdf(sprintf("%s/mks.%s.%s.all.pdf", out_fn, cti, spi), width = 20, height = 4)
		layout(matrix(1:5, nrow = 1))
		
		# volcano plot sps i
		mks_i = mks_l[[spi]]
		gen_i = gen_l[[spi]]
		mks_i = mks_i [ order(-mks_i$p_val_adj), ]
		# cap padj
		mks_i_logpadj = -log10(mks_i$p_val_adj)
		mks_i_logpadj [ mks_i_logpadj > 12 ] = 12
		# hist(mks_i_logpadj)
		mks_i_logpadj_v = cut(mks_i_logpadj, c(0, 3, 6, 9, 12), include.lowest = TRUE)
		color_palette = colorRampPalette(interpolate="l",c("gray90", "#9ae445","#44aa2c","#1a622d"))
		mks_i_logpadj_c = color_palette(nlevels(mks_i_logpadj_v))  [ mks_i_logpadj_v ]
		# cap fc
		mks_i_logfc = mks_i$avg_log2FC
		mks_i_logfc [ mks_i_logfc > 12 ] = 12
		mks_i_logfc [ mks_i_logfc < -12 ] = -12
		# cap umi
		mks_i_umi = mks_i$umi
		mks_i_umi [ mks_i_umi > 1e4 ] = 1e4
		plot(mks_i_logfc, mks_i_umi, col = mks_i_logpadj_c, pch = 19, ylim = c(1,1e4), xlim = c(-12,12), log = "y", las = 1)
		title(main = sprintf("%s\nFC v UMI in %s", spi, cti), cex.main = 1, font.main = 1)
		abline(v = c(-log2(1.1),0,log2(1.1)), lty = 2, col = "darkred")
		
		# plot fc v pvalue
		plot(mks_i$avg_log2FC, mks_i_logpadj, col = mks_i_logpadj_c, pch = 19, las = 1, xlim = c(-12,12), ylim = c(0,12))
		title(main = sprintf("%s\nFC v UMI in %s", spi, cti), cex.main = 1, font.main = 1)
		abline(v = c(-log2(1.1),0,log2(1.1)), lty = 2, col = "darkred")
		
		layout(matrix(1:5, nrow = 1))
		for (spj in sps_list) {
			if (spi != spj) {
				message(sprintf("evo | %s | compare %s-%s...", cti, spi,spj))
				# load
				gen_i = gen_l[[spi]]
				gen_j = gen_l[[spj]]
				mks_i = mks_l[[spi]]
				mks_j = mks_l[[spj]]
				# ogs_i_only = 
				# ogs_j_only = 
				# sort
				gen_i = gen_i [ !duplicated(gen_i$orthogroup), ]
				mks_i = mks_i [ order(mks_i$orthogroup, mks_i$p_val_adj, -mks_i$fc), ]
				mks_i = mks_i [ !duplicated(mks_i$orthogroup), ]
				# mks_i = mks_i [ mks_i$orthogroup %in% c(gen_i$orthogroup, gen_j$orthogroup),  ]
				# same for j
				gen_j = gen_j [ !duplicated(gen_j$orthogroup), ]
				mks_j = mks_j [ order(mks_j$orthogroup, mks_j$p_val_adj, -mks_j$fc), ]
				mks_j = mks_j [ !duplicated(mks_j$orthogroup), ]
				# mks_j = mks_j [ mks_j$orthogroup %in% c(gen_i$orthogroup, gen_j$orthogroup),  ]
				# merge
				# gen_m = merge(gen_i[,c("orthogroup","avg_log2FC","p_val_adj","gene")],  gen_j[,c("orthogroup","avg_log2FC","p_val_adj","gene")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c("_i","_j"))
				mks_m = merge(mks_i[,c("orthogroup","avg_log2FC","p_val_adj","gene")],  mks_j[,c("orthogroup","avg_log2FC","p_val_adj","gene")], by = "orthogroup", all.x = TRUE, all.y = TRUE, suffixes = c("_i","_j"))
				# OG classification
				mks_m$og_class = factor(NA, levels = c("i-only","j-only","shared"))
				mks_m$og_class [ mks_m$orthogroup %in% setdiff(gen_i$orthogroup, gen_j$orthogroup) ] = "i-only"
				mks_m$og_class [ mks_m$orthogroup %in% setdiff(gen_j$orthogroup, gen_i$orthogroup) ] = "j-only"
				mks_m$og_class [ mks_m$orthogroup %in% intersect(gen_i$orthogroup, gen_j$orthogroup) ] = "shared"
				mks_m_og_class_c = c("mediumpurple1","plum2","blue4")
				mks_m = mks_m [ order(mks_m$og_class), ]
				mks_m = mks_m [ mks_m$gene_i %in% gen_i$gene | mks_m$gene_j %in% gen_j$gene, ]
				names(mks_m_og_class_c) = levels(mks_m$og_class)				
				plot(mks_m$avg_log2FC_i, mks_m$avg_log2FC_j, col = mks_m_og_class_c[mks_m$og_class], pch = 19, xlab = spi, ylab = spj, cex = 0.8)
				abline(v = c(0,log2(1.1)), h = c(0,log2(1.1)), lty = 2, col = "darkred")
				title(main = sprintf("%s\nFC %s-%s", cti, spi, spj), cex.main = 1, font.main = 1)
			} else {
				plot(0)
				legend("topleft", c("sps i-only","sps j-only","shared"), col = c("mediumpurple1","plum2","blue4"), pch = 19, bty = "n")
			}
		}
		dev.off()
	}
	
}

message("All done!")
