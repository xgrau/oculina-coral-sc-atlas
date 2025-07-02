# libraries
suppressMessages(library("Seurat"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
graphics.off()
suppressMessages(library("ape"))

# output
out_fn = "results_ct_evolution_v2/"
dir.create(out_fn, showWarnings = FALSE)

# data index
sps_tree_s = c("(((((Ocupat,Ocuarb)Oculinidae,Spis)Robusta,Amil)Scleractinia,Nvec)Hexacorallia,Xesp)Anthozoa;")
sps_tree = ape::read.tree(text = sps_tree_s)
# ape::write.tree(sps_tree, sprintf("%s/small_tree.newick", out_fn))
# sps_tree$edge.length = rep(1e-5, nrow(sps_tree$edge)) # add tiny edge lengths, to replicate dollo
# sps_tree_sansnodes = sps_tree
# sps_tree_sansnodes$node.label = NULL
descendants_from_nodes = adephylo::listTips(sps_tree)
nod_list = c(sps_tree$tip.label, sps_tree$node.label)
sps_list = sps_tree$tip.label

# list of comparisons to perform: for each cell type, select which clusters in each species are meant to be included in the reconstruction
comp_list = list(
	"gastrodermis_alga_hosting_symbio" = list(Ocupat = "gastrodermis_alga_hosting", Ocuarb = "gastrodermis_alga_hosting", Spin = "gastrodermis_alga_hosting", Amil = "gastrodermis_alga_hosting"),
	"gastrodermis" = list(Ocupat = "gastrodermis", Ocuarb = "gastrodermis", Spin = "gastrodermis", Amil = c("gastrodermis_1","gastrodermis_2"), Nvec = "gastrodermis", Xesp = "gastrodermis"),
	"glanddigestive" = list(Ocupat = "gland_3", Ocuarb = "gland_3", Spin = "gland_2", Amil = "gland_4", Nvec = "gland_17"),
	"glandmucin" = list(Ocupat = "gland_8", Ocuarb = "gland_8", Spin = "gland_3", Amil = "gland_2", Nvec = "gland_4"),
	"digestfil" = list(Ocupat = "digestive_filaments", Ocuarb = "digestive_filaments", Spin = "digestive_filaments", Amil = "digestive_filaments", Nvec = "digestive_filaments")
)

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

# load evolutionary histories of all ogs
pr_ge_m = as.matrix(read.table(sprintf("%s/small_tree_plus.orthogroups.possvm.dollo_pres.csv", out_fn)))
ga_ge_m = as.matrix(read.table(sprintf("%s/small_tree_plus.orthogroups.possvm.dollo_gain.csv", out_fn)))
lo_ge_m = as.matrix(read.table(sprintf("%s/small_tree_plus.orthogroups.possvm.dollo_loss.csv", out_fn)))

# binarise
pr_ge_m = (pr_ge_m > 0) * 1
ga_ge_m = (ga_ge_m > 0) * 1
lo_ge_m = (lo_ge_m > 0) * 1

# loop over marker sets
for (nn in 1:length(comp_list)) {

	cti = names(comp_list)[[nn]]
	ctv = comp_list[[nn]]
	message(sprintf("evo | %s | init...", cti))
	
	# species included in comparison
	sps_list_in_comp = names(comp_list[[nn]])
	
	# load set of markers for each species
	pr_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_pres.csv", out_fn, cti)))
	ga_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_gain.csv", out_fn, cti)))
	lo_ex_m = as.matrix(read.table(sprintf("%s/anc.%s.mat.possvm.dollo_loss.csv", out_fn, cti)))
	
	# output
	dir.create(sprintf("%s/tests_%s/", out_fn, cti), showWarnings = FALSE)
	
	# binarise
	pr_ex_m = (pr_ex_m > 0) * 1
	ga_ex_m = (ga_ex_m > 0) * 1
	lo_ex_m = (lo_ex_m > 0) * 1

	# restrict genome conservation info to ogs of interest
	pr_ge_m_i = pr_ge_m [ rownames(ga_ex_m), ]
	ga_ge_m_i = ga_ge_m [ rownames(ga_ex_m), ]
	lo_ge_m_i = lo_ge_m [ rownames(ga_ex_m), ]
	
	# matrix of complex events
	ev_ex_a = as.matrix(data.frame(
		expressed = colSums(pr_ex_m),
		newexp_oldgene = colSums(ga_ex_m==1 & ga_ge_m_i==0, na.rm = TRUE),
		newexp_newgene = colSums(ga_ex_m==1 & ga_ge_m_i==1, na.rm = TRUE),
		lossexp_presgene = colSums(lo_ex_m==1 & pr_ge_m_i==1, na.rm = TRUE),
		lossexp_lossgene = colSums(lo_ex_m==1 & pr_ge_m_i==0, na.rm = TRUE)
	))
	# drop gains of expression/conservation in ancestral node
	ev_ex_a [ nrow(ev_ex_a), "newexp_newgene" ] = 0
	
	# dataframes for plotting over the phylogeny
	df_ga = ev_ex_a[,c("newexp_newgene","newexp_oldgene","expressed")]
	df_ga[,"expressed"] = df_ga[,"expressed"] - df_ga[,"newexp_newgene"] - df_ga[,"newexp_oldgene"]
	df_ga = df_ga + 0.01
	co_ga = c("chartreuse4","olivedrab2","lightcyan2")
	df_lo = ev_ex_a[,c("lossexp_presgene","lossexp_lossgene")]+0.01
	co_lo = c("indianred1","indianred4")
	
	pdf(sprintf("%s/anc.%s.mat.possvm.tree.pdf", out_fn, cti), width = 6, height = 5)
	ape::plot.phylo(sps_tree, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex = 1, show.tip.label = FALSE, show.node.label = FALSE)
	ape::nodelabels(pie = df_ga[sps_tree$node.label,], piecol = co_ga, cex = sqrt(rowSums(df_ga[sps_tree$node.label,]) / 100))
	ape::tiplabels(pie = df_ga[sps_tree$tip.label,],   piecol = co_ga, cex = sqrt(rowSums(df_ga[sps_tree$tip.label,]) / 100))
	ape::nodelabels(pie = df_lo[sps_tree$node.label,], piecol = co_lo, cex = sqrt(rowSums(df_lo[sps_tree$node.label,]) / 100), adj = 0.1)
	ape::tiplabels(pie = df_lo[sps_tree$tip.label,],  piecol = co_lo, cex = sqrt(rowSums(df_lo[sps_tree$tip.label,]) / 100), adj = 0.1)
	ape::nodelabels(
		text = sprintf("%s\np=%i\n+%i(%i)\n-%i(%i)", sps_tree$node.label, ev_ex_a[sps_tree$node.label,"expressed"], colSums(ga_ex_m[,sps_tree$node.label]), floor(df_ga[sps_tree$node.label,"newexp_newgene"]), colSums(lo_ex_m[,sps_tree$node.label]), floor(df_lo[sps_tree$node.label,"lossexp_lossgene"])),
		col = scales::alpha("gray20", 0.9), frame="none", cex = 0.7)
	ape::tiplabels(
		text = sprintf("%s\np=%i\n+%i(%i)\n-%i(%i)", sps_tree$tip.label, ev_ex_a[sps_tree$tip.label,"expressed"], colSums(ga_ex_m[,sps_tree$tip.label]), floor(df_ga[sps_tree$tip.label,"newexp_newgene"]), colSums(lo_ex_m[,sps_tree$tip.label]), floor(df_lo[sps_tree$tip.label,"lossexp_lossgene"])),
		col = scales::alpha("gray20", 0.9), frame="none", cex = 0.7)
	title(main = sprintf("%s", cti))
	dev.off()
	
	
	# list of OGs involved in simple events
	for (noi in nod_list) {
		
		# ogs present, gained and lost
		ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
		ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
		ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
		write.table(ogs_pr_i, sprintf("%s/tests_%s/tst.%s.%s.all.pres.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		write.table(ogs_ga_i, sprintf("%s/tests_%s/tst.%s.%s.all.gain.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		write.table(ogs_lo_i, sprintf("%s/tests_%s/tst.%s.%s.all.loss.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
		
		# define species to use as ref to analyse ogs
		if (noi %in% sps_tree$tip.label) {
			spi_v = noi
		} else {
			spi_v = names(descendants_from_nodes[[noi]])
		}
		
		# functional enrichments
		for (spi_w in spi_v) {
			
			# functional enrichment of each gene set		
			go_annot = gsa_topgo_load_emapper(emapper_fn = sprintf("../data/reference/%s_ensembl.GO.csv", spi_w), index_col_GOs = 2)
			pf_annot = gsa_enrichment_load_pfam_list(pfam_architecture_file = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi_w))
			names(pf_annot) = oga_gtv[names(pf_annot)]
			names(go_annot) = oga_gtv[names(go_annot)]
			annotated_genes = unique(c(names(pf_annot), names(go_annot)))

			# load KO info in species of interest
			keg_ll = read.table(sprintf("results_pathways/kegg_map_to_%s.csv",  spi_w), header = TRUE, sep = "\t")
			
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
			write.table(keg_ldl, sprintf("results_pathways/kegg_gene_to_map_%s.csv", spi_w), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
			keg_lll = gsa_enrichment_load_pfam_list(sprintf("results_pathways/kegg_gene_to_map_%s.csv", spi_w), architecture_sep = ",")
			names(keg_lll) = gg_v [ names(keg_lll) ]


			# expand og to genes using ref sps
			gen_pr_d = oga [ oga$orthogroup %in% ogs_pr_i & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_ga_d = oga [ oga$orthogroup %in% ogs_ga_i & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_lo_d = oga [ oga$orthogroup %in% ogs_lo_i & oga$species %in% spi_w, c("gene","orthogroup") ]
			gen_pr_d = gen_pr_d [ gen_pr_d$gene %in% annotated_genes, ]
			gen_ga_d = gen_ga_d [ gen_ga_d$gene %in% annotated_genes, ]
			gen_lo_d = gen_lo_d [ gen_lo_d$gene %in% annotated_genes, ]
			# gen_pr_i = gen_pr_d [ !duplicated(gen_pr_d$orthogroup),"gene" ]
			# gen_ga_i = gen_ga_d [ !duplicated(gen_ga_d$orthogroup),"gene" ]
			# gen_lo_i = gen_lo_d [ !duplicated(gen_lo_d$orthogroup),"gene" ]
			
			# list of tests
			list_tests = list(
				gain = gen_ga_d, pres = gen_pr_d, loss = gen_lo_d
			)
			
			for (ttt in 1:length(list_tests)) {
			
				# fetch
				ttlab = names(list_tests)[ttt]
				gen_x_d = list_tests[[ttt]]
				message(sprintf("evo | %s | test %s on %s branch, from %s...", cti, ttlab, noi, spi_w))
				
				# list of genes
				gen_x_i = gen_x_d [ !duplicated(gen_x_d$orthogroup), "gene" ]
				ttgen = gen_x_i
				
				# test domain enrichment
				if (length(ttgen) > 0) {
					
					# output tables
					gen_x_d$orthogroup = factor(gen_x_d$orthogroup)
					gen_x_d$orthogroup_name = ogm_gv [ gen_x_d$gene ]
					gen_x_d$orthogroup_name [ gen_x_d$orthogroup_name == "" ] = NA
					# presence of OGs
					gen_x_d_a1 = aggregate(gene ~ orthogroup, gen_x_d, function(vv) paste(vv, collapse = ","), drop = FALSE)
					gen_x_d_a2 = aggregate(orthogroup_name ~ orthogroup, gen_x_d, function(vv) paste(unique(sort(vv)), collapse = ","), drop = FALSE)
					gen_x_d_a  = data.frame(orthogroup = gen_x_d_a2[,1], orthogroup_name = gen_x_d_a2[,2], gene = gen_x_d_a1[,2])
					gen_x_d_a  = gen_x_d_a [ order(gen_x_d_a$orthogroup_name), ]
					write.table(gen_x_d_a, sprintf("%s/tests_%s/tst.%s.%s.%s.%s.oglist.tsv", out_fn, cti, cti, noi, spi_w, ttlab), quote = FALSE, sep = "\t", row.names = FALSE)

					# KEGG enrichment
					kegg_enrichment = gsa_enrichment_hypergeometric_v2(
						annotation = keg_lll,
						genes_fg = ttgen,
						out_fn = sprintf("%s/tests_%s/tst.%s.%s.%s.%s.kegg.csv", out_fn, cti, cti, noi, spi_w, ttlab)
					)
					keg_ld_i = keg_ld [ keg_ld$gene %in% ttgen, ]
					keg_ld_i$kegg_name = factor(keg_ld_i$kegg_name, levels = kegg_enrichment[[1]]$annot)
					keg_ld_i = keg_ld_i [ order(keg_ld_i$kegg_name, keg_ld_i$ko, keg_ld_i$gene), ]
					
					# KEGG enrichment pvalue
					kegg_enrichment_pv_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$pval_adj)			
					kegg_enrichment_fg_v = dic_from_vecs(kegg_enrichment[[1]]$annot, kegg_enrichment[[1]]$freq_in_fg)	
					keg_ld_i$kegg_pvalue = kegg_enrichment_pv_v[keg_ld_i$kegg_name]		
					keg_ld_i$kegg_n_in_fg = kegg_enrichment_fg_v[keg_ld_i$kegg_name]		
					write.table(keg_ld_i, sprintf("%s/tests_%s/tst.%s.%s.%s.%s.keggsum.csv", out_fn, cti, cti, noi, spi_w, ttlab), sep = "\t", quote = FALSE, row.names = FALSE)
					
					# message(sprintf("pathways %s | %s test %s pfam...", noi, cti, ttlab))
					pf_enrichment_l = gsa_enrichment_hypergeometric_v2(
						annotation = pf_annot,
						genes_fg = ttgen,
						out_fn = sprintf("%s/tests_%s/tst.%s.%s.%s.%s.pfam.csv", out_fn, cti, cti, noi, spi_w, ttlab)
					)

					# test GOs
					# message(sprintf("pathways %s | %s test %s GOs...", noi, cti, ttlab))
					go_enrichment_l = gsa_topgo_enrichment_v2(
						annotation = go_annot,
						genes_fg = ttgen,
						out_fn = sprintf("%s/tests_%s/tst.%s.%s.%s.%s.topgo.csv", out_fn, cti, cti, noi, spi_w, ttlab),
						ontologyset = c("BP","MF","CC"),
						tg_test = "fisher",
						tg_algorithm = "elim",
						top_markers = 30,
						nodesize = 10
					)
					
					# integrate GOs and pfam in a single table
					# message(sprintf("pathways %s | %s merge %s GO+pfam...", noi, cti, ttlab))
					gsa_merge_tables_v2(
						main_enrichment = go_enrichment_l[[1]],
						main_mapping = go_enrichment_l[[2]],
						pfam_enrichment = pf_enrichment_l[[1]],
						pfam_mapping = pf_enrichment_l[[2]],
						output_fn = sprintf("%s/tests_%s/tst.%s.%s.%s.%s.summary.csv", out_fn, cti, cti, noi, spi_w, ttlab),
						extra_annot_dict = list(ogm_gv, oga_gov)
					) 
					
				}
			
			}

		}
		
	}
	
	
	# # load expression values of markers
	# mks_l = list()
	# umi_l = list()
	# mks_t = data.frame()
	# for (spi in sps_list_in_comp) {
		
	# 	message(sprintf("load %s", spi))
	# 	mks_l[[spi]] = read.table(sprintf("%s/mks.%s.%s.all.tsv", out_fn, cti, spi), sep = "\t", header = TRUE)
	# 	umi_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umicount.rds", spi, spi))
	# 	umi_i_v = rowSums(as.matrix(umi_l[[spi]][,ctv[[spi]]]))
	# 	names(umi_i_v) = gg_v [ names(umi_i_v) ]
	# 	mks_f = mks_l[[spi]] [ , c("gene","orthogroup","umi","umifrac","fc", "p_val_adj")]
	# 	# add nonexpressed markers
	# 	mks_f_nonexp = data.frame(gene = names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ], orthogroup = oga_gov [ names(umi_i_v) [ !names(umi_i_v) %in% mks_f$gene ] ])
	# 	mks_f_nonexp$umi = 0
	# 	mks_f_nonexp$umifrac = 0
	# 	mks_f_nonexp$p_val_adj = 1
	# 	mks_f_nonexp$fc = 0
	# 	mks_f = rbind(mks_f, mks_f_nonexp)
	# 	# is significant?
	# 	mks_f$significant = mks_f$p_val_adj < thr_pv & mks_f$fc >= thr_fps
	
	# 	# add og (as gene name or og if absent)
	# 	mks_f$gene_name = ogm_gv [ mks_f$gene ]
	# 	mks_f = mks_f [ !is.na(mks_f$gene_name), ]
	# 	mks_f = mks_f [ order(mks_f$orthogroup, -as.numeric(mks_f$significant), mks_f$p_val_adj, -mks_f$fc) , ]
	# 	# mks_f = mks_f [ !duplicated(mks_f$orthogroup), ]
	# 	mks_f$species = spi
	# 	mks_t = rbind(mks_t, mks_f)
	
	# }
	# mks_t$species = factor(mks_t$species, levels = sps_list_in_comp)
	# mks_t_f = mks_t 
	# mks_t_f = mks_t_f [ order(mks_t_f$orthogroup, mks_t_f$species), ]
	# mks_t_f [ mks_t_f$gene_name == "Gale", ]
	# # mks_t_f [ mks_t_f$gene_name == "Galk2", ]
	# # mks_t_f [ mks_t_f$gene_name == "Galt", ]

		
	# }
	# # list of OGs involved in simple events
	# for (noi in nod_list) {
		
	# 	# ogs present, gained and lost
	# 	ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
	# 	ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
	# 	ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
	# 	write.table(ogs_pr_i, sprintf("%s/tests_%s/tst.%s.%s.all.pres.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	# 	write.table(ogs_ga_i, sprintf("%s/tests_%s/tst.%s.%s.all.gain.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
	# 	write.table(ogs_lo_i, sprintf("%s/tests_%s/tst.%s.%s.all.loss.oglist.txt", out_fn, cti, cti, noi), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


	# }
		
	# # list of OGs involved in simple events
	# for (noi in nod_list) {
		
	# 	# ogs present, gained and lost
	# 	ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
	# 	ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
	# 	ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
		
	# 	message(sprintf("aggregated expression | %s | %s, load...", cti, noi))
	# 	# expression profile of genes involved in each event (gain/loss/presence)
	# 	gen_pr_l = list()
	# 	gen_ga_l = list()
	# 	gen_lo_l = list()
	# 	ufs_l = list()
	# 	for (spi in sps_list_in_comp) {
			
	# 		if (spi == "Spin") {
	# 			spi_w = "Spis"
	# 		} else {
	# 			spi_w = spi
	# 		}
			
	# 		ufs_i = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umifrac.rds", spi, spi))
	# 		fps_i = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi, spi))
	# 		fcc_i = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_fracell.rds", spi, spi))
	# 		fps_i_v = apply(as.matrix(fps_i [ ,ctv[[spi]] ]), 1, function(vv) { max(vv)  } )
	# 		ufs_i_v = apply(as.matrix(ufs_i [ ,ctv[[spi]] ]), 1, function(vv) { max(vv)  } )
	# 		fcc_i_v = apply(as.matrix(fcc_i [ ,ctv[[spi]] ]), 1, function(vv) { max(vv)  } )
	# 		# ufs_l[[spi]] = ufs_i_v

	# 		# expand ogs to genes
	# 		gen_pr_l[[spi]] = oga [ oga$orthogroup %in% ogs_pr_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		gen_ga_l[[spi]] = oga [ oga$orthogroup %in% ogs_ga_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		gen_lo_l[[spi]] = oga [ oga$orthogroup %in% ogs_lo_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		# store expression: get gene ids
	# 		gen_pr_l[[spi]]$gene = gsub("_", "-", gen_pr_l[[spi]]$gene)
	# 		gen_ga_l[[spi]]$gene = gsub("_", "-", gen_ga_l[[spi]]$gene)
	# 		gen_lo_l[[spi]]$gene = gsub("_", "-", gen_lo_l[[spi]]$gene)
	# 		# store expression: fps
	# 		gen_pr_l[[spi]]$fp = fps_i_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$fp = fps_i_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$fp = fps_i_v [ gen_lo_l[[spi]]$gene ]
	# 		# store expression: umifrac
	# 		gen_pr_l[[spi]]$umifrac = ufs_i_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$umifrac = ufs_i_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$umifrac = ufs_i_v [ gen_lo_l[[spi]]$gene ]
	# 		# store expression: fracell
	# 		gen_pr_l[[spi]]$fracell = ufs_i_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$fracell = ufs_i_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$fracell = ufs_i_v [ gen_lo_l[[spi]]$gene ]
	# 		# # drop genes detected in few cells
	# 		# gen_pr_l[[spi]] = gen_pr_l[[spi]] [ gen_pr_l[[spi]]$fracell >= 5, ]
	# 		# gen_ga_l[[spi]] = gen_ga_l[[spi]] [ gen_ga_l[[spi]]$fracell >= 5, ]
	# 		# gen_lo_l[[spi]] = gen_lo_l[[spi]] [ gen_lo_l[[spi]]$fracell >= 5, ]
	# 		# keep one gene per OG (always max?)
	# 		gen_pr_l[[spi]] = gen_pr_l[[spi]] [ order(gen_pr_l[[spi]]$orthogroup, -gen_pr_l[[spi]]$umifrac) , ] 
	# 		gen_ga_l[[spi]] = gen_ga_l[[spi]] [ order(gen_ga_l[[spi]]$orthogroup, -gen_ga_l[[spi]]$umifrac) , ] 
	# 		gen_lo_l[[spi]] = gen_lo_l[[spi]] [ order(gen_lo_l[[spi]]$orthogroup, -gen_lo_l[[spi]]$umifrac) , ] 
	# 		gen_pr_l[[spi]] = gen_pr_l[[spi]] [ !duplicated(gen_pr_l[[spi]]$orthogroup) , ] 
	# 		gen_ga_l[[spi]] = gen_ga_l[[spi]] [ !duplicated(gen_ga_l[[spi]]$orthogroup) , ] 
	# 		gen_lo_l[[spi]] = gen_lo_l[[spi]] [ !duplicated(gen_lo_l[[spi]]$orthogroup) , ] 
	
	# 	}
		
	# 	# plot fps
	# 	message(sprintf("aggregated expression | %s | %s, fps...", cti, noi))
	# 	gen_pr_l_v = sapply(1:length(gen_pr_l), function(vv) { gen_pr_l[[vv]]$fp })
	# 	gen_ga_l_v = sapply(1:length(gen_ga_l), function(vv) { gen_ga_l[[vv]]$fp })
	# 	gen_lo_l_v = sapply(1:length(gen_lo_l), function(vv) { gen_lo_l[[vv]]$fp })
	# 	names(gen_pr_l_v) = names(gen_pr_l)
	# 	names(gen_ga_l_v) = names(gen_ga_l)
	# 	names(gen_lo_l_v) = names(gen_lo_l)
	# 	# pdf
	# 	pdf(sprintf("%s/exp.%s.%s.fps.pdf", out_fn, cti, noi), width = 3, height = 10)
	# 	layout(1:3)
	# 	boxplot(gen_pr_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "Footprint", names = sprintf("%s\n%i genes", names(gen_pr_l_v), lengths(gen_pr_l_v)))
	# 	title(main = sprintf("%s pres\nn=%i sps", noi, length(gen_pr_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_ga_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "Footprint", names = sprintf("%s\n%i genes", names(gen_ga_l_v), lengths(gen_ga_l_v)))
	# 	title(main = sprintf("%s gain\nn=%i sps", noi, length(gen_ga_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_lo_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "Footprint", names = sprintf("%s\n%i genes", names(gen_lo_l_v), lengths(gen_lo_l_v)))
	# 	title(main = sprintf("%s loss\nn=%i sps", noi, length(gen_lo_l_v)), cex.main = 1, font.main = 1)
	# 	dev.off()
		
	# 	# plot umifrac
	# 	message(sprintf("aggregated expression | %s | %s, umifrac...", cti, noi))
	# 	gen_pr_l_v = sapply(1:length(gen_pr_l), function(vv) { gen_pr_l[[vv]]$umifrac })
	# 	gen_ga_l_v = sapply(1:length(gen_ga_l), function(vv) { gen_ga_l[[vv]]$umifrac })
	# 	gen_lo_l_v = sapply(1:length(gen_lo_l), function(vv) { gen_lo_l[[vv]]$umifrac })
	# 	names(gen_pr_l_v) = names(gen_pr_l)
	# 	names(gen_ga_l_v) = names(gen_ga_l)
	# 	names(gen_lo_l_v) = names(gen_lo_l)
	# 	# pdf
	# 	pdf(sprintf("%s/exp.%s.%s.umifrac.pdf", out_fn, cti, noi), width = 3, height = 10)
	# 	layout(1:3)
	# 	boxplot(gen_pr_l_v, outline = FALSE, ylim = c(0,10), las = 2, ylab = "UMI/10k", names = sprintf("%s\n%i genes", names(gen_pr_l_v), lengths(gen_pr_l_v)))
	# 	title(main = sprintf("%s pres\nn=%i sps", noi, length(gen_pr_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_ga_l_v, outline = FALSE, ylim = c(0,10), las = 2, ylab = "UMI/10k", names = sprintf("%s\n%i genes", names(gen_ga_l_v), lengths(gen_ga_l_v)))
	# 	title(main = sprintf("%s gain\nn=%i sps", noi, length(gen_ga_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_lo_l_v, outline = FALSE, ylim = c(0,10), las = 2, ylab = "UMI/10k", names = sprintf("%s\n%i genes", names(gen_lo_l_v), lengths(gen_lo_l_v)))
	# 	title(main = sprintf("%s loss\nn=%i sps", noi, length(gen_lo_l_v)), cex.main = 1, font.main = 1)
	# 	dev.off()
		
	# 	# scatter plots
	# 	pdf(sprintf("%s/exp.%s.%s.scatterfp.pdf", out_fn, cti, noi), width = 8, height = 16)

	# 	message(sprintf("aggregated expression | %s | %s, scatter, pres...", cti, noi))
	# 	for (spi in names(gen_pr_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_pr_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_pr_l[[spi]][,c("orthogroup","gene","fp")], gen_pr_l[[spj]][,c("orthogroup","gene","fp")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0,5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}

	# 	message(sprintf("aggregated expression | %s | %s scatter, gain...", cti, noi))
	# 	for (spi in names(gen_ga_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_ga_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_ga_l[[spi]][,c("orthogroup","gene","fp")], gen_ga_l[[spj]][,c("orthogroup","gene","fp")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
		
	# 	message(sprintf("aggregated expression | %s | %s scatter, loss...", cti, noi))
	# 	for (spi in names(gen_lo_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_lo_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_lo_l[[spi]][,c("orthogroup","gene","fp")], gen_lo_l[[spj]][,c("orthogroup","gene","fp")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
	# 	dev.off()
		
		
	# 	# scatter plots
	# 	pdf(sprintf("%s/exp.%s.%s.scatteruf.pdf", out_fn, cti, noi), width = 8, height = 16)

	# 	message(sprintf("aggregated expression | %s | %s, scatter, pres...", cti, noi))
	# 	for (spi in names(gen_pr_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_pr_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_pr_l[[spi]][,c("orthogroup","gene","umifrac")], gen_pr_l[[spj]][,c("orthogroup","gene","umifrac")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 20 ] = 20
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 20), xlim = c(0, 20), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}

	# 	message(sprintf("aggregated expression | %s | %s scatter, gain...", cti, noi))
	# 	for (spi in names(gen_ga_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_ga_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_ga_l[[spi]][,c("orthogroup","gene","umifrac")], gen_ga_l[[spj]][,c("orthogroup","gene","umifrac")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 20 ] = 20
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 20), xlim = c(0, 20), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
		
	# 	message(sprintf("aggregated expression | %s | %s scatter, loss...", cti, noi))
	# 	for (spi in names(gen_lo_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_lo_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_lo_l[[spi]][,c("orthogroup","gene","umifrac")], gen_lo_l[[spj]][,c("orthogroup","gene","umifrac")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 20 ] = 20
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 20), xlim = c(0, 20), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 20), breaks = seq(0, 20, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
	# 	dev.off()
		
	# }

	# # based on markers fc
	# seu_l = list()
	# ufs_l = list()
	# # load seurat and umifrac
	# for (spi in sps_list_in_comp) {
			
	# 	if (spi == "Spin") {
	# 		spi_w = "Spis"
	# 	} else {
	# 		spi_w = spi
	# 	}
		
	# 	ufs_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.cts_umifrac.rds", spi, spi))
	# 	seu_l[[spi]] = readRDS(sprintf("results_metacell_%s_filt/dat.%s.seurat_final.rds", spi, spi))
		
	# }
	
	# # get markers
	# mks_l = list()
	# for (spi in sps_list_in_comp) {
	# 	SeuratObject::Idents(seu_l[[spi]]) = seu_l[[spi]]@meta.data[,"cell_type"]
	# 	mks_i = Seurat::FindMarkers(seu_l[[spi]], ident.1 = ctv[[spi]], layer = "RNA", min.pct = 0.05, test.use = "wilcox")
	# 	mks_i$gene = gg_v [ rownames(mks_i) ]
	# 	mks_i = mks_i [ !is.na(mks_i$gene), ]
	# 	mks_i$orthogroup = oga_gov [ mks_i$gene ]
	# 	mks_i$gene_name = ogm_gv [ mks_i$gene ]
	# 	mks_l[[spi]] = mks_i
	# }

	# for (noi in nod_list) {
		
	# 	# ogs present, gained and lost
	# 	ogs_pr_i = names(which(pr_ex_m[,noi] == 1))
	# 	ogs_ga_i = names(which(ga_ex_m[,noi] == 1))
	# 	ogs_lo_i = names(which(lo_ex_m[,noi] == 1))
		
	# 	message(sprintf("aggregated expression | %s | %s, load...", cti, noi))

	# 	# expression profile of genes involved in each event (gain/loss/presence)
	# 	gen_pr_l = list()
	# 	gen_ga_l = list()
	# 	gen_lo_l = list()
	# 	for (spi in sps_list_in_comp) {
			
	# 		if (spi == "Spin") {
	# 			spi_w = "Spis"
	# 		} else {
	# 			spi_w = spi
	# 		}
			
	# 		# expand ogs to genes
	# 		gen_pr_l[[spi]] = oga [ oga$orthogroup %in% ogs_pr_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		gen_ga_l[[spi]] = oga [ oga$orthogroup %in% ogs_ga_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		gen_lo_l[[spi]] = oga [ oga$orthogroup %in% ogs_lo_i & oga$species %in% spi_w, c("gene","orthogroup") ]
	# 		# store expression: from mks
	# 		mks_i_fc_v = dic_from_vecs(mks_l[[spi]]$gene, 2^mks_l[[spi]]$avg_log2FC)
	# 		mks_i_fr_v = dic_from_vecs(mks_l[[spi]]$gene, mks_l[[spi]]$pct.1)
	# 		mks_i_pv_v = dic_from_vecs(mks_l[[spi]]$gene, mks_l[[spi]]$p_val_adj)
	# 		gen_pr_l[[spi]]$fc         = mks_i_fc_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_pr_l[[spi]]$fracell    = mks_i_fr_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_pr_l[[spi]]$pval       = mks_i_pv_v [ gen_pr_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$fc         = mks_i_fc_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$fracell    = mks_i_fr_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_ga_l[[spi]]$pval       = mks_i_pv_v [ gen_ga_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$fc         = mks_i_fc_v [ gen_lo_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$fracell    = mks_i_fr_v [ gen_lo_l[[spi]]$gene ]
	# 		gen_lo_l[[spi]]$pval       = mks_i_pv_v [ gen_lo_l[[spi]]$gene ]
	# 		# keep one gene per OG (always max?)
	# 		gen_pr_l[[spi]]$fc [ is.na(gen_pr_l[[spi]]$fc) ] = 0
	# 		gen_ga_l[[spi]]$fc [ is.na(gen_ga_l[[spi]]$fc) ] = 0
	# 		gen_lo_l[[spi]]$fc [ is.na(gen_lo_l[[spi]]$fc) ] = 0
	# 		gen_pr_l[[spi]]$pval [ is.na(gen_pr_l[[spi]]$pval) ] = 1
	# 		gen_ga_l[[spi]]$pval [ is.na(gen_ga_l[[spi]]$pval) ] = 1
	# 		gen_lo_l[[spi]]$pval [ is.na(gen_lo_l[[spi]]$pval) ] = 1
	# 		gen_pr_l[[spi]] = gen_pr_l[[spi]] [ order(gen_pr_l[[spi]]$orthogroup, -gen_pr_l[[spi]]$fc) , ] 
	# 		gen_ga_l[[spi]] = gen_ga_l[[spi]] [ order(gen_ga_l[[spi]]$orthogroup, -gen_ga_l[[spi]]$fc) , ] 
	# 		gen_lo_l[[spi]] = gen_lo_l[[spi]] [ order(gen_lo_l[[spi]]$orthogroup, -gen_lo_l[[spi]]$fc) , ] 
	# 		gen_pr_l[[spi]] = gen_pr_l[[spi]] [ !duplicated(gen_pr_l[[spi]]$orthogroup) , ] 
	# 		gen_ga_l[[spi]] = gen_ga_l[[spi]] [ !duplicated(gen_ga_l[[spi]]$orthogroup) , ] 
	# 		gen_lo_l[[spi]] = gen_lo_l[[spi]] [ !duplicated(gen_lo_l[[spi]]$orthogroup) , ] 

	# 	}
		
		
	# 	# plot fps
	# 	message(sprintf("aggregated expression | %s | %s, fcs...", cti, noi))
	# 	gen_pr_l_v = sapply(1:length(gen_pr_l), function(vv) { gen_pr_l[[vv]]$fc })
	# 	gen_ga_l_v = sapply(1:length(gen_ga_l), function(vv) { gen_ga_l[[vv]]$fc })
	# 	gen_lo_l_v = sapply(1:length(gen_lo_l), function(vv) { gen_lo_l[[vv]]$fc })
	# 	names(gen_pr_l_v) = names(gen_pr_l)
	# 	names(gen_ga_l_v) = names(gen_ga_l)
	# 	names(gen_lo_l_v) = names(gen_lo_l)
	# 	# pdf
	# 	pdf(sprintf("%s/exp.%s.%s.fcs.pdf", out_fn, cti, noi), width = 3, height = 10)
	# 	layout(1:3)
	# 	boxplot(gen_pr_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "FC", names = sprintf("%s\n%i genes", names(gen_pr_l_v), lengths(gen_pr_l_v)))
	# 	title(main = sprintf("%s pres\nn=%i sps", noi, length(gen_pr_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_ga_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "FC", names = sprintf("%s\n%i genes", names(gen_ga_l_v), lengths(gen_ga_l_v)))
	# 	title(main = sprintf("%s gain\nn=%i sps", noi, length(gen_ga_l_v)), cex.main = 1, font.main = 1)
	# 	boxplot(gen_lo_l_v, outline = FALSE, ylim = c(0,6), las = 2, ylab = "FC", names = sprintf("%s\n%i genes", names(gen_lo_l_v), lengths(gen_lo_l_v)))
	# 	title(main = sprintf("%s loss\nn=%i sps", noi, length(gen_lo_l_v)), cex.main = 1, font.main = 1)
	# 	dev.off()
		
	# 	# scatter plots
	# 	pdf(sprintf("%s/exp.%s.%s.scatterfc.pdf", out_fn, cti, noi), width = 8, height = 16)

	# 	message(sprintf("aggregated expression | %s | %s, scatter, pres...", cti, noi))
	# 	for (spi in names(gen_pr_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_pr_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_pr_l[[spi]][,c("orthogroup","gene","fc")], gen_pr_l[[spj]][,c("orthogroup","gene","fc")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("pres %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}

	# 	message(sprintf("aggregated expression | %s | %s scatter, gain...", cti, noi))
	# 	for (spi in names(gen_ga_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_ga_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_ga_l[[spi]][,c("orthogroup","gene","fc")], gen_ga_l[[spj]][,c("orthogroup","gene","fc")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("gain %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
		
	# 	message(sprintf("aggregated expression | %s | %s scatter, loss...", cti, noi))
	# 	for (spi in names(gen_lo_l)) {
	# 		layout(matrix(1:8, ncol = 2, byrow = 2))
	# 		for (spj in names(gen_lo_l)) {
	# 			if (spi != spj) {
	# 				mm = merge(gen_lo_l[[spi]][,c("orthogroup","gene","fc")], gen_lo_l[[spj]][,c("orthogroup","gene","fc")], by = "orthogroup", all.x = FALSE, all.y = FALSE, suffixes = c(spi, spj))
	# 				mmm = as.matrix(mm[,c(3,5)])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, col = "blue", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 					# histogram
	# 					mmm1 = mmm[,1]
	# 					mmm1 [ mmm1 > 10 ] = 10
	# 					hist(mmm1, col = scales::alpha("orchid2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1))
	# 					mmm2 = mmm[,2]
	# 					mmm2 [ mmm2 > 10 ] = 10
	# 					hist(mmm2, col = scales::alpha("orange2", 0.5), xlim = c(0, 5), breaks = seq(0, 5, 0.1), add = TRUE)
	# 					legend("topright", c(spi, spj), fill = scales::alpha(c("orchid2","orange2"), 0.5), cex = 1, bty = "n")
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("loss %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 			}
	# 		}
	# 	}
	# 	dev.off()


	# 	# scatter plots
	# 	pdf(sprintf("%s/exp.%s.%s.scatterfcall.pdf", out_fn, cti, noi), width = 9, height = 10)

	# 	message(sprintf("aggregated expression | %s | %s, scatter, pres...", cti, noi))
	# 	for (spi in names(gen_pr_l)) {
	# 		layout(matrix(1:9, ncol = 3, byrow = TRUE))
	# 		for (spj in names(gen_pr_l)) {
	# 			if (spi != spj) {
	# 				# mm = merge(mks_l[[spi]], mks_l[[spj]], by = "orthogroup", all.x = FALSE, all.y = FALSE)
	# 				# mmm = as.matrix(mm[,c("avg_log2FC.x","avg_log2FC.y")])
	# 				# mmm = 2 ^ mmm
	# 				# mmm [ is.na(mmm) ] = 0
	# 				# mmm [ mmm > 10 ] = 10
	# 				# mms = -1 * log10(mm$p_val_adj.x)
	# 				# mmc = cut(mms, c(0, 3, 6, 9, Inf), include.lowest = TRUE)
	# 				# # mmc = c(0.5, 1, 1.5, 2) [ mmc ]
	# 				# mmc = 1
	# 				# if (nrow(mmm)>1) {
	# 				# 	plot(mmm, xlab = spi, ylab = spj, pch = 19, cex = mmc, col = "lightcyan3", ylim = c(0, max(mmm)), xlim = c(0, max(mmm)), las = 1)
	# 				# 	abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 				# 	title(main = sprintf("p/g/l %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				# } else {
	# 				# 	plot(0)
	# 				# 	title(main = sprintf("p/g/l %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				# }
	# 				mm = merge(gen_ga_l[[spi]], gen_ga_l[[spj]], by = "orthogroup", all.x = FALSE, all.y = FALSE)
	# 				mmm = as.matrix(mm[,c("fc.x","fc.y")])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				mms = -1 * log10(mm$pval.x)
	# 				mmc = cut(mms, c(0, 3, 6, 9, Inf), include.lowest = TRUE)
	# 				mmc = c(0.5, 1, 1, 1) [ mmc ]
	# 				if (nrow(mmm)>1) {
	# 					plot(mmm, xlab = spi, ylab = spj, pch = 19, cex = mmc, col = "chartreuse3", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					abline(v = 1.1, h = 1.1, lty = 2, col = "orangered4")
	# 					title(main = sprintf("g/l %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				} else {
	# 					plot(0)
	# 					title(main = sprintf("g/l %s\n%s~%s\nn=%i genes", noi, spi, spj, nrow(mm)), cex.main = 1, font.main = 1)
	# 				}
	# 				mm = merge(gen_lo_l[[spi]], gen_lo_l[[spj]], by = "orthogroup", all.x = FALSE, all.y = FALSE)
	# 				mmm = as.matrix(mm[,c("fc.x","fc.y")])
	# 				mmm [ is.na(mmm) ] = 0
	# 				mmm [ mmm > 5 ] = 5
	# 				mms = -1 * log10(mm$pval.x)
	# 				mmc = cut(mms, c(0, 3, 6, 9, Inf), include.lowest = TRUE)
	# 				mmc = c(0.5, 1, 1, 1) [ mmc ]
	# 				if (nrow(mmm)>1) {
	# 					points(mmm, xlab = spi, ylab = spj, pch = 19, cex = mmc, col = "indianred1", ylim = c(0, 5), xlim = c(0, 5), las = 1)
	# 					abline(a = 0, b = 1, lty = 2, col = "darkred")
	# 					abline(v = 1.1, h = 1.1, lty = 2, col = "orangered4")
	# 				} else {
	# 					plot(0)
	# 				}
	# 			}
	# 		}
	# 		plot(0)
	# 		legend("topright", col = c("chartreuse3","indianred1"), c(sprintf("gain in %s", noi), sprintf("loss in %s", noi)), pch = 19, bty = "n")
	# 		legend("bottomright", col = c("blue"), c(sprintf("padj>1e-3 %s", spi), sprintf("padj<1e-3 %s", spi)), pch = 19, cex = c(0.5, 1), bty = "n")
	# 	}
	# 	dev.off()
		
	# }

}

message("All done!")
