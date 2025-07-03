# libraries
suppressMessages(library("scales"))
suppressMessages(source("../scripts/Gene_module_functions.R"))
suppressMessages(source("../scripts/Cross_species_functions.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(library("UpSetR"))
graphics.off()

# data index, alga v gastrodermis modules
list_comparisons = list(
	Ocupat = list(focus = c("brown2","darkmagenta"), other = c("brown")),
	Amil = list(focus = c("brown"), other = c("black")),
	Spin = list(focus = c("blue"), other = c("greenyellow")),
	Xesp = list(focus = c("white","darkolivegreen"), other = c("turquoise")),
	Ocuarb = list(focus = c("firebrick4","lightsteelblue","yellowgreen"), other = c("ivory","maroon","red","turquoise","grey"))
)

# load orthology
oga = read.table("../data/orthology_Anthozoa_plus/orthogroup_conservation.csv", sep = "\t", header = TRUE, quote = "")
oga_gv = dic_from_vecs(names = oga$gene, terms = oga$orthogroup_name)
oga_tv = dic_from_vecs(names = oga$transcript, terms = oga$orthogroup_name)
oga_gtv = dic_from_vecs(names = oga$transcript, terms = oga$gene)
oga_v = dic_from_vecs(names = oga$gene, terms = oga$orthogroup)
gg_v = dic_from_vecs(gsub("_","-",oga$gene), oga$gene)

# load gene names
gna_fn = "../data/orthology_Metazoa_plus/orthogroup_conservation.csv"
gna = read.table(gna_fn, sep = "\t", header = TRUE)

out_fn = "results_host_cells/"
dir.create(out_fn, showWarnings = FALSE)

list_ogs_i = vector(mode = "list", length = length(list_comparisons))
list_ogs_j = vector(mode = "list", length = length(list_comparisons))
names(list_ogs_i) = names(list_ogs_j) = names(list_comparisons)

for (spi in names(list_comparisons)) {
	
	# input
	inp_fn = sprintf("results_metacell_%s_filt/", spi)

	# reload
	message(sprintf("wgcna | %s | reload WGCNA...", spi))
	mc_wgcna_gmods = readRDS(sprintf("%s/gmod.%s.wgcna.memberships.rds", inp_fn, spi))
	mc_wgcna_gmods_annot = read.table(sprintf("%s/gmod.%s.wgcna.memberships.csv", inp_fn, spi), sep = "\t", header = TRUE)
	mc_wgcna_me = readRDS(sprintf("%s/gmod.%s.wgcna.ME.rds", inp_fn, spi))
	colnames(mc_wgcna_me) = gsub("^ME","", colnames(mc_wgcna_me))
	mc_wgcna_gmods_annot$gene = gsub("_", "-", mc_wgcna_gmods_annot$gene)

	message(sprintf("wgcna | %s | load annotations...", spi))
	ctt = read.table(sprintf("%s/annot.%s.leiden.csv", inp_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	ctm = read.table(sprintf("%s/annot.%s.mcs.csv", inp_fn, spi), sep = "\t", header = TRUE, comment.char = "")
	ctm_gastrodermis_bool = grepl("gastrodermis", ctm$metacell_top_cell_type)
	ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
	ctt_bct_v = dic_from_vecs(ctt$cell_type, ctt$broad_cell_type)

	message(sprintf("wgcna | %s | load expression...", spi))
	umf = readRDS(sprintf("results_metacell_%s_filt/dat.%s.expression.mcs_umifrac.rds", spi, spi))

	pdf(sprintf("%s/alga_hosting_module_activity.%s.pdf", out_fn, spi), width = 12, height = 4)
	for (moi in list_comparisons[[spi]]$focus) {

		# eigenvalues i		
		ev_i = t(mc_wgcna_me)[moi,]
		
		# umifrac i
		ge_i = mc_wgcna_gmods_annot$gene [ mc_wgcna_gmods_annot$module == moi ]
		uf_i = colSums(umf[ge_i,])
		
		# store away ogs in module
		ge_i_ogs = as.character(oga_v[as.character(gg_v[ge_i])])
		list_ogs_i[[spi]] = unique(c(list_ogs_i[[spi]], ge_i_ogs))

		for (moj in list_comparisons[[spi]]$other) {
		
			message(sprintf("wgcna | %s | compare 1st:%s - 2nd:%s...", spi, moi, moj))
			
			# eigenvalues j		
			ev_j = t(mc_wgcna_me)[moj,]

			# umifrac i
			ge_j = mc_wgcna_gmods_annot$gene [ mc_wgcna_gmods_annot$module == moj ]
			uf_j = colSums(umf[ge_j,])

			# store away ogs in module
			ge_j_ogs = as.character(oga_v[as.character(gg_v[ge_j])])
			list_ogs_j[[spi]] = unique(c(list_ogs_j[[spi]], ge_j_ogs))
			
			layout(matrix(1:3, nrow = 1))

			# plot eigenval
			plot(
				ev_i,
				ev_j,
				col = ctm_mcs_col_v[rownames(mc_wgcna_me)],
				pch = 19, xlab = moi, ylab = moj, 
				main = sprintf("%s eigenvalues", spi),
				xlim = c(-0.1,max(c(ev_j, ev_i))),
				ylim = c(-0.1,max(c(ev_j, ev_i)))
			)
			abline(a=0, b=1, lty=2, col = "darkred")
			text(ev_i, ev_j, rownames(mc_wgcna_me), col = colorspace::darken(ctm_mcs_col_v[rownames(mc_wgcna_me)]), cex = 0.2)

			# plot umifrac
			plot(
				uf_i,
				uf_j,
				col = ctm_mcs_col_v[names(uf_i)],
				pch = 19, xlab = moi, ylab = moj, 
				main = sprintf("%s sum umifrac", spi),
				xlim = c(0,max(c(uf_i, uf_j))),
				ylim = c(0,max(c(uf_i, uf_j)))
			)
			abline(a=0, b=1, lty=2, col = "darkred")
			text(uf_i, uf_j, names(uf_i), col = colorspace::darken(ctm_mcs_col_v[names(uf_i)]), cex = 0.2)

			# plot eigenval in gastrodermis
			barplot(ev_i[ctm_gastrodermis_bool], col = ctm_mcs_col_v[rownames(mc_wgcna_me)[ctm_gastrodermis_bool]], border = NA, space = 0, main = sprintf("%s: %s", spi, moi), names.arg = NULL)
		
		}
	}
	dev.off()

}

xtd = UpSetR::fromList(list_ogs_i)
# plot
uu1 = UpSetR::upset(
	xtd,
	nsets = ncol(xtd),
	keep.order = TRUE,
	sets = colnames(xtd),
	cutoff = 0,
	order.by = "freq",
	mainbar.y.label = "shared families host cells",
	sets.x.label = "num families",
	matrix.color = "black",
	main.bar.color = "lightblue3", 
	sets.bar.color = "lightblue4"
)

pdf(sprintf("%s/overlaps.gmod.host.pdf", out_fn), width = 10, height = 5)
print(uu1)
dev.off()

xtd = UpSetR::fromList(list_ogs_j)
# plot
uu1 = UpSetR::upset(
	xtd,
	nsets = ncol(xtd),
	keep.order = TRUE,
	sets = colnames(xtd),
	cutoff = 0,
	order.by = "freq",
	mainbar.y.label = "shared families gastrodermal cells",
	sets.x.label = "num families",
	matrix.color = "black",
	main.bar.color = "lightblue3", 
	sets.bar.color = "lightblue4"
)

pdf(sprintf("%s/overlaps.gmod.gast.pdf", out_fn), width = 10, height = 5)
print(uu1)
dev.off()
