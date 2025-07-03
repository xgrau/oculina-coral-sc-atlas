library("rtracklayer")
library("zoo")
library("Matrix")
library("viridis")
source("../scripts/helper.R")
source("../scripts/Cross_species_functions.R")
graphics.off()

# out
out_fn = "results_macrosynteny_plus"

# define species lists
# sps_list = read.table("data/species_list_synteny_blocks.txt")[,1]
sps_tree = ape::read.tree("data/species_list_synteny_blocks_plus.newick")
sps_list = sps_tree$tip.label


# homology groups
hom_d = read.table(sprintf("%s/mcl.all.out.txt", out_fn), header = FALSE, col.names = c("homology_group", "transcript"))
hom_d$species = gsub("_.*","", hom_d$transcript)
# hom_d$species [ hom_d$species == "Hvul" ] = "Hvul_v3"
hom_d = hom_d [ hom_d$species %in% sps_list , ]
hom_d = hom_d [ !duplicated(hom_d$transcript), ]
rownames(hom_d) = hom_d$transcript
hom_d_v = dic_from_vecs(hom_d$transcript, hom_d$homology_group)

# compare
list_comparisons = list(
	hsrhxe = c("Hvul","Rhoesc","Xesp")
)


# color palettes
categorical_colors      = colorspace::darken(c("magenta4","violet","firebrick1","orange","khaki1","palegreen2","springgreen3","darkgreen"), 0.1)
categorical_colors_dark = colorspace::darken(c("deepskyblue","paleturquoise1","mediumblue","darkviolet","orchid1","thistle1"), 0.1)
color_palette           = colorRampPalette(categorical_colors)
color_palette_dark      = colorRampPalette(categorical_colors_dark)

dir.create(sprintf("%s/pairwise/", out_fn), showWarnings = FALSE)

list_gff = list()
for (spi in sps_list) {
	
	message(sprintf("preload | %s GFF...", spi))
	list_gff[[spi]] = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/%s_long.annot.gtf", spi))
	
}

for (n in 1:length(list_comparisons)) {
	
	cmpnm = names(list_comparisons)[n]

	# load ancestral linkage groups
	alg_d = read.table(sprintf("%s/alg.%s.csv", out_fn, cmpnm), sep = "\t", header = TRUE)
	list_algs = unique(alg_d$ancestral_linkage_group)
	
	# informative HGs
	list_hgs = unique(alg_d$homology_group)
	dict_hgs = dic_from_vecs(alg_d$homology_group, alg_d$ancestral_linkage_group)

	# color ALGs
	spr_colors = color_palette(n=length(list_algs))
	names(spr_colors) = unique(list_algs)

	# loop over queries	
	
	for (spi in sps_list) {
		
		# load
		message(sprintf("ALG %s | %s | parse...", cmpnm, spi))
		gen_i = list_gff[[spi]]
		gen_i = gen_i [ gen_i$type == "transcript" ]
		gen_i_d = as.data.frame(gen_i)[,c("seqnames","start","end","strand","transcript_id")]
		gen_i_d$homology_group = hom_d_v [ gen_i_d$transcript_id ]
		gen_i_d = gen_i_d [ !is.na(gen_i_d$homology_group),  ]
		gen_i_d$seqnames = paste(spi, gen_i_d$seqnames, sep = ":")
		# gen_i_d$seqnames = factor(gen_i_d$seqnames, levels = names(gix_i_v))
		gen_i_d = gen_i_d [ order(gen_i_d$seqnames, gen_i_d$start), ]

		# filter out chromosomes with few genes
		gen_i_d = gen_i_d [ gen_i_d$seqnames %in% names(which(table(gen_i_d$seqnames) > 200)),  ]

		# filter informative HGs
		gen_i_d = gen_i_d [ gen_i_d$homology_group %in% list_hgs,  ]

		pdf(sprintf("%s/pairwise/%s.ref_is_%s.pdf", out_fn, cmpnm, spi), width = 12, height = 3 * length(sps_list))
		layout(1:length(sps_list))
		for (nnj in 1:length(sps_list)) {

			# load			
			spj = sps_list[nnj]
			message(sprintf("ALG %s | %s | compare to %s...", cmpnm, spi, spj))
			gen_j = list_gff[[spj]]
			gen_j = gen_j [ gen_j$type == "transcript" ]
			gen_j_d = as.data.frame(gen_j)[,c("seqnames","start","end","strand","transcript_id")]
			gen_j_d$homology_group = hom_d_v [ gen_j_d$transcript_id ]
			gen_j_d = gen_j_d [ !is.na(gen_j_d$homology_group),  ]
			gen_j_d$seqnames = paste(spj, gen_j_d$seqnames, sep = ":")
			# gen_j_d$seqnames = factor(gen_j_d$seqnames, levels = names(gix_j_v))
			gen_j_d = gen_j_d [ order(gen_j_d$seqnames, gen_j_d$start), ]

			# filter out chromosomes with few genes
			gen_j_d = gen_j_d [ gen_j_d$seqnames %in% names(which(table(gen_j_d$seqnames) > 200)),  ]

			# filter informative HGs
			gen_j_d = gen_j_d [ gen_j_d$homology_group %in% list_hgs, ]
			
			# filter informative HGs
			message(sprintf("ALG %s | %s | compare to %s, filter HGs...", cmpnm, spi, spj))
			list_hgs_good_ij = intersect(
				names(which(table(gen_i_d$homology_group)<=1)),
				names(which(table(gen_j_d$homology_group)<=1))
			)
			gen_i_d_f = gen_i_d [ gen_i_d$homology_group %in% list_hgs_good_ij, ]
			gen_j_d_f = gen_j_d [ gen_j_d$homology_group %in% list_hgs_good_ij, ]
			
			# add ALG call
			gen_i_d_f$ALG = dict_hgs [ gen_i_d_f$homology_group ]
			gen_j_d_f$ALG = dict_hgs [ gen_j_d_f$homology_group ]
			gen_i_d_f$ALG = factor(gen_i_d_f$ALG, levels = list_algs)
			gen_j_d_f$ALG = factor(gen_j_d_f$ALG, levels = list_algs)
			

			# reorder chromosomes based on ALG alignment, i
			message(sprintf("ALG %s | %s | compare to %s, reorder i...", cmpnm, spi, spj))
			mat_i = table_to_matrix(table(gen_i_d_f$seqnames, gen_i_d_f$homology_group))
			ovs_i = matrix(nrow = length(list_algs), ncol = nrow(mat_i))
			for (ali in 1:length(list_algs)) {
				alg = list_algs[ali]
				hgs_in_alg = alg_d [ alg_d$ancestral_linkage_group == alg, "homology_group" ]
				hgs_in_alg = hgs_in_alg [ hgs_in_alg %in% hom_d [ hom_d$species == spi, "homology_group" ] ]
				mat_t_a = mat_i [ , hgs_in_alg [ hgs_in_alg %in% colnames(mat_i) ]]
				ovs_i [ ali , ] = rowSums(mat_t_a) / ncol(mat_t_a)
			}
			rownames(ovs_i) = list_algs
			colnames(ovs_i) = rownames(mat_i)
			# only do this once, so that all plots are in the same order
			if (nnj == 1) {
				ovs_i_f_chrom_order = unlist(unique(apply(ovs_i, 1, function(r) names(which.max(r)))))
				ovs_i_f_chrom_order = c(ovs_i_f_chrom_order, setdiff(colnames(ovs_i), ovs_i_f_chrom_order) )
			}

			# reorder chromosomes based on ALG alignment, j
			message(sprintf("ALG %s | %s | compare to %s, reorder j...", cmpnm, spi, spj))
			mat_j = table_to_matrix(table(gen_j_d_f$seqnames, gen_j_d_f$homology_group))
			ovs_j = matrix(nrow = length(list_algs), ncol = nrow(mat_j))
			for (ali in 1:length(list_algs)) {
				alg = list_algs[ali]
				hgs_in_alg = alg_d [ alg_d$ancestral_linkage_group == alg, "homology_group" ]
				hgs_in_alg = hgs_in_alg [ hgs_in_alg %in% hom_d [ hom_d$species == spi, "homology_group" ] ]
				mat_t_a = mat_j [ , hgs_in_alg [ hgs_in_alg %in% colnames(mat_j) ]]
				ovs_j [ ali , ] = rowSums(mat_t_a) / ncol(mat_t_a)
			}
			rownames(ovs_j) = list_algs
			colnames(ovs_j) = rownames(mat_j)
			ovs_j_f_chrom_order = unlist(unique(apply(ovs_j, 1, function(r) names(which.max(r)))))
			ovs_j_f_chrom_order = c(ovs_j_f_chrom_order, setdiff(colnames(ovs_j), ovs_j_f_chrom_order) )

			# match chromosomes by ALG sharedness relative to i
			chralg_m_i = table_to_matrix(table(gen_i_d_f$ALG, gen_i_d_f$seqnames))
			chralg_m_j = table_to_matrix(table(gen_j_d_f$ALG, gen_j_d_f$seqnames))
			chralg_m_ij_c = cor(chralg_m_i, chralg_m_j)
			chralg_m_ij_v = apply(chralg_m_ij_c, 1, function(vv) colnames(chralg_m_ij_c) [ which.max(vv) ])
			chralg_m_ji_v = apply(chralg_m_ij_c, 2, function(vv) rownames(chralg_m_ij_c) [ which.max(vv) ])

			# use this order as factor
			gen_i_d_f$seqnames = factor(gen_i_d_f$seqnames, levels = ovs_i_f_chrom_order)
			gen_j_d_f$seqnames = factor(gen_j_d_f$seqnames, levels = ovs_j_f_chrom_order)
			gen_i_d_f = gen_i_d_f [ order(gen_i_d_f$seqnames, gen_i_d_f$start), ]
			gen_j_d_f = gen_j_d_f [ order(gen_j_d_f$seqnames, gen_j_d_f$start), ]

			# # add chr length to position
			# gen_i_d_f$start_new = gen_i_d_f$start + gix_i_v [gen_i_d_f$seqnames]
			# gen_j_d_f$start_new = gen_j_d_f$start + gix_j_v [gen_j_d_f$seqnames]
			
			# # instead of position, calculate rank within chromosome
			# gen_i_d_f$rank_new = order(gen_i_d_f$start_new)
			# gen_j_d_f$rank_new = order(gen_j_d_f$start_new)

			# # sort by ALG and position
			# gen_i_d_f = gen_i_d_f [ order(gen_i_d_f$ALG, gen_i_d_f$start), ]
			# gen_j_d_f = gen_j_d_f [ order(gen_j_d_f$ALG, gen_j_d_f$start), ]
			
			# instead of position, calculate rank within chromosome
			gen_i_d_f$rank = NA
			gen_j_d_f$rank = NA
			mxr = 0
			for (chr in levels(gen_i_d_f$seqnames)) {
				mxo = order(gen_i_d_f$start [ gen_i_d_f$seqnames == chr ])
				mxo = scales::rescale(mxo, c(0, 100))
				gen_i_d_f$rank [ gen_i_d_f$seqnames == chr ] = mxr + mxo
				mxr = max(gen_i_d_f$rank, na.rm = TRUE) + 50
			}
			mxr = 0
			for (chr in levels(gen_j_d_f$seqnames)) {
				mxo = order(gen_j_d_f$start [ gen_j_d_f$seqnames == chr ])
				hmo = gen_j_d_f [ gen_j_d_f$seqnames == chr , "homology_group" ]
				mxo = scales::rescale(mxo, c(0, 100))
				# do we need to flip chromosome j to match i orientation?
				if (chr %in% names(chralg_m_ji_v)) {
					chr_i = chralg_m_ji_v[chr]
					mxo_i = gen_i_d_f [ gen_i_d_f$seqnames == chr_i , c("rank","homology_group") ]
					mxo_i_v = dic_from_vecs(mxo_i$homology_group, mxo_i$rank)
					mxo_i_c = stats::cor(mxo, mxo_i_v[hmo], method = "spearman", use = "complete.obs")
					if (!is.na(mxo_i_c) & mxo_i_c < -0.2) {
						mxo = scales::rescale(mxo, c(100, 0))
					} 
				}
				gen_j_d_f$rank [ gen_j_d_f$seqnames == chr ] = mxr + mxo
				mxr = max(gen_j_d_f$rank, na.rm = TRUE) + 50
			}
			
			message(sprintf("ALG %s | %s | compare to %s, merge...", cmpnm, spi, spj))
			gen_m_d_f = merge(gen_i_d_f, gen_j_d_f, by = "homology_group", suffixes = c("_i","_j"))
			gen_m_d_f = gen_m_d_f [ order(gen_m_d_f$seqnames_i, gen_m_d_f$start_i), ]
			
			# rescale ranks globally?
			gen_m_d_f$rank_i = scales::rescale(gen_m_d_f$rank_i, c(0,2000))
			gen_m_d_f$rank_j = scales::rescale(gen_m_d_f$rank_j, c(0,2000))
			
			message(sprintf("ALG %s | %s | compare to %s, aggregate...", cmpnm, spi, spj))
			gen_m_d_f_a_i = aggregate(rank_i ~ seqnames_i, data = gen_m_d_f[,c("seqnames_i","rank_i")], min)
			gen_m_d_f_a_j = aggregate(rank_j ~ seqnames_j, data = gen_m_d_f[,c("seqnames_j","rank_j")], min)
			colnames(gen_m_d_f_a_i)[2] = "start_rank_chr"
			colnames(gen_m_d_f_a_j)[2] = "start_rank_chr"
			gen_m_d_f_a_i$end_rank_chr = aggregate(rank_i ~ seqnames_i, data = gen_m_d_f[,c("seqnames_i","rank_i")], max)[,2]
			gen_m_d_f_a_j$end_rank_chr = aggregate(rank_j ~ seqnames_j, data = gen_m_d_f[,c("seqnames_j","rank_j")], max)[,2]

			# complex labels, i
			gen_m_d_f_a_i_t = aggregate(ALG_i ~ seqnames_i, data = gen_m_d_f[,c("seqnames_i","ALG_i")], function(vv) { 
				vt = sort(table(vv), decreasing = TRUE)
				vt = vt [ vt >= 5]
				vl = paste(sprintf("%s (%i)", names(vt), vt), collapse = "\n")
			})
			gen_m_d_f_a_i_tt = apply(gen_m_d_f_a_i_t, 1, function(vv) paste(vv, collapse = "\n"))

			# complex labels, j
			gen_m_d_f_a_j_t = aggregate(ALG_j ~ seqnames_j, data = gen_m_d_f[,c("seqnames_j","ALG_j")], function(vv) { 
				vt = sort(table(vv), decreasing = TRUE)
				vt = vt [ vt >= 5 ]
				vl = paste(sprintf("%s (%i)", names(vt), vt), collapse = "\n")
			})
			gen_m_d_f_a_j_tt = apply(gen_m_d_f_a_j_t, 1, function(vv) paste(vv, collapse = "\n"))

			# plot
			message(sprintf("ALG %s | %s | compare to %s, plot...", cmpnm, spi, spj))
			plot(NA, xlim = c(0,max(c(2400,gen_m_d_f$rank_i, gen_m_d_f$rank_j))), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
			title(main = sprintf("%s\n%s (%i) - %s (%i)", cmpnm, spi, nlevels(gen_m_d_f$seqnames_i), spj, nlevels(gen_m_d_f$seqnames_j)), cex.main = 1, font.main = 1)
			segments(x0 = gen_m_d_f$rank_i, y0 = 0.8, x1 = gen_m_d_f$rank_j, y1 = 0.2, col = spr_colors[gen_m_d_f$ALG_i])
			segments(x0 = gen_m_d_f_a_i$start_rank_chr, x1 = gen_m_d_f_a_i$end_rank_chr, y0 = 0.8, y1 = 0.8, lwd = 4, col = "gray50")
			segments(x0 = gen_m_d_f_a_j$start_rank_chr, x1 = gen_m_d_f_a_j$end_rank_chr, y0 = 0.2, y1 = 0.2, lwd = 4, col = "gray50")
			text(x = gen_m_d_f_a_i$start_rank_chr, y = 1.00, gen_m_d_f_a_i_tt, col = "gray5", cex = 0.5, adj = c(0, 1))
			text(x = gen_m_d_f_a_j$start_rank_chr, y = 0.15, gen_m_d_f_a_j_tt, col = "gray5", cex = 0.5, adj = c(0, 1))
			if (nnj == 1) {
				legend("topright", names(spr_colors), col = spr_colors, lty = 1, bty = "n", ncol = 3, cex = 0.5)
			}
			# text(x = gen_m_d_f_a_i$start_rank_chr, y = 0, gen_m_d_f_a_i_tt, col = "gray5", cex = 0.6, pos = 4, offset = 0)
			# text(x = gen_m_d_f_a_j$start_rank_chr, y = 1, gen_m_d_f_a_j_tt, col = "gray5", cex = 0.6, pos = 4, offset = 0)
			# text(x = gen_m_d_f_a_i$start_rank_chr, y = 0.05, gen_m_d_f_a_i_t[,2], col = "gray15", cex = 0.3, pos = 4, offset = 0)
			# text(x = gen_m_d_f_a_j$start_rank_chr, y = 0.95, gen_m_d_f_a_j_t[,2], col = "gray15", cex = 0.3, pos = 4, offset = 0)

		}
		dev.off()

	}
	
}