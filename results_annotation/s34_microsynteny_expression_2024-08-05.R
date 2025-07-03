# load libs
suppressMessages(require("ape"))
suppressMessages(require("igraph"))
suppressMessages(require("rtracklayer"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
graphics.off()

## Define input ##

# parse info (will restrict pairs to this distance)
max_dis = 1

# species
sps_list = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp")

# input
dat_fn = "../data/orthology_Metazoa_plus/orthogroup_conservation.csv"
daa_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"
gff_fo = "../data/reference/" # where can I find gff files?

# output
out_fn = "results_microsynteny_plus/"
dir.create(out_fn, showWarnings = FALSE)

# bins
bin_width = 1
promoter_len = 1000

# read gene classification, taxonomy, etc
message("microsynteny | read orthology")
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
daa = read.table(daa_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ort_names_v = dic_from_vecs(dat$orthogroup, dat$orthogroup_name)
ort_genes_v = dic_from_vecs(dat$gene, dat$orthogroup)
ort_genon_v = dic_from_vecs(dat$gene, dat$orthogroup_name)
ora_names_v = dic_from_vecs(daa$orthogroup, daa$orthogroup_name)
ora_genes_v = dic_from_vecs(daa$gene, daa$orthogroup)
ora_genon_v = dic_from_vecs(daa$gene, daa$orthogroup_name)
oga_gtv = dic_from_vecs(daa$transcript, daa$gene)

# heatmap col
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))

# function to find gene pairs along chr
retrieve_pair_info = function(gff, dis) {

	# find neighbouring gene at distance dis
	gff$gene_b   = ifelse(gff$seqid == data.table::shift(gff$seqid, n=-dis), data.table::shift(gff$gene_id, n=-dis), NA)

	# find its coordinates
	gff$start_b  = ifelse(gff$seqid == data.table::shift(gff$seqid, n=-dis), data.table::shift(gff$start, n=-dis), NA)
	gff$end_b    = ifelse(gff$seqid == data.table::shift(gff$seqid, n=-dis), data.table::shift(gff$end, n=-dis), NA)
	gff$strand_b = ifelse(gff$seqid == data.table::shift(gff$seqid, n=-dis), data.table::shift(gff$strand, n=-dis), NA)

	# pair info: arrangement
	gff$orientation = NA
	gff$orientation [ gff$strand == "-" & gff$strand_b == "+" ] = "head-to-head"
	gff$orientation [ gff$strand == "+" & gff$strand_b == "-" ] = "tail-to-tail"
	gff$orientation [ gff$strand == gff$strand_b ] = "collinear"

	# pair distance: num genes
	gff$gene_distance = dis

	# return   
	return(gff)

}

# function to create matrix from table
table_to_matrix = function(table) {
    mat = matrix(table, nrow = nrow(table))
    rownames(mat) = rownames(table)
    colnames(mat) = colnames(table)
    return(mat)
}

# loop over species
fgbg_d_t = data.frame()
for (spi in sps_list) {
	
	
	if (spi == "Spis") {
		spi_dd = "Spin"
	} else {
		spi_dd = spi
	}

	
	# loop other sps
	for (spj in sps_list) {


		if (spj == "Spis") {
			spj_dd = "Spin"
		} else {
			spj_dd = spj
		}
		
		if (spi != spj) {
			
			# load GFF
			message(sprintf("synteny | %s - %s load GFF...", spi, spj))
			gi = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo, spi))
			gj = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo, spj))
			gig = gi[gi$type == "transcript",]
			gjg = gj[gj$type == "transcript",]
			gig = gig[order(gig$seqid, gig$start),]
			gjg = gjg[order(gjg$seqid, gjg$start),]

			# all possible pairs, spi
			gig_all = data.frame()
			for (dii in 1:max_dis) {
				gig_i = retrieve_pair_info(gig, dis = dii)
				gig_all = rbind(gig_all, gig_i)
			}
			gig_all$gene_pair = gsub("_", "-", paste(gig_all$gene_id, gig_all$gene_b))
			gig_all$gene_id = gsub("_", "-", gig_all$gene_id)
			gig_all$gene_b = gsub("_", "-", gig_all$gene_b)
			
			# all possible pairs, spj
			gjg_all = data.frame()
			for (dii in 1:max_dis) {
				gjg_i = retrieve_pair_info(gjg, dis = dii)
				gjg_all = rbind(gjg_all, gjg_i)
			}
			gjg_all$gene_pair = gsub("_", "-", paste(gjg_all$gene_id, gjg_all$gene_b))
			gjg_all$gene_id = gsub("_", "-", gjg_all$gene_id)
			gjg_all$gene_b = gsub("_", "-", gjg_all$gene_b)
			
			# load pairs
			message(sprintf("synteny | %s - %s synteny pairs...", spi, spj))
			msp = read.table(sprintf("%s/pairs_%s-%s.csv", out_fn, spi, spj), header = TRUE, sep = "\t")
			msp = msp [ msp$gene_distance_i <= max_dis , ]
			msp$gene_a_i = gsub("_", "-", oga_gtv [ msp$transcript_a_i ])
			msp$gene_b_i = gsub("_", "-", oga_gtv [ msp$transcript_b_i ])
			msp$gene_a_j = gsub("_", "-", oga_gtv [ msp$transcript_a_j ])
			msp$gene_b_j = gsub("_", "-", oga_gtv [ msp$transcript_b_j ])

			# find background (pairs not in synteny files)
			message(sprintf("synteny | %s - %s background pairs...", spi, spj))
			msp$pair_i = paste(msp$gene_a_i, msp$gene_b_i)
			msp$pair_j = paste(msp$gene_a_j, msp$gene_b_j)
			gig_bg = gig_all [ !gig_all$gene_pair %in% msp$pair_i, ]
			gjg_bg = gjg_all [ !gjg_all$gene_pair %in% msp$pair_j, ]

			# orientation tag (fix)
			msp$is_head_to_head_i [ msp$strand_a_i == "+" & msp$strand_a_i == "-" ] = TRUE
			msp$is_head_to_head_j [ msp$strand_a_i != msp$strand_b_i ] = TRUE
			msp$orientation_i = NA
			msp$orientation_j = NA
			msp$orientation_i [ msp$strand_a_i == "-" & msp$strand_b_i == "+" ] = "head-to-head"
			msp$orientation_j [ msp$strand_a_j == "-" & msp$strand_b_j == "+" ] = "head-to-head"
			msp$orientation_i [ msp$strand_a_i == "+" & msp$strand_b_i == "-" ] = "tail-to-tail"
			msp$orientation_j [ msp$strand_a_j == "+" & msp$strand_b_j == "-" ] = "tail-to-tail"
			msp$orientation_i [ msp$strand_a_i == msp$strand_b_i ] = "collinear"
			msp$orientation_j [ msp$strand_a_j == msp$strand_b_j ] = "collinear"
			
			# load expression
			message(sprintf("synteny | %s - %s load mcfp...", spi, spj))
			mcfp_i = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi_dd, spi_dd))
			mcfp_j = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spj_dd, spj_dd))
			# mcuc_i = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_umicount.rds", spi_dd, spi_dd))
			# mcuc_j = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_umicount.rds", spj_dd, spj_dd))

			message(sprintf("synteny | %s - %s evaluate n=%i gene pairs at max d=%i...", spi, spj, nrow(msp), max_dis))
			set.seed(1)
			fgbg_l = sapply(1:nrow(msp), function(nn) {
				
				# foreground
				gai = msp$gene_a_i[nn]
				gbi = msp$gene_b_i[nn]
				gaj = msp$gene_a_j[nn]
				gbj = msp$gene_b_j[nn]
				
				# distance and orientation in foreground
				ddi = msp$gene_distance_i[nn]
				ddj = msp$gene_distance_j[nn]
				ori = msp$orientation_i[nn]
				orj = msp$orientation_j[nn]
				
				# distance- and orientation-matched background
				bg_i = gig_bg[gig_bg$gene_distance == ddi & gig_bg$orientation == ori, c("gene_id","gene_b")]
				# bg_j = gjg_bg[gjg_bg$gene_distance == ddi & gjg_bg$orientation == ori, c("gene_id","gene_b")]
				bg_i = bg_i [ !is.na(bg_i[,1]), ]
				bg_i = bg_i [ !is.na(bg_i[,2]), ]
				# bg_j = bg_j [ !is.na(bg_j[,1]), ]
				# bg_j = bg_j [ !is.na(bg_j[,2]), ]
				bg_i_ix = sample(1:nrow(bg_i), 1)
				# bg_j_ix = sample(1:nrow(bg_j), 1)
				bg_gai = bg_i [ bg_i_ix, "gene_id" ]
				bg_gbi = bg_i [ bg_i_ix, "gene_b" ]
				# bg_gaj = bg_j [ bg_j_ix, "gene_id" ]
				# bg_gbj = bg_j [ bg_j_ix, "gene_b" ]
				
				genes_i = c(gai, gbi, bg_gai, bg_gbi)
				# genes_j = c(gaj, gbj, bg_gaj, bg_gbj)
				
				if (all(genes_i %in% rownames(mcfp_i))) {
					
					# foreground
					fg_eai = mcfp_i[gai,]
					fg_ebi = mcfp_i[gbi,]
					# fg_eaj = mcfp_j[gaj,]
					# fg_ebj = mcfp_j[gbj,]
					fg_ggci = stats::cor(fg_eai, fg_ebi, method = "pearson")
					# fg_ggcj = stats::cor(fg_eaj, fg_ebj, method = "pearson")
					
					# background
					bg_eai = mcfp_i[bg_gai,]
					bg_ebi = mcfp_i[bg_gbi,]
					# bg_eaj = mcfp_j[bg_gaj,]
					# bg_ebj = mcfp_j[bg_gbj,]
					bg_ggci = stats::cor(bg_eai, bg_ebi, method = "pearson")
					# bf_ggcj = stats::cor(bg_eaj, bg_ebj, method = "pearson")
					
					fgbgl = list(fg = fg_ggci, bg = bg_ggci, fg_ori = ori, fg_orj = orj, fg_a = gai, fg_b = gbi, bg_a = bg_gai, bg_b = bg_gbi)
					
				} else {

					fgbgl = list(fg = NA, bg = NA, fg_ori = ori, fg_orj = orj, fg_a = gai, fg_b = gbi, bg_a = bg_gai, bg_b = bg_gbi)
					
				}
				
			})
			
			# turn to dataframe
			fgbg_d = data.frame(
				gene_fg_a = unlist(fgbg_l["fg_a",]),
				gene_fg_b = unlist(fgbg_l["fg_b",]),
				gene_bg_a = unlist(fgbg_l["bg_a",]),
				gene_bg_b = unlist(fgbg_l["bg_b",]),
				correlation_fg = unlist(fgbg_l["fg",]),
				correlation_bg = unlist(fgbg_l["bg",]),
				orientation_i = unlist(fgbg_l["fg_ori",]),
				orientation_j = unlist(fgbg_l["fg_orj",])
			)
			fgbg_d = fgbg_d [ !is.na(fgbg_d$correlation_fg) & !is.na(fgbg_d$correlation_bg), ]
			fgbg_d$spi = spi
			fgbg_d$spj = spj
			fgbg_d_t = rbind(fgbg_d_t, fgbg_d)
			# fgbg_ll = list(fgbg_l[1,])
			
			# plot similarity of expression
			pdf(sprintf("%s/expcon.%s-%s.boxplots.pdf", out_fn, spi, spj), width = 6, height = 4)
			layout(matrix(1:3, nrow = 1))
			
			# all pairs
			fgbg_d_f = fgbg_d
			ks_t = ks.test(fgbg_d_f$correlation_fg, fgbg_d_f$correlation_bg)
			boxplot(list(fg = fgbg_d_f$correlation_fg, bg = fgbg_d_f$correlation_bg), col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
			title(main = sprintf("correlation of %s gene pairs\nsyntenic in %s with matched bgs", spi, spj), cex.main = 0.8, font.main = 1)
			title(sub = sprintf("%s\nn=%i pairs\np = %.2E | D = %.2f", ks_t$method, nrow(fgbg_d), ks_t$p.value, ks_t$statistic), cex.sub = 0.6)
			text(x = 1, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_fg)), col = "darkblue", cex = 0.6)
			text(x = 2, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_bg)), col = "darkblue", cex = 0.6)
			abline(h = 0, lty = 2, col = "darkred")
			
			# pairs that are h2h in both species
			fgbg_d_f = fgbg_d [ fgbg_d$orientation_i == "head-to-head" & fgbg_d$orientation_j == "head-to-head", ]
			ks_t = ks.test(fgbg_d_f$correlation_fg, fgbg_d_f$correlation_bg)
			boxplot(list(fg = fgbg_d_f$correlation_fg, bg = fgbg_d_f$correlation_bg), col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
			title(main = sprintf("correlation of %s h-to-h pairs\nsyntenic in %s with matched bgs", spi, spj), cex.main = 0.8, font.main = 1)
			title(sub = sprintf("%s\nn=%i pairs\np = %.2E | D = %.2f", ks_t$method, nrow(fgbg_d), ks_t$p.value, ks_t$statistic), cex.sub = 0.6)
			text(x = 1, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_fg)), col = "darkblue", cex = 0.6)
			text(x = 2, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_bg)), col = "darkblue", cex = 0.6)
			abline(h = 0, lty = 2, col = "darkred")
			
			# pairs that have not been inverted in any species
			fgbg_d_f = fgbg_d [ fgbg_d$orientation_i == fgbg_d$orientation_j , ]
			ks_t = ks.test(fgbg_d_f$correlation_fg, fgbg_d_f$correlation_bg)
			boxplot(list(fg = fgbg_d_f$correlation_fg, bg = fgbg_d_f$correlation_bg), col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
			title(main = sprintf("correlation of %s noninverted pairs\nsyntenic in %s, compare to matched bgs", spi, spj), cex.main = 0.8, font.main = 1)
			title(sub = sprintf("%s\nn=%i pairs\np = %.2E | D = %.2f", ks_t$method, nrow(fgbg_d), ks_t$p.value, ks_t$statistic), cex.sub = 0.6)
			text(x = 1, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_fg)), col = "darkblue", cex = 0.6)
			text(x = 2, y = 0, sprintf("%.3f", median(fgbg_d_f$correlation_bg)), col = "darkblue", cex = 0.6)
			abline(h = 0, lty = 2, col = "darkred")
			
			# close plot
			dev.off()
			
			
			# compare orientation conservation
			xtb = table(msp$orientation_i, msp$orientation_j)
			xtf = xtb / rowSums(xtb)
			
			# plot
			pdf(sprintf("%s/expcon.%s-%s.orientation.pdf", out_fn, spi, spj), width = 6, height = 4)

			# counts
			hm1 = pheatmap::pheatmap(
				table_to_matrix(xtb), 
				color = col_blue(20), 
				breaks = seq(1,800,length.out = 20 + 1), 
				cellwidth = 10, cellheight = 10, na_col = "grey", 
				cluster_cols = FALSE, cluster_rows = FALSE, 
				number_format = "%i", fontsize = 5, number_color = "orange",
				border_color = "white", 
				display_numbers = TRUE, 
				main = sprintf("orientation combinations\n(row=%s, col=%s)\nn=%i", spi, spj, sum(xtb)),
				main.font = 1)
			print(hm1)

			# fraction
			hm2 = pheatmap::pheatmap(
				table_to_matrix(xtf), 
				color = col_blue(20), 
				breaks = seq(0,1,length.out = 20 + 1), 
				cellwidth = 10, cellheight = 10, na_col = "grey", 
				cluster_cols = FALSE, cluster_rows = FALSE, 
				number_format = "%.2f", fontsize = 5, number_color = "orange",
				border_color = "white", 
				display_numbers = TRUE, 
				main = sprintf("orientation combinations\n(row=%s, col=%s)\nn=%i", spi, spj, sum(xtb)),
				main.font = 1)
			print(hm2)

			dev.off()
			
			
			# compare sequence similarity in gene vicinity?
			if (spi %in% c("Ocupat")) {
				
				# load
				message(sprintf("synteny | conservation vicinity | %s phastcons load...", spi))
				phc_r = rtracklayer::import.bw(sprintf("results_wga_phastcons_plus/phastcons.%s.phastcons_scores.bw", spi))	
				message(sprintf("synteny | conservation vicinity | %s gff load...", spi))
				gig_r = rtracklayer::import.gff(sprintf("%s/%s_long.annot.gtf", gff_fo, spi))	
				gig_r = gig_r [ gig_r$type == "transcript" ]
				gig_r = gig_r [ grepl("chr", as.character(seqnames(gig_r))) ]
				seqlevels(gig_r) = unique(as.character(seqnames(gig_r)))
				gig_r$gene_id = gsub("_","-", gig_r$gene_id)
				
				# seqlengths
				message(sprintf("synteny | conservation vicinity | %s seqlen load...", spi))
				gix_fn = sprintf("../data/reference/%s_gDNA.fasta.fai", spi)
				gix = read.table(gix_fn)
				gix_v = gix[,2]
				names(gix_v) = gix[,1]
				# gix_r = GenomicRanges::GRanges(seqnames = gix[,1], IRanges::IRanges(start = 1, end = gix[,2]))
				# giw_r = GenomicRanges::slidingWindows(gix_r, win_lenghth, step = win_step)
				# names(giw_r) = unique(gix[,1])
				seqlengths(gig_r) = gix_v [ as.character(seqlevels(gig_r)) ]

				# strand dic
				strand_v = dic_from_vecs(as.character(gig_r$gene_id), as.character(strand(gig_r)))

				# fg and bg (restrict to shared h2h)
				message(sprintf("synteny | conservation vicinity | %s set fg and bg ranges", spi))
				fgbg_d_f = fgbg_d [ fgbg_d$orientation_i == "head-to-head" & fgbg_d$orientation_j == "head-to-head", ] 
				gig_fg_r = gig_r [ gig_r$gene_id %in% c(fgbg_d_f$gene_fg_a, fgbg_d_f$gene_fg_b), ]
				gig_bg_r = gig_r [ gig_r$gene_id %in% c(fgbg_d_f$gene_bg_a, fgbg_d_f$gene_bg_b), ]
				
				# promoter regions
				message(sprintf("synteny | conservation vicinity | %s find promoters of fg and bg", spi))
				gig_fg_r_pro = GenomicRanges::promoters(gig_fg_r, upstream = promoter_len, downstream = 0)
				gig_bg_r_pro = GenomicRanges::promoters(gig_bg_r, upstream = promoter_len, downstream = 0)
				
				# sliding windows in fg
				message(sprintf("synteny | conservation vicinity | %s sliding windows, fg...", spi))
				gig_fg_r_pro_t = GenomicRanges::slidingWindows(gig_fg_r_pro, width = bin_width, step = bin_width)
				gig_fg_r_pro_t_u = unlist(gig_fg_r_pro_t)
				gig_fg_r_pro_t_u$gene_id = sapply(gig_fg_r_pro$gene_id, function(vv) { rep(vv, promoter_len) } )
				gig_fg_r_pro_t_u$bin     = sapply(gig_fg_r_pro$gene_id, function(vv) { 
					if (strand_v[vv] == "+") { 
						bb = (promoter_len/bin_width):1
					} else {
						bb = 1:(promoter_len/bin_width)
					}
				} )
				
				# sliding windows in bg
				message(sprintf("synteny | conservation vicinity | %s sliding windows, bg...", spi))
				gig_bg_r_pro_t = GenomicRanges::slidingWindows(gig_bg_r_pro, width = bin_width, step = bin_width)
				gig_bg_r_pro_t_u = unlist(gig_bg_r_pro_t)
				gig_bg_r_pro_t_u$gene_id = sapply(gig_bg_r_pro$gene_id, function(vv) { rep(vv, promoter_len) } )
				gig_bg_r_pro_t_u$bin     = sapply(gig_bg_r_pro$gene_id, function(vv) { 
					if (strand_v[vv] == "+") { 
						bb = (promoter_len/bin_width):1
					} else {
						bb = 1:(promoter_len/bin_width)
					}
				} )
				
				# coverage track
				message(sprintf("synteny | conservation vicinity | %s phastcons rle object...", spi))
				phc_r = phc_r [ as.character(seqnames(phc_r)) %in% seqlevels(gig_r) ]
				seqlevels(phc_r) = seqlevels(gig_r)
				seqlengths(phc_r) = seqlengths(gig_r)
				phc_i_r_rle = GenomicRanges::coverage(phc_r, weight = "score")
				# shared_chrs = intersect(as.character(seqlevels(gig_fg_r_pro_t_u)), names(phc_i_r_rle))
				# gig_fg_r_pro_t_u = gig_fg_r_pro_t_u [ as.character(seqnames(gig_fg_r_pro_t_u)) %in% shared_chrs ]
				# gig_bg_r_pro_t_u = gig_bg_r_pro_t_u [ as.character(seqnames(gig_bg_r_pro_t_u)) %in% shared_chrs ]
				# seqlevels(gig_fg_r_pro_t_u) = shared_chrs
				# seqlevels(gig_bg_r_pro_t_u) = shared_chrs

				message(sprintf("synteny | conservation vicinity | %s phastcons binned average...", spi))
            	gig_fg_r_pro_t_u = GenomicRanges::binnedAverage(gig_fg_r_pro_t_u, phc_i_r_rle, "phastcons")
            	gig_bg_r_pro_t_u = GenomicRanges::binnedAverage(gig_bg_r_pro_t_u, phc_i_r_rle, "phastcons")


				message(sprintf("synteny | conservation vicinity | %s phastcons aggregate, all...", spi))
				gig_fg_r_pro_d = data.frame(gig_fg_r_pro_t_u)
				gig_bg_r_pro_d = data.frame(gig_bg_r_pro_t_u)
				gig_fg_r_pro_da = aggregate(phastcons ~ bin, data = gig_fg_r_pro_d, function(vv) { mean(vv) })
				gig_bg_r_pro_da = aggregate(phastcons ~ bin, data = gig_bg_r_pro_d, function(vv) { mean(vv) })
				# gig_fg_r_pro_da = aggregate(phastcons ~ bin, data = gig_fg_r_pro_d, function(vv) { mean(vv) })
				# gig_bg_r_pro_da = aggregate(phastcons ~ bin, data = gig_bg_r_pro_d, function(vv) { mean(vv) })
				
				# open
				pdf(sprintf("%s/expcon.%s-%s.phastcons_promoters_new.pdf", out_fn, spi, spj), width = 4, height = 4)
				message(sprintf("synteny | conservation vicinity | %s phastcons aggregate, plot...", spi))

				plot(gig_fg_r_pro_da$bin, gig_fg_r_pro_da$phastcons, pch = 19, col = "purple", lwd = 1.5, ylim = c(0.2,0.6), type = "l", ylab = "phastcons", xlab = "bp from TSS", las = 2)
				lines(gig_bg_r_pro_da$bin, gig_bg_r_pro_da$phastcons, pch = 19, col = "gray50", lwd = 1.5)
				legend("topright", c("syntenic","nonsyntenic bg"), col = c("purple", "gray50"), pch = 19, bty = "n", cex = 0.7)
				title(main = sprintf("conservation in upstream of %s\nsyntenic in %s, compare to matched bgs", spi, spj), cex.main = 0.8, font.main = 1)
				title(sub = sprintf("fg = %i, bg = %i upstream regions", length(gig_fg_r_pro), length(gig_bg_r_pro)), cex.sub = 0.8)
				
				# stratify by correlation of foreground
				fgbg_d_f$correlation_fg_corbin = cut(fgbg_d_f$correlation_fg, breaks = c(-1,-0.2,0.2,1), include.lowest = TRUE)
				fgbg_d_corbin_v = dic_from_vecs(fgbg_d_f$gene_fg_a, fgbg_d_f$correlation_fg_corbin)
				col_v = viridisLite::viridis(nlevels(fgbg_d_corbin_v), begin = 0.7, end = 0.2)
				names(col_v) = levels(fgbg_d_corbin_v) 
				plot(gig_fg_r_pro_da$bin, gig_fg_r_pro_da$phastcons, pch = 19, col = "purple", lwd = 1.5, ylim = c(0.2,0.6), type = "l", ylab = "phastcons", xlab = "bp from TSS", las = 2)
				lines(gig_bg_r_pro_da$bin, gig_bg_r_pro_da$phastcons, pch = 19, col = "gray50", lwd = 1.5)
				for (bii in levels(fgbg_d_corbin_v)) {
					genes_bii = names(fgbg_d_corbin_v) [ fgbg_d_corbin_v == bii ] 
					gig_fg_r_pro_da = aggregate(phastcons ~ bin, data = gig_fg_r_pro_d[gig_fg_r_pro_d$gene_id %in% genes_bii,], function(vv) { mean(vv) })
					lines(gig_fg_r_pro_da$bin, gig_fg_r_pro_da$phastcons, pch = 19, col = col_v[bii], lwd = 1.5)
				}
				title(main = sprintf("conservation in upstream of %s\nsyntenic in %s, stratified by gene-gene correlation", spi, spj), cex.main = 0.8, font.main = 1)
				legend("topright", sprintf("%s,n=%i", names(col_v), table(fgbg_d_corbin_v)), col = col_v, pch = 19, bty = "n", cex = 0.7, title = "gene-gene coexpression, rho")
				dev.off()
				
			}
			

		}
	}
	
	# cross-species comparisons centered on each ref species
	fgbg_d_i = fgbg_d_t [ fgbg_d_t$spi == spi, ]

	# order of species to compare to
	if (spi == "Ocupat") {
		sps_list_ordered = c("Ocuarb","Spis","Amil","Nvec","Xesp")
		sps_list_orderix = c(1,2,3,4,5)
	} else if (spi == "Ocuarb") {
		sps_list_ordered = c("Ocupat","Spis","Amil","Nvec","Xesp")
		sps_list_orderix = c(1,2,3,4,5)
	} else if (spi == "Spis") {
		sps_list_ordered = c("Ocupat","Ocuarb","Amil","Nvec","Xesp")
		sps_list_orderix = c(1,1,2,3,4)
	} else if (spi == "Amil") {
		sps_list_ordered = c("Ocupat","Ocuarb","Spis","Nvec","Xesp")
		sps_list_orderix = c(1,1,1,2,3)
	} else if (spi == "Nvec") {
		sps_list_ordered = c("Ocupat","Ocuarb","Spis","Amil","Xesp")
		sps_list_orderix = c(1,1,1,1,2)
	} else if (spi == "Xesp") {
		sps_list_ordered = c("Ocupat","Ocuarb","Spis","Amil","Nvec")
		sps_list_orderix = c(1,1,1,1,1)
	}
	
	# apply order
	fgbg_d_i$spi = factor(fgbg_d_i$spi, levels = sps_list_ordered)
	fgbg_d_i$spj = factor(fgbg_d_i$spj, levels = sps_list_ordered)
	
	# open plot
	pdf(sprintf("%s/expcon.%s-all.pdf", out_fn, spi), width = 6, height = 8)
	layout(matrix(1:6, nrow = 2))

	# all pairs
	fgbg_d_i_ll = list()
	for (spj in sps_list_ordered [ !sps_list_ordered %in% spi ]) {
		fgbg_d_i_ll[[spj]] = fgbg_d_i$correlation_fg [ fgbg_d_i$spj == spj ]
	}
	boxplot(fgbg_d_i_ll, col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
	abline(h = 0, lty = 2, col = "darkred")
	text(x = 1:length(fgbg_d_i_ll), y = -1, sprintf("md=%.3f", unlist(lapply(fgbg_d_i_ll, median))), col = "darkblue", cex = 0.6, srt = 90)
	title(main = sprintf("correlation of %s pairs\nsyntenic in other species", spi), cex.main = 0.8, font.main = 1)
	# test with increasing distance
	fgbg_d_i_lld = data.frame(
		cor = unlist(fgbg_d_i_ll),
		spj = rep(names(fgbg_d_i_ll), lengths(fgbg_d_i_ll))
	)
	fgbg_d_i_lld$spj = factor(fgbg_d_i_lld$spj, levels = sps_list_ordered)
	fgbg_d_i_lld$spj_dist = sps_list_orderix [ fgbg_d_i_lld$spj ]
	sp_t = cor.test(fgbg_d_i_lld$cor, fgbg_d_i_lld$spj_dist, method = "spearman")
	title(sub = sprintf("spearman cor with phy dist\np=%.1E rho=%.4f", sp_t$p.value, sp_t$estimate), cex.sub = 0.6)
	# counts
	b = barplot(lengths(fgbg_d_i_ll), las = 2, col = "lightblue2", ylab = "Num pairs")
	text(b, y = 0, lengths(fgbg_d_i_ll), srt = 90, pos = 4, cex = 0.7, col = "darkblue")

	# head-to-head pairs
	fgbg_d_i_ll = list()
	for (spj in sps_list_ordered [ !sps_list_ordered %in% spi ]) {
		fgbg_d_i_ll[[spj]] = fgbg_d_i$correlation_fg [ fgbg_d_i$spj == spj & fgbg_d_i$orientation_i == "head-to-head" & fgbg_d_i$orientation_j == "head-to-head" ]
	}
	boxplot(fgbg_d_i_ll, col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
	abline(h = 0, lty = 2, col = "darkred")
	text(x = 1:length(fgbg_d_i_ll), y = -1, sprintf("md=%.3f", unlist(lapply(fgbg_d_i_ll, median))), col = "darkblue", cex = 0.6, srt = 90)
	title(main = sprintf("correlation of %s pairs\nsyntenic in other species", spi), cex.main = 0.8, font.main = 1)
	# test with increasing distance
	fgbg_d_i_lld = data.frame(
		cor = unlist(fgbg_d_i_ll),
		spj = rep(names(fgbg_d_i_ll), lengths(fgbg_d_i_ll))
	)
	fgbg_d_i_lld$spj = factor(fgbg_d_i_lld$spj, levels = sps_list_ordered)
	fgbg_d_i_lld$spj_dist = sps_list_orderix [ fgbg_d_i_lld$spj ]
	sp_t = cor.test(fgbg_d_i_lld$cor, fgbg_d_i_lld$spj_dist, method = "spearman")
	title(sub = sprintf("spearman cor with phy dist\np=%.1E rho=%.4f", sp_t$p.value, sp_t$estimate), cex.sub = 0.6)
	# counts
	b = barplot(lengths(fgbg_d_i_ll), las = 2, col = "lightblue2")
	text(b, y = 0, lengths(fgbg_d_i_ll), srt = 90, pos = 4, cex = 0.7, col = "darkblue")

	# noninverted pairs
	fgbg_d_i_ll = list()
	for (spj in sps_list_ordered [ !sps_list_ordered %in% spi ]) {
		fgbg_d_i_ll[[spj]] = fgbg_d_i$correlation_fg [ fgbg_d_i$spj == spj & fgbg_d_i$orientation_i == fgbg_d_i$orientation_j ]
	}
	boxplot(fgbg_d_i_ll, col = "lightblue3", ylim = c(-1,1), las = 2, ylab = "Pearson")
	abline(h = 0, lty = 2, col = "darkred")
	text(x = 1:length(fgbg_d_i_ll), y = -1, sprintf("md=%.3f", unlist(lapply(fgbg_d_i_ll, median))), col = "darkblue", cex = 0.6, srt = 90)
	title(main = sprintf("correlation of %s noninverted pairs\nsyntenic in other species", spi), cex.main = 0.8, font.main = 1)
	# test with increasing distance
	fgbg_d_i_lld = data.frame(
		cor = unlist(fgbg_d_i_ll),
		spj = rep(names(fgbg_d_i_ll), lengths(fgbg_d_i_ll))
	)
	fgbg_d_i_lld$spj = factor(fgbg_d_i_lld$spj, levels = sps_list_ordered)
	fgbg_d_i_lld$spj_dist = sps_list_orderix [ fgbg_d_i_lld$spj ]
	sp_t = cor.test(fgbg_d_i_lld$cor, fgbg_d_i_lld$spj_dist, method = "spearman")
	title(sub = sprintf("spearman cor with phy dist\np=%.1E rho=%.4f", sp_t$p.value, sp_t$estimate), cex.sub = 0.6)
	# counts
	b = barplot(lengths(fgbg_d_i_ll), las = 2, col = "lightblue2")
	text(b, y = 0, lengths(fgbg_d_i_ll), srt = 90, pos = 4, cex = 0.7, col = "darkblue")

	# close
	dev.off()
	
}