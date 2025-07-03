# libraries
library("rphast")
library("ape")
library("GenomicRanges")
source("../scripts/helper.R")
source("../scripts/mta_downstream_functions.R")

# reference species
spi  = "Ocuarb"
# spi  = "Amil"
# spi  = "Spis"

# files
gix_fn = sprintf("../data/reference/%s_gDNA.fasta.fai", spi)
gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", spi)
out_fn = "results_wga_phastcons_plus/"
dir.create(out_fn, showWarnings = FALSE)

# get chrom lengths
gix = read.table(gix_fn)
gix_v = gix[,2]
names(gix_v) = gix[,1]

# # gff file with TEs
# tes_fn = sprintf("../data/reference/%s_EDTA.TEanno.original_scaffolds.gff", spi)
# tes_r = rtracklayer::import(tes_fn)
# tes_f = rphast::feat(seq = seqnames(tes_r), start = start(tes_r), end = end(tes_r), feature = tes_r$type , attribute = tes_r$ID)

## Data ##

# reload training
message(sprintf("phast | train based on %s | load training data...", spi))
pha_mod_neutral_pc = readRDS(sprintf("%s/model.phastcons.ref%s.rds", out_fn, spi))

# create genomic tiles
win_lenghth = 200
win_step = 100
gix_r = GenomicRanges::GRanges(seqnames = gix[,1], IRanges::IRanges(start = 1, end = gix[,2]))
giw_r = GenomicRanges::slidingWindows(gix_r, win_lenghth, step = win_step)
names(giw_r) = unique(gix[,1])


## Scoring ##

# first, init empty ranges
phr_r = GenomicRanges::GRanges()
pha_r = GenomicRanges::GRanges()
php_r = GenomicRanges::GRanges()
ppp_r = GenomicRanges::GRanges()
ppp_r_giw = GenomicRanges::GRanges()
# ppp_r_tes = GenomicRanges::GRanges()

# list of chromosomes
msa_fo = "results_wga/results_alignments_plus/"
msa_fl = list.files(msa_fo, pattern = ".maf$", full.names = TRUE)
msa_fl = msa_fl [ grepl(sprintf("cactus.%s", spi), msa_fl) ]

for (msa_fi in msa_fl) {

	# log
	chr_i = basename(msa_fi)
	chr_i = gsub(".maf$","", chr_i)
	chr_i = gsub(sprintf("^cactus\\.%s\\.", spi),"",chr_i)
	message(sprintf("phast | %s | chr %s\t| load...", spi, chr_i))

	# read alignments for chromosome i
	msa_i = rphast::read.msa(msa_fi)

	# phastcons
	message(sprintf("phast | %s | chr %s\t| phastCons...", spi, chr_i))
	pha_pc = rphast::phastCons(msa = msa_i, mod = pha_mod_neutral_pc$tree.models, viterbi = TRUE, quiet = FALSE)

	# phyloP
	message(sprintf("phast | %s | chr %s\t| phyloP...", spi, chr_i))
	pha_pp = rphast::phyloP(msa = msa_i, mod = pha_mod_neutral_pc$tree.models$noncons.mod, method = "LRT", mode = "CONACC")

	# phyloP on windows
	# select TEs from first chromosome
	giw_f_i = rphast::feat(seq = seqnames(giw_r[[chr_i]]), start = start(giw_r[[chr_i]]), end = end(giw_r[[chr_i]]), feature = "window")
	giw_f_i = giw_f_i [ order(giw_f_i$start, giw_f_i$end), ]
	if (length(giw_f_i$seqname) > 0 ) giw_f_i$seqname = spi
	message(sprintf("phast | %s | chr %s\t| phyloP on windows...", spi, chr_i))
	pha_pp_giw = rphast::phyloP(msa = msa_i, mod = pha_mod_neutral_pc$tree.models$noncons.mod, method = "LRT", mode = "CONACC", features = giw_f_i)

	# # phyloP on TEs
	# # select TEs from first chromosome
	# tes_f_i = tes_f [ tes_f$seqname == chr_i, ]
	# tes_f_i = tes_f_i [ order(tes_f_i$start, tes_f_i$end), ]
	# if (length(tes_f_i$seqname) > 0 ) tes_f_i$seqname = spi
	# message(sprintf("phast | %s | chr %s\t| phyloP on TE features...", spi, chr_i))
	# pha_pp_tes = rphast::phyloP(msa = msa_i, mod = pha_mod_neutral_pc$tree.models$noncons.mod, method = "LRT", mode = "CONACC", features = tes_f_i)
	# pha_pp_tes$name = tes_f_i$attribute

	# save as granges
	if (nrow(pha_pc$most.conserved) > 0) {
		phr_r = c(phr_r, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pc$most.conserved$start, end = pha_pc$most.conserved$end), strand = pha_pc$most.conserved$strand, score = pha_pc$most.conserved$score))
	} 
	if (nrow(pha_pc$post.prob.wig) > 0) {
		pha_r = c(pha_r, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pc$post.prob.wig$coord, width = 1), phastcons_pp = pha_pc$post.prob.wig$post.prob, phylop_score = pha_pp$score))
		php_r = c(php_r, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pc$post.prob.wig$coord, width = 1), score = pha_pc$post.prob.wig$post.prob))
		ppp_r = c(ppp_r, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pp$coord, width = 1), score = pha_pp$score))
		ppp_r_giw = c(ppp_r_giw, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pp_giw$start, end = pha_pp_giw$end), score = pha_pp_giw$score))
		# ppp_r_tes = c(ppp_r_tes, GenomicRanges::GRanges(seqnames = chr_i, IRanges::IRanges(start = pha_pp_tes$start, end = pha_pp_tes$end), score = pha_pp_tes$score, name = pha_pp_tes$name, class = pha_pp_tes$feature))
	}

}

# bed
message(sprintf("phast | %s | save bed...", spi))
mta_granges_to_bed(gr = phr_r, bed_fn = sprintf("%s/phastcons.%s.phastcons_top_regions.bed", out_fn, spi))
mta_granges_to_bed(gr = ppp_r_giw, bed_fn = sprintf("%s/phastcons.%s.phylop_windows.bed", out_fn, spi))
# mta_granges_to_bed(gr = ppp_r_tes, bed_fn = sprintf("%s/phastcons.%s.phylop_transposons.bed", out_fn, spi))

# bigwigs (requires seqlengths)
message(sprintf("phast | %s | save bigwig...", spi))
seqlengths(php_r) = gix_v [ levels(seqnames(php_r)) ]
seqlengths(ppp_r) = gix_v [ levels(seqnames(ppp_r)) ]
seqlengths(ppp_r_giw) = gix_v [ levels(seqnames(ppp_r_giw)) ]
rtracklayer::export.bw(php_r, con = sprintf("%s/phastcons.%s.phastcons_scores.bw",  out_fn, spi))
rtracklayer::export.bw(ppp_r, con = sprintf("%s/phastcons.%s.phylop_scores.bw", out_fn, spi))
ppp_r_giw_to_export = ppp_r_giw
width(ppp_r_giw_to_export) = win_step
rtracklayer::export.bw(ppp_r_giw_to_export, con = sprintf("%s/phastcons.%s.phylop_windows.bw", out_fn, spi))

pdf(sprintf("%s/phastcons.%s.phylop_windows.pdf", out_fn, spi), width = 7, height = 5)
hist(ppp_r_giw$score[ppp_r_giw$score!=0], breaks = 60, xlab = "PhyloP score", main = "PhyloP regions", col = "lightblue2", las = 1)
title(sub=sprintf("%i slow- and %i fast-evolving regions (P<1e-6), %i aligned out of %i total", sum(ppp_r_giw$score >= 6), sum(ppp_r_giw$score <= -6 ), sum(ppp_r_giw$score != 0), length(ppp_r_giw) ))
abline(v=c(-6,0,6), lty = 2, col = "red")
dev.off()

# pdf(sprintf("%s/phastcons.%s.phylop_transposons.pdf", out_fn, spi), width = 7, height = 5)
# hist(ppp_r_tes$score[ppp_r_tes$score!=0], breaks = 60, xlab = "PhyloP score", main = "PhyloP regions", col = "lightblue2", las = 1)
# title(sub=sprintf("%i slow- and %i fast-evolving regions (P<1e-6), %i aligned out of %i total", sum(ppp_r_tes$score >= 6), sum(ppp_r_tes$score <= -6 ), sum(ppp_r_tes$score != 0), length(ppp_r_tes) ))
# abline(v=c(-6,0,6), lty = 2, col = "red")
# dev.off()

message("all done!")
