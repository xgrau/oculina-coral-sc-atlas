# load libraries
library("scales")
library("GenomicRanges")
source("../scripts/helper.R")
graphics.off()

# input files
out_fn = "results_blobtools/"
dir.create(out_fn, showWarnings = FALSE)
gix_fn = "../data/reference/Ocupat_gDNA.all_scaffolds.fasta.fai"
gen_fn = "../results_scatlas/annotate_3p/reannotate_Ocupat_genes.gtf"
cov_fn = "../data/sequencing_nanopore_Oculina/covcalc.paf.gz"

# read genome info
message("read: genome")
gen_r = rtracklayer::import(gen_fn)
gen_r = gen_r [ gen_r$type == "gene" ]
gix_d = read.table(gix_fn, col.names = c("chromosome","length","OFFSET","LINEBASES","LINEWIDTH"))
gix_d_lengths = dic_from_vecs(gix_d$chromosome, gix_d$length)
seqlevels(gen_r) = names(gix_d_lengths)
seqlengths(gen_r) = gix_d_lengths [ seqlevels(gen_r) ]

# read coverage
message("read: seq coverage")
cov = data.table::fread(cov_fn, sep = "\t", select = 1:12, fill = TRUE)
colnames(cov) = c("query","qlen","qstart","qend","strand","target","tlen","tstart","tend","num_matches","alignment_block_len","mapq")
cov_d = as.data.frame(cov)
cov_df = cov_d [ cov_d$mapq >= 20, ]
cov_r = GenomicRanges::makeGRangesFromDataFrame(
  cov_df,
  seqnames.field = "target",
  start.field = "tstart",
  end.field = "tend",
  na.rm=TRUE
)
seqlevels(cov_r) = names(gix_d_lengths)
seqlengths(cov_r) = gix_d_lengths [ seqlevels(cov_r) ]

# create genomic bins
message("do: bin genome")
gix_r = GenomicRanges::makeGRangesFromDataFrame(
  gix_d,
  seqnames.field = "chromosome",
  start.field = "length",
  end.field = "length"
)
start(gix_r) = 1
bin_r = GenomicRanges::slidingWindows(gix_r, width = 10000, step = 5000)
seqlevels(bin_r) = names(gix_d_lengths)
seqlengths(bin_r) = gix_d_lengths [ seqlevels(bin_r) ]
bin_r = GenomicRanges::trim(bin_r)
bin_r = unlist(bin_r)

# coverage per gene
message("do: ranges coverage")
cov_r_per_base_coverage = GenomicRanges::coverage(cov_r)
gen_r = GenomicRanges::binnedAverage(gen_r, cov_r_per_base_coverage, "coverage")

# output table
gencov_d = data.frame(gene = mcols(gen_r)$transcript_id, coverage = as.integer(mcols(gen_r)$coverage), bases = width(gen_r))

# header
ff = file(sprintf("%s/blobgene.coverage.cov", out_fn))
writeLines(
	c(
		sprintf("## blobtools v1.0"),
		sprintf("## Total Reads = %i", nrow(cov_d)),
		sprintf("## Mapped Reads = %i", nrow(cov_df)),
		sprintf("## Unmapped Reads = 0"),
		sprintf("## Source(s) : nye.bam"),
		sprintf("# contig_id	read_cov	base_cov")
	), 
	ff)
close(ff)

# table
write.table(gencov_d, sprintf("%s/blobgene.coverage.cov", out_fn), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)

message("all done!")

