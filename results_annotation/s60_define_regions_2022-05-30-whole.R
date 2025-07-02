suppressMessages(library("rtracklayer"))
suppressMessages(library("zoo"))
suppressMessages(library("Matrix"))

# out
out_fn = "results_macrosynteny_plus"

# define species lists
sps_list = read.table("data/species_list_synteny_blocks_plus.txt")[,1]

# bin length and step
# applies to genes belonging to homology groups (others are ignored)
bin_w = 200
bin_s = 50

# iterate over regions & find counts of constituent homolgy groups, using fixed length

# homology file
hom_d = read.table(sprintf("%s/mcl.all.out.txt", out_fn), header = FALSE, col.names = c("homology_group", "transcript"))
hom_d = hom_d [ !duplicated(hom_d$transcript), ]
hom_d$species = gsub("_.*","", hom_d$transcript)
rownames(hom_d) = hom_d$transcript

# # retain only homology groups in >1 species
# hh = table(hom_d$homology_group, hom_d$species)
# hh_f = hh [ apply(hh, 1, function(v) sum(v>1) > 1 ), ]
# hh_v = rownames(hh_f)
# hom_d = hom_d [ hom_d$homology_group %in% hh_v, ]

for (spi in sps_list) {

	message(sprintf("%s | load %s...", "all", spi))
	gen = rtracklayer::readGFFAsGRanges(sprintf("../data/reference/%s_long.annot.gtf", spi))
	gen = gen [ gen$type == "transcript" ]
	mcols(gen)$homology_group = hom_d [ gen$transcript_id , "homology_group" ]
	list_hg = sort(unique(mcols(gen)$homology_group))
	list_hg = list_hg [ !is.na(list_hg) ]
	
	# remove genes without homology group (useless for cross-species comparisons)
	gen = gen [ !is.na(mcols(gen)$homology_group) ]
	
	# list of chromosomes
	chr_d = data.frame(chr = seqnames(gen), max_end = end(gen))
	chr_d = chr_d [ order(chr_d$chr, chr_d$max_end, decreasing = TRUE), ]
	chr_d = chr_d [ !duplicated(chr_d$chr) , ]
	chr_d = chr_d [ order(chr_d$max_end, decreasing = TRUE), ]
	
	# ignore small chromosomes
	# long_chromosomes = stringr::str_sort(names(which(sort(table(seqnames(gen)), decreasing=TRUE) > 500)), numeric = TRUE)
	# chr_d = chr_d [ chr_d$chr %in% long_chromosomes, ]
	
	# list of chrs
	chr_l = chr_d$chr 
	
	# init output
	met_d = data.frame()
	mat_d = data.frame()
	skip_c = 0
	skip_g = 0
	for (chr_i in chr_l) {
		
		gen_i = gen [ seqnames(gen) == chr_i ]
		hg_valid = length(table(mcols(gen_i)$homology_group)) > 0
		
		if (length(gen_i) >= bin_w & hg_valid)  {
			
			message(sprintf("%s | %s %s, n=%i genes", "all", spi, chr_i, length(gen_i)))
			mav_i = zoo::rollapply(mcols(gen_i)$homology_group,   width = bin_w, by = bin_s, function(v) v)
			mat_i = t(apply(mav_i, 1, function(r) { 
				r = factor(r, levels = list_hg) 
				table(r)
			} ))
			
			# bin names
			ixv_i = sprintf("%s:%s:%s%i", spi, chr_i, "bin", 1:nrow(mav_i))
			rownames(mat_i) = ixv_i

			# metadata for each bin
			# coordinates
			cmi_i = zoo::rollapply(start(gen_i), width = bin_w, by = bin_s, function(v) min(v))
			cma_i = zoo::rollapply(end(gen_i),   width = bin_w, by = bin_s, function(v) max(v))
			# list genes
			gnl_i = zoo::rollapply(mcols(gen_i)$transcript_id,   width = bin_w, by = bin_s, function(v) v)
			gnl_i = apply(gnl_i, 1, function(r) { paste(r, collapse = ",") })
			# list homology groups
			mvv_i = apply(mav_i, 1, function(r) { paste(r, collapse = ",") })
			# store metadata
			met_i = data.frame(
				chromosome = chr_i,
				start = cmi_i,
				end = cma_i,
				genes = gnl_i,
				homology_groups = mvv_i
			)
			rownames(met_i) = ixv_i
			
			# concatenate
			met_d = rbind(met_d, met_i)
			mat_d = rbind(mat_d, mat_i)
		
		} else {
			
			skip_c = skip_c + 1
			skip_g = skip_g + length(gen_i)
			
		}
		
	}
	
	# log
	message(sprintf("%s | %s %i chromosomes with n=%i genes have been skipped (too few, too discontiguous)", "all", spi, skip_c, skip_g))

	# save output	
	mat_d_s = Matrix::Matrix(as.matrix(mat_d), sparse = TRUE)
	Matrix::writeMM(mat_d_s, file = sprintf("%s/bins/bin.all-%s.mat.csv", out_fn, spi))
	write.table(met_d,              sprintf("%s/bins/bin.all-%s.dat.csv", out_fn, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
	write.table(colnames(mat_d),    sprintf("%s/bins/bin.all-%s.hgs.txt", out_fn, spi), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
	
}


message("all done!")
