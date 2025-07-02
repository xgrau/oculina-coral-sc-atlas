library("rtracklayer")
library("zoo")
library("Matrix")
library("viridis")
source("../scripts/helper.R")
source("../scripts/Cross_species_functions.R")

# out
out_fn = "results_macrosynteny_plus"

# define species lists
sps_list = read.table("data/species_list_synteny_blocks_plus.txt")[,1]

# homology groups
hom_d = read.table(sprintf("%s/mcl.all.out.txt", out_fn), header = FALSE, col.names = c("homology_group", "transcript"))
hom_d$species = gsub("_.*","", hom_d$transcript)
# hom_d$species [ hom_d$species == "Hvul" ] = "Hvul_v3"
hom_d = hom_d [ hom_d$species %in% sps_list , ]
hom_d = hom_d [ !duplicated(hom_d$transcript), ]
rownames(hom_d) = hom_d$transcript

# color palettes
categorical_colors      = colorspace::darken(c("magenta4","violet","firebrick1","orange","khaki1","palegreen2","springgreen3","darkgreen"), 0.1)
categorical_colors_dark = colorspace::darken(c("deepskyblue","paleturquoise1","mediumblue","darkviolet","orchid1","thistle1"), 0.1)
color_palette           = colorRampPalette(categorical_colors)
color_palette_dark      = colorRampPalette(categorical_colors_dark)

dir.create(sprintf("%s/overlaps/", out_fn), showWarnings = FALSE)

# iterate over regions & find counts of constituent homolgy groups, using fixed length
for (spr in sps_list) {

	message(sprintf("%s | load reference %s...", spr, spr))
	mat_r = Matrix::readMM(file = sprintf("%s/bins/bin.all-%s.mat.csv", out_fn, spr))
	met_r = read.table(sprintf("%s/bins/bin.all-%s.dat.csv", out_fn, spr))
	hgs_r = read.table(sprintf("%s/bins/bin.all-%s.hgs.txt", out_fn, spr))[,1]
	colnames(mat_r) = hgs_r
	rownames(mat_r) = rownames(met_r)
	
	# color palette
	spr_colors = color_palette(n=length(unique(met_r$chromosome)))
	names(spr_colors) = unique(met_r$chromosome)
	
	pdf(sprintf("%s/overlaps/ovs.%s.mat.pdf", out_fn, spr), height = 12, width = 12)
	for (spi in sps_list [ sps_list != spr ]) {

		message(sprintf("%s | %s load query", spr, spi))
		mat_i = Matrix::readMM(file = sprintf("%s/bins/bin.all-%s.mat.csv", out_fn, spi))
		met_i = read.table(sprintf("%s/bins/bin.all-%s.dat.csv", out_fn, spi))
		hgs_i = read.table(sprintf("%s/bins/bin.all-%s.hgs.txt", out_fn, spi))[,1]
		colnames(mat_i) = hgs_i
		rownames(mat_i) = rownames(met_i)
		
		# filter out HGs that have undergone expansions in either species (noise)
		hom_i = hom_d [ hom_d$species %in% c(spi, spr), ]
		hom_i_m = as.matrix(table(hom_i$homology_group, hom_i$species))
		# ratio between number of copies in each species (highest to second highest)
		hom_i_ratio = apply(hom_i_m, 1, function(h) { 
			rs = sort(h) 
			ra = rs[2] / rs[1]
		})
		keep_hgs = names(which(hom_i_ratio < 4))
		
		# subset to comparable hgs
		joint_hgs = intersect(intersect(hgs_i, hgs_r), keep_hgs)
		mat_r_f = t(mat_r[, joint_hgs])
		mat_i_f = t(mat_i[, joint_hgs])
		mat_r_f [ mat_r_f > 1 ] = 1
		mat_i_f [ mat_i_f > 1 ] = 1
		
		# jaccard similarity
		message(sprintf("%s | %s correlation", spr, spi))
		mat_c = jaccard(as.matrix(mat_r_f), as.matrix(mat_i_f))
		mat_c_f = mat_c
		# mat_c_f = mat_c_f [ apply(mat_c_f,1,max) > 0.1 , ]
		# mat_c_f = mat_c_f [ , apply(mat_c_f,2,max) > 0.1 ]
		mat_o_f = mat_c_f
		
		# are there enough datapoints?
		do_heatmap = !is.null(dim(mat_o_f))
		if (do_heatmap) {
			do_heatmap = ncol(mat_o_f) > 1 & nrow(mat_o_f) > 1
		}

		if ( do_heatmap ) {
		
			message(sprintf("%s | %s plot overlaps", spr, spi))
		
			# color strings
			chrs_row = stringr::str_split(rownames(mat_o_f), ":", simplify = TRUE)[,2]
			chrs_col = stringr::str_split(colnames(mat_o_f), ":", simplify = TRUE)[,2]
			chrs_row = factor(chrs_row, levels = unique(stringr::str_split(rownames(mat_r), ":", simplify = TRUE)[,2]))
			chrs_col = factor(chrs_col, levels = unique(stringr::str_split(rownames(mat_i), ":", simplify = TRUE)[,2]))
			spi_colors = color_palette_dark(n=nlevels(chrs_col))
			colors_row = spr_colors [ chrs_row ]
			colors_col = spi_colors [ chrs_col ]
			names(colors_row) = chrs_row
			names(colors_col) = chrs_col
			
			hm = plot_complex_heatmap(
				mat_o_f, 
				name = "jac",
				color_mat = c("gray95","#d6e72e","#6fb600","#003f4d"),
				color_min = 0.05, 
				color_max = 0.25 ,
				categories_row = chrs_row, 
				categories_col = chrs_col,
				colors_row = colors_row,
				colors_col = colors_col,
				cluster_row = FALSE,
				cluster_col = FALSE,
				name_row_show = TRUE,
				name_col_show = TRUE,
				separate_row = TRUE,
				separate_col = TRUE,
				show_legend_row = TRUE,
				show_legend_col = TRUE,
				title_col = sprintf("%s gene blocks", spi), title_row = sprintf("%s gene blocks", spr)
			)
			hm@row_names_param$anno@width     = unit(1.5, "cm")
			hm@column_names_param$anno@height = unit(1.5, "cm")
			hm@heatmap_param$width  = unit(18, "cm")
			hm@heatmap_param$height = unit(18, "cm")
			print(hm)
		
		} else {
			
			message(sprintf("%s | not enough overlaps with %s", spr, spi))
			
		}
	
	}
	dev.off()

}

message("all done!")
