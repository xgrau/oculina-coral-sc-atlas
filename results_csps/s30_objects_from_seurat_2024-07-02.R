# libraries
suppressMessages(source("../scripts/Seurat_functions.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(require("Seurat"))
suppressMessages(require("SeuratWrappers"))
graphics.off()


# list of datasets (first item is dataset id, second item is subsetting criteria, e.g. subset to a specific set of batches)
dat_list = list(
	Ocupat = list("Ocupat", NULL),
	Ocuarb = list("Ocuarb", NULL),
	Amil = list("Amil", NULL),
	Spin = list("Spin", NULL),
	Nvec = list("Nvec", NULL),
	Xesp = list("Xesp", NULL),
	opasym = list("Ocupat", list("batch_method", c("Ocupat02_OPUB","Ocupat04_Opat02U"))),
	opaapo = list("Ocupat", list("batch_method", c("Ocupat01_OPB","Ocupat03_Opat02B"))),
	oarsym = list("Ocuarb", list("batch_method", c("Ocuarb01_sym_SRR29367137"))),
	oarapo = list("Ocuarb", list("batch_method", c("Ocuarb02_apo_SRR29367138")))
)

# output
out_fn = sprintf("results_samap/data/")
dir.create(out_fn, showWarnings = FALSE)

for (nn in 1:length(dat_list)) {

	## Load Seurat ##

	spi = dat_list[[nn]][[1]]
	fil_field = dat_list[[nn]][[2]][[1]]
	fil_value = dat_list[[nn]][[2]][[2]]
	ddi = names(dat_list)[nn]

	# seurat object with batch info
	message(sprintf("metacell | %s | load Seurat for %s...", ddi, spi))
	seu = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.seurat_final.rds", spi, spi))
	
	# save cell type annotations
	message(sprintf("metacell | %s | save cell info...", ddi))
	cell_metadata = data.frame(cell = rownames(seu@meta.data), cell_type = seu@meta.data$cell_type)
	if (!is.na(dat_list[[nn]][2])) {
		cell_metadata = cell_metadata [ seu@meta.data[,fil_field] %in% fil_value, ]
	}
	write.table(cell_metadata, sprintf("%s/sam.%s.cell_to_cts.csv", out_fn, ddi), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)	
	
	# save genes
	message(sprintf("metacell | %s | save gene info...", ddi))
	gene_list = rownames(seu)
	gene_list = gsub("-","_", gene_list)
	write.table(gene_list, sprintf("%s/sam.%s.genes.txt", out_fn, ddi), sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)	
	
	# save counts matrix
	message(sprintf("metacell | %s | save mtx...", ddi))
	mm = seu@assays$RNA$counts [ rownames(seu) , cell_metadata$cell ]
	Matrix::writeMM(mm, sprintf("%s/sam.%s.counts.mtx", out_fn, ddi))
	
}

message("All done!")
