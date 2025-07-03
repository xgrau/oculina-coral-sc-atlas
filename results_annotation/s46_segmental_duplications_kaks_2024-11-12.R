# load libs
suppressMessages(require("ape"))
suppressMessages(require("igraph"))
suppressMessages(require("rtracklayer"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))
suppressMessages(source("../scripts/Downstream_functions.R"))
graphics.off()

## Scope of analysis ##

list_tests = list(
	"allsclerplusplus" = list(sps_list = c("Ocupat","Ocuarb","Gasp","Fspp","Spis","Pocdam","Amil","Adig","Gfas"), anc_list = c("Scleractinia","Hexacorallia","Anthozoa"))
)

sps_focus = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp")


## Define input ##

# input
daa_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"
phy_fn = "../data/species_tree.Anthozoa_plus.newick"

# output
out_fn = "results_segmental_duplications_plus/"
dir.create(out_fn, showWarnings = FALSE)


# read gene classification, taxonomy, etc
message("posteriors | read orthology")
daa = read.table(daa_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ora_names_v = dic_from_vecs(daa$orthogroup, daa$orthogroup_name)
ora_genes_v = dic_from_vecs(daa$gene, daa$orthogroup)
ora_genon_v = dic_from_vecs(daa$gene, daa$orthogroup_name)
oga_gtv = dic_from_vecs(daa$transcript, daa$gene)

# read species tree
message("posteriors | read tree")
phy = ape::read.tree(phy_fn)
sps_list = phy$tip.label
nod_list = c(phy$tip.label, phy$node.label)
root_label = phy$node.label[1]


for (nn in 1:length(list_tests)) {
	
	tid = names(list_tests)[nn]
	anc_list = list_tests[[nn]]$anc_list
	sps_list = list_tests[[nn]]$sps_list
	sps_focus_f = sps_focus [ sps_focus %in% sps_list ]
	message(sprintf("posteriors | %s | ka/ks analysis in %i species in %i ancestors", tid, length(sps_focus_f), length(anc_list)))

	# reload
	segments_d = read.table(sprintf("%s/summary_segmental_duplications.%s.csv", out_fn, tid), sep = "\t", header = TRUE)
	segments_d = segments_d [ segments_d$species %in% sps_focus_f, ]
	segments_d$species = factor(segments_d$species, levels = sps_focus_f)
	message(sprintf("posteriors | %s | ka/ks analysis of %i segment genes across all species (%i families)", tid, nrow(segments_d), length(unique(segments_d$orthogroup))))
	
	# create output folder for cds alignments
	dir.create(sprintf("%s/alignments_%s/", out_fn, tid), showWarnings = FALSE)
	
	for (spi in sps_focus_f) { 
		
		if (spi == "Spis") {
			spi_dd = "Spin"
		} else {
			spi_dd = spi
		}
		
		# log
		message(sprintf("posteriors | %s | %s: load data...", tid, spi))

		# load ct annot
		ctt = read.table(sprintf("../results_scatlas/results_metacell_%s_filt/annot.%s.leiden.csv", spi_dd, spi_dd), sep = "\t", header = TRUE, comment.char = "")
		ctt$cell_type = factor(ctt$cell_type, levels = unique(ctt$cell_type))
		ctt_lou_cts_v = dic_from_vecs(ctt$cluster, ctt$cell_type)
		ctt_lou_col_v = dic_from_vecs(ctt$cluster, ctt$color)
		ctt_cts_col_v = dic_from_vecs(ctt$cell_type, ctt$color)
	
		# metacell
		ctm = read.table(sprintf("../results_scatlas/results_metacell_%s_filt/annot.%s.mcs.csv", spi_dd, spi_dd), sep = "\t", header = TRUE, comment.char = "")
		ctm$metacell = factor(ctm$metacell, levels = unique(ctm$metacell))
		ctm_mcs_col_v = dic_from_vecs(ctm$metacell, ctm$color)
		ctm_mcs_cts_v = dic_from_vecs(ctm$metacell, ctm$metacell_annotation)
	
		# load expression
		ctfp = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.cts_fp.rds", spi_dd, spi_dd))
		ctfp_v = apply(ctfp, 1, function(vv) { zz = which.max(vv); sprintf("%s(%.2f)", names(zz), max(vv)) })
		ctfp_vn = apply(ctfp, 1, function(vv) { zz = which.max(vv); sprintf("%s", names(zz)) })
		# mcfp = readRDS(sprintf("../results_scatlas/results_metacell_%s_filt/dat.%s.expression.mcs_fp.rds", spi_dd, spi_dd))
		# mcfp_v = apply(mcfp, 1, function(vv) { zz = which.max(vv); sprintf("%s(%.2f)", ctm_mcs_cts_v[names(zz)], max(vv)) })

		# load CDS fasta
		cds = seqinr::read.fasta(sprintf("../data/reference/%s_long.cds.fasta", spi))
		names(cds) = gsub("_", "-", as.character(oga_gtv [ names(cds) ]))
		
		# get segments spi
		segments_i = segments_d [ segments_d$species == spi, ]
		message(sprintf("posteriors | %s | %s: ka/ks analysis of %i segment genes (%i families)", tid, spi, nrow(segments_i), length(unique(segments_i$orthogroup))))
		
		# set gene-level table
		dupgenes_i = segments_i [ c("gene","orthogroup","duplicate_set","gene_name") ]
		dupgenes_i$has_redundant_copy = NA
		dupgenes_i$most_similar_copy  = NA
		dupgenes_i$corr_similar_copy  = NA
		dupgenes_i$pid_similar_copy   = NA
		dupgenes_i$status             = NA
		rownames(dupgenes_i) = gsub("_","-", dupgenes_i$gene)

		# family-wise ka/ks distributions
		for (ogi in unique(dupgenes_i$orthogroup)) {
			
			genes_ogi = dupgenes_i$gene [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = dupgenes_i$gene_name [ dupgenes_i$orthogroup == ogi ]
			ognam_ogi = paste(unique(sort(ognam_ogi)), collapse = " / ")
			genes_ogi_v = gsub("_","-", genes_ogi)
			
			# get expression matrix
			if (all(genes_ogi_v %in% rownames(ctfp))) {
				gge = ctfp[genes_ogi_v,]
			} else {
				genes_ogi_vv = genes_ogi_v [ genes_ogi_v %in% rownames(ctfp) ]
				genes_ogi_vn = setdiff(genes_ogi_v, rownames(ctfp))
				ctfp_pp = ctfp
				for (ggin in genes_ogi_vn) {
					ctfp_pp = rbind(ctfp, rep(NA, ncol(ctfp)))
				}
				rownames(ctfp_pp) = c(rownames(ctfp), genes_ogi_vn)
				gge = ctfp_pp[genes_ogi_v,]
			}
			
			# gene-gene correlation
			ggc_i = cor(t(gge), method = "pearson", use = "everything")
			ggc_i [ is.na(ggc_i) ] = 0
			
			# gene-gene distance using ka/ks (omit, slow...)
			# compute codon alignments clustalw
			cds_i = cds [ genes_ogi_v ]
			seqinr::write.fasta(cds_i, names = names(cds_i), file.out = sprintf("%s/alignments_%s/%s.%s.cds.fasta", out_fn, tid, ogi, spi))
			
			
			if (!all(unlist(lapply(1:(length(cds_i)-1), function(nnn) all(cds_i[[nnn]] == cds_i[[nnn+1]]))))) {
				invisible(seqinr::reverse.align(
					nucl.file =  sprintf("%s/alignments_%s/%s.%s.cds.fasta", out_fn, tid, ogi, spi),
					out.file =   sprintf("%s/alignments_%s/%s.%s.cds.l.fasta", out_fn, tid, ogi, spi),
					align.prot = TRUE,
					clustal.path = "/home/xavi/mambaforge/bin/clustalw"
				))
				# read alignment
				ali_i = suppressMessages(seqinr::read.alignment(sprintf("%s/alignments_%s/%s.%s.cds.l.fasta", out_fn, tid, ogi, spi), format = "fasta"))
				ali_i_kaks = suppressMessages(seqinr::kaks(ali_i))
				if (length(ali_i_kaks) > 1) {
					ali_i_m_ka = as.matrix(ali_i_kaks$ka)
					ali_i_m_ks = as.matrix(ali_i_kaks$ks)
					ali_i_m_kaks = ali_i_m_ka / ali_i_m_ks
					write.table(ali_i_m_ka, sprintf("%s/alignments_%s/%s.%s.kaks.ka.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
					write.table(ali_i_m_ks, sprintf("%s/alignments_%s/%s.%s.kaks.ks.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
					write.table(ali_i_m_kaks, sprintf("%s/alignments_%s/%s.%s.kaks.kaks.tsv", out_fn, tid, ogi, spi), sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
					
					# similarity
					ali_i_m_ks_sim = 1 - ali_i_m_ks / 10
					rownames(ali_i_m_ks_sim) = rownames(ali_i_m_ks)
					colnames(ali_i_m_ks_sim) = colnames(ali_i_m_ks)
				} else {
					ali_i_m_ks_sim = ggc_i
					ali_i_m_ks_sim [ ali_i_m_ks_sim != 0 ] = 0
					ali_i_m_ks = ali_i_m_ks_sim
				}
			}
			
		}
		
	}
	
	
}