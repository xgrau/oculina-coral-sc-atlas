# load libs
require("ape")
require("scales")
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/gene-set-analysis.R"))


## Define input ##

# input
prob_thr = 0.5
dat_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"
ort_fn = "results_gene_family_evolution_anthozoa_plus/orthogroup_counts.posteriors.G2.csv"
phy_fn = "../data/species_tree.Anthozoa_plus.newick"
out_fn = "results_gene_family_evolution_anthozoa_plus/"
dol_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.possvm.dollo_summary.csv"

# read gene classification, taxonomy, etc
message("posteriors | read orthology")
dat = read.table(dat_fn, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
ort_names_v = dic_from_vecs(dat$orthogroup, dat$orthogroup_name)
ort_genes_v = dic_from_vecs(dat$gene, dat$orthogroup)

# read dollo reconstruction
message("posteriors | read dollo")
dol = read.table(dol_fn, sep = "\t", header = TRUE)
rownames(dol) = dol$orthogroup


# read species tree
message("posteriors | read tree")
phyl = ape::read.tree(phy_fn)
sps_list = phyl$tip.label
phyl$edge.length = rep(1, nrow(phyl$edge))
# prepare species dataframe:
# dataframe of edges
phyl_edge					 = as.data.frame(phyl$edge)
colnames(phyl_edge) = c("edge_start","edge_end")
phyl_edge$ix_edges = as.numeric(rownames(phyl_edge))
phyl_edge$ends_in_tip = phyl_edge$edge_end <= length(phyl$tip.label)
# dataframe of nodes
phyl_nods = data.frame(taxa = c(phyl$tip.label, phyl$node.label))
phyl_nods$edge_end = as.numeric(rownames(phyl_nods))
phyl_nods$is_tip	 = phyl_nods$edge_end <= length(phyl$tip.label)
# merge them
phyl_edge = merge(phyl_edge, phyl_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")

# heatmap col
col_blue = colorRampPalette(interpolate="l",c("gray90", "deepskyblue","dodgerblue3","midnightblue"))


## Read Posteriors ##

# read orthogroup data
message("posteriors | read posteriors data")
ort = as.data.frame(data.table::fread(ort_fn, header = TRUE))
# drop ABSENT entry (useless here)
ort = ort[ort$Family != "ABSENT", ]
ort_taxa = stringr::str_remove(colnames(ort)[grepl("\\:1$", colnames(ort))], "\\:1")
rownames(ort) = ort$Family

# # drop singletons?
# drop_singletons_bool = sapply(ort$Pattern, function(v) sum(as.numeric(unlist(stringr::str_split(v, pattern = " |(|)"))), na.rm = TRUE) == 1 )
# ort = ort [ !drop_singletons_bool, ]

# get presence, gains and loss matrices
message("posteriors | prepare posteriors g/l/e/r matrices")
mat_pre1 = as.matrix(ort[,grepl("\\:1$", colnames(ort))])
mat_prem = as.matrix(ort[,grepl("\\:m$", colnames(ort))])
mat_gain = as.matrix(ort[,grepl("\\:gain$", colnames(ort))])
mat_loss = as.matrix(ort[,grepl("\\:loss$", colnames(ort))])
mat_expa = as.matrix(ort[,grepl("\\:expansion$", colnames(ort))])
mat_redu = as.matrix(ort[,grepl("\\:reduction$", colnames(ort))])
# sum one/more presence probs
mat_pres = mat_pre1 + mat_prem
# add empty ANCESTRAL node to gains and losses
mat_gain = cbind(mat_gain, 0)
mat_loss = cbind(mat_loss, 0)
mat_expa = cbind(mat_expa, 0)
mat_redu = cbind(mat_redu, 0)
# colnames
colnames(mat_pres) = ort_taxa
colnames(mat_pre1) = ort_taxa
colnames(mat_prem) = ort_taxa
colnames(mat_gain) = ort_taxa
colnames(mat_loss) = ort_taxa
colnames(mat_expa) = ort_taxa
colnames(mat_redu) = ort_taxa
# reorder to match species tree
taxa_lis = c(phyl$tip.label, rev(phyl$node.label))
mat_pres = mat_pres[,taxa_lis]
mat_pre1 = mat_pre1[,taxa_lis]
mat_prem = mat_prem[,taxa_lis]
mat_gain = mat_gain[,taxa_lis]
mat_loss = mat_loss[,taxa_lis]
mat_expa = mat_expa[,taxa_lis]
mat_redu = mat_redu[,taxa_lis]

# get summary matrix
mat_summ = ort[,c("Family","Pattern")]
colnames(mat_summ) = c("orthogroup","pattern")

# add summary info
message("posteriors | summary info")
mat_summ$sum_gain = apply(mat_gain, 1, sum)
mat_summ$sum_loss = apply(mat_loss, 1, sum)
mat_summ$sum_redu = apply(mat_redu, 1, sum)
mat_summ$sum_expa = apply(mat_expa, 1, sum)
mat_summ$sum_pres_total = apply(mat_pres, 1, sum)
mat_summ$sum_pres_extant = apply(mat_pres[,colnames(mat_pres) %in% sps_list], 1, sum)

# add og name
mat_summ$name = ort_names_v [ rownames(mat_summ) ]


## Output matrices ##

# create presence, gain and loss matrices per node
message("posteriors | write matrices")
write.table(mat_pres, sprintf("%s.pres.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_pre1, sprintf("%s.pres_singlecopy.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_prem, sprintf("%s.pres_multicopy.csv",  gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_loss, sprintf("%s.loss.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_gain, sprintf("%s.gain.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_redu, sprintf("%s.redu.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_expa, sprintf("%s.expa.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)
write.table(mat_summ, sprintf("%s.summ.csv", gsub("\\.csv$","",ort_fn)), sep = "\t", quote = FALSE)


## Plot tree ##

# read species tree
message("posteriors | prepare cladogram")
phy = ape::read.tree(phy_fn)
nod_list = c(phy$tip.label, phy$node.label)

# prepare species dataframe
# dataframe of edges
message("posteriors | prepare cladogram, edges data")
phy_edge = as.data.frame(phy$edge)
colnames(phy_edge) = c("edge_start","edge_end")
phy_edge$ix_edges = as.numeric(rownames(phy_edge))
phy_edge$ends_in_tip = phy_edge$edge_end <= length(phy$tip.label)
# dataframe of nodes
phy_nods = data.frame(taxa = c(phy$tip.label, phy$node.label))
phy_nods$edge_end = as.numeric(rownames(phy_nods))
phy_nods$is_tip   = phy_nods$edge_end <= length(phy$tip.label)
# merge them
phy_edge = merge(phy_edge, phy_nods, all.x = TRUE, all.y = TRUE, by.x = "edge_end", by.y = "edge_end")
# which is the root node?
root_label = phy$node.label[1]

# summarise gains, losses and presences per node, summing probabilities?
message("posteriors | prepare cladogram, summarise metadata")
mss_f = data.frame(
	row.names = nod_list,
	taxa = nod_list,
	pres = colSums(mat_pres,na.rm=TRUE)[nod_list],
	gain = colSums(mat_gain,na.rm=TRUE)[nod_list],
	loss = colSums(mat_loss,na.rm=TRUE)[nod_list],
	expa = colSums(mat_expa,na.rm=TRUE)[nod_list],
	redu = colSums(mat_redu,na.rm=TRUE)[nod_list]
)

# summarise gains, losses and presences per node, with probability threshold?
mss_f = data.frame(
	row.names = nod_list,
	taxa = nod_list,
	pres = colSums(mat_pres>=0.5,na.rm=TRUE)[nod_list],
	gain = colSums(mat_gain>=0.5,na.rm=TRUE)[nod_list],
	loss = colSums(mat_loss>=0.5,na.rm=TRUE)[nod_list],
	expa = colSums(mat_expa>=0.5,na.rm=TRUE)[nod_list],
	redu = colSums(mat_redu>=0.5,na.rm=TRUE)[nod_list]
)


# info to plot
message("posteriors | prepare cladogram, conform")
phy_data = merge(phy_edge, mss_f, by.x = "taxa", by.y = "taxa", all.x = TRUE)
rownames(phy_data) = phy_data$taxa
phy_data = phy_data [ order(phy_data$ix_edges), ]
phy_data = phy_data [ !is.na(phy_data$ix_edges), ]
max_data = 100

# dataframes of gains/losses, to build pies ("pres" category includes all)
message("posteriors | prepare cladogram, dataframes for pies")
gain_df = data.frame(
	row.names = rownames(phy_data),
	gain = phy_data$gain,
	pres = phy_data$pres - phy_data$gain
)
gain_df = gain_df + 1e-9
exre_df = data.frame(
	row.names = rownames(phy_data),
	expa = phy_data$expa,
	redu = phy_data$redu,
	pres = phy_data$pres - phy_data$expa - phy_data$redu
)
exre_df = exre_df + 1e-9
loss_df = data.frame(
	row.names = rownames(phy_data),
	loss = phy_data$loss + 1e-9,
	nothing = 0
)

# plot tree
message("posteriors | prepare cladogram, plot tree")
pdf(sprintf("%s.tree.pdf", gsub(".csv$", "", ort_fn)), width = 12, height = 16)

# plot of gains and losses
ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)
ape::edgelabels(pie = gain_df, piecol = c("olivedrab2","lightcyan2"), cex = sqrt(phy_data$pres / max_data) * 0.1)
ape::edgelabels(pie = loss_df, piecol = c("indianred1","indianred3"), cex = sqrt(phy_data$loss / max_data) * 0.1, adj = 0)
ape::edgelabels(
	text = sprintf("%s\np=%.1f\n+%.1f | -%.1f", phy_data$taxa, phy_data$pres, phy_data$gain, phy_data$loss),
	col = scales::alpha("gray20", 0.9), frame="none", cex = 0.5)
title(sub  = sprintf("root present n=%.1f", mss_f[root_label,"pres"]))
title(main  = "gains and losses")

# plot of increases and contractions
ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)
ape::edgelabels(pie = exre_df, piecol = c("dodgerblue2","orchid2","lightcyan2"), cex = sqrt(phy_data$pres / max_data) * 0.1)
ape::edgelabels(
	text = sprintf("%s\np=%.1f\nup:%.1f | down:%.1f", phy_data$taxa, phy_data$pres, phy_data$expa, phy_data$redu),
	col = scales::alpha("gray20", 0.9), frame="none", cex = 0.5)
title(sub  = sprintf("root present n=%.1f", mss_f[root_label,"pres"]))
title(main  = "within present, expanded and contracted")

# plot of expansions
ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)
ape::edgelabels(pie = cbind(exre_df[,"expa"],0), piecol = c("dodgerblue2","black"), cex = sqrt(phy_data$expa / max_data) * 0.1)
ape::edgelabels(
	text = sprintf("%s\nup:%.1f", phy_data$taxa, phy_data$expa),
	col = scales::alpha("gray20", 0.9), frame="none", cex = 0.5)
title(sub  = sprintf("root present n=%.1f", mss_f[root_label,"pres"]))
title(main  = "expanded")

# plot of contractions
ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1)
ape::edgelabels(pie = cbind(exre_df[,"redu"],0), piecol = c("orchid2","black"), cex = sqrt(phy_data$redu / max_data) * 0.1)
ape::edgelabels(
	text = sprintf("%s\ndown:%.1f", phy_data$taxa, phy_data$redu),
	col = scales::alpha("gray20", 0.9), frame="none", cex = 0.5)
title(sub  = sprintf("root present n=%.1f", mss_f[root_label,"pres"]))
title(main  = "contracted")

# plot node size legend
ape::plot.phylo(phy, font=1, type="phylogram", label.offset = 0.1, root.edge = TRUE, align.tip.label = TRUE, cex=1, edge.color = FALSE, show.tip.label = FALSE, main = "dotsize legend")
legend_vec = c(5,10,50,100,200,500,1000,2000,5000,10000,20000,50000)
ape::tiplabels(pie = c(rep(1, length(legend_vec)), rep(0, length(phy$tip.label) - length(legend_vec))), piecol = c("violet"), cex = sqrt(legend_vec / max_data) * 0.1, offset = -2)
ape::tiplabels(c(legend_vec, rep(0, length(phy$tip.label) - length(legend_vec))), col = c("darkorchid"), frame = "none", offset = -3)

dev.off()


## Gene ages ##

# write reports of likely gene age per species
message("posteriors | gene ages")
for (spi in c("Ocupat","Ocuarb","Xesp","Nvec","Amil","Spis")) {

	message(sprintf("posteriors | gene ages report for: %s", spi))
	nod_path_i = nod_list[ ape::nodepath(phy, which(nod_list == spi), which(nod_list == root_label))]
	keep_ogs_i = unique(dat [ dat$species == spi, "orthogroup" ])
	keep_ogs_i = keep_ogs_i [ keep_ogs_i %in% rownames(mat_summ) ]
	
	hist_i = data.frame(
		gain_node = apply(mat_gain[keep_ogs_i,nod_path_i], 1, function(v) { nod_path_i[which.max(v)] }), 
		gain_prob = apply(mat_gain[keep_ogs_i,nod_path_i], 1, function(v) { max(v,na.rm=TRUE) })
	)
	hist_i$gain_node = factor(hist_i$gain_node, levels = nod_path_i)
	
	# dollo subset
	dol_i = dol[ keep_ogs_i,]
	
	# which are unreliable and should be assigned to another node based on dollo?
	hist_i_ancestral_boo = hist_i$gain_prob < 0.5
	hist_i [ hist_i_ancestral_boo, "gain_node" ] = dol_i[hist_i_ancestral_boo,"gain"]
	hist_i [ hist_i_ancestral_boo, "gain_prob" ] = 0
	
	# dictionaries
	hist_i_p = dic_from_vecs(rownames(hist_i), hist_i$gain_prob)
	hist_i_v = dic_from_vecs(rownames(hist_i), hist_i$gain_node)
	
	gist_i = data.frame(
		gene = dat[dat$species==spi, "gene"],
		orthogroup = dat[dat$species==spi, "orthogroup"],
		orthogroup_name = dat[dat$species==spi, "orthogroup_name"]
	)
	gist_i$gain_node = hist_i_v [ gist_i$orthogroup ]
	gist_i$gain_prob = hist_i_p [ gist_i$orthogroup ]
	gist_i$gain_node_age = as.numeric(gist_i$gain_node)
	
	write.table(gist_i, sprintf("%s/geneages.%s.posteriors.csv", out_fn, spi), sep = "\t", row.names= FALSE, quote = FALSE)
	
}

message("done!")

