# load libs
library("stringr")
library("ape")
source("../scripts/gene-set-analysis.R")

# input
ort_fn = "../data/orthology_Anthozoa_plus/orthogroup_conservation.csv"
out_fn = "results_gene_family_evolution_anthozoa_plus/"
phy_fn = "../data/species_tree.Anthozoa_plus.newick"
dir.create(out_fn, showWarnings = FALSE)

# read phylogeny data
phy = read.tree(phy_fn)
sps_list = phy$tip.label


#### PFAM-based training dataset ####

# input
# dat_fl = list.files("../data/reference/", pattern = "pep.pfamscan_archs.csv", full.names = TRUE)
# dat = data.frame()
# for (fi in dat_fl) {
#   dai = read.table(fi, header = F, col.names = c("gene","pfamscan"), sep = "\t")  
#   dat_i_l = gsa_enrichment_load_pfam_list(fi, do_unique = TRUE)
#   dat_i_d = data.frame(
#     gene = unlist(sapply(1:length(dat_i_l), function(vv) { rep(names(dat_i_l)[vv], lengths(dat_i_l)[vv] )  } )),
#     pfam = unlist(dat_i_l)
#   )
#   dat = rbind(dat, dat_i_d)
# }

dat = data.frame()

for (spi in sps_list) {
  fi = sprintf("../data/reference/%s_long.pep.pfamscan_archs.csv", spi)
  message(fi, " ", file.exists(fi), " ",spi)
  dai = read.table(fi, header = F, col.names = c("gene","pfamscan"), sep = "\t") 
    dat_i_l = gsa_enrichment_load_pfam_list(fi, do_unique = TRUE)
  dat_i_d = data.frame(
    gene = unlist(sapply(1:length(dat_i_l), function(vv) { rep(names(dat_i_l)[vv], lengths(dat_i_l)[vv] )  } )),
    pfam = unlist(dat_i_l)
  )
  dat = rbind(dat, dat_i_d)
 
}


dat$species = stringr::str_split(dat$gene, pattern = "_", simplify = TRUE)[,1]
dat$species = factor(dat$species, levels = sps_list)

# output counts (all)
ort_counts = table(orthogroup = dat$pfam, species = dat$species)
write.table(ort_counts, file = sprintf("%s/pfam_domain_counts.csv", out_fn), sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
ort_counts = read.table(sprintf("%s/pfam_domain_counts.csv", out_fn))

# keep only domains present in at least 5% of the dataset
filt_rows = apply(ort_counts, 1, function(x) sum(x>0)) > ncol(ort_counts) * 0.05
ort_counts_filt = ort_counts[filt_rows,]

# get 1000 random rows
set.seed(1)
ort_counts_train = ort_counts_filt[sample(nrow(ort_counts_filt), 2000, replace = FALSE),]
# ort_counts_train[ort_counts_train>1] = 1

# output training set
ott_fn = sprintf("%s/pfam_domain_counts.train.csv", out_fn)
write.table(ort_counts_train, file = ott_fn, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
system(command = sprintf("sed -i 1's/^/\t/' %s", ott_fn))



#### POSSVM output ####

# read orthogroup data
ort = read.table(ort_fn, header = TRUE, sep = "\t")
ort$species = factor(ort$species, levels = sps_list)

ort_counts = table(orthogroup = ort$orthogroup, species = ort$species)

# output counts (all)
write.table(ort_counts, file = sprintf("%s/orthogroup_counts.csv", out_fn), sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
system(command = sprintf("sed -i 1's/^/\t/' %s/orthogroup_counts.csv", out_fn))

