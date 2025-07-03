# libraries
library("rphast")
library("ape")
source("../scripts/helper.R")
# source("scripts/mta_downstream_functions.R")

# reference species
sps_ref  = "Ocupat"
msa_fn = "results_wga/results_alignments_plus/cactus.Ocupat.chroc10.maf"
chr_ref = "chroc10"
# sps_ref  = "Amil"
# msa_fn = "results_wga/results_alignments/cactus.Amil.chr10.maf"
# chr_ref = "chr10"
# sps_ref  = "Spis"
# msa_fn = "results_wga/results_alignments/cactus.Spis.NW_019217785.1.maf"
# chr_ref = "NW_019217785.1"

# files
gix_fn = sprintf("../data/reference/%s_gDNA.fasta.fai", sps_ref)
gtf_fn = sprintf("../data/reference/%s_long.annot.gtf", sps_ref)
out_fn = "results_wga_phastcons_plus/"
dir.create(out_fn, showWarnings = FALSE)

# get chrom lengths
gix = read.table(gix_fn)
gix_v = gix[,2]
names(gix_v) = gix[,1]

# # gff file with TEs
# tes_fn = "results_repeats/results_edta_Tadh/Tadh.fasta.mod.EDTA.TEanno.original_scaffolds.gff"
# tes_r = rtracklayer::import(tes_fn)
# tes_f = rphast::feat(seq = seqnames(tes_r), start = start(tes_r), end = end(tes_r), feature = tes_r$type , attribute = tes_r$ID)

message(sprintf("phast | train based on %s | load training data...", sps_ref))
# species tree
tre_fn = "../data/species_tree.Scleractinia2.newick"
tre = scan(tre_fn, what="character", sep=NULL)
# sps_tree = "((Tadh, TrH2), (Hhon, HoiH23));"

# load alignment
msa = rphast::read.msa(msa_fn)

# load features
fet_f = rphast::read.feat(gtf_fn)

# add introns, get genes and CDSs
fet_f = rphast::add.introns.feat(fet_f)
gen_f = fet_f [ fet_f$feature == "transcript", ]
cds_f = fet_f [ fet_f$feature == "CDS", ]

# select CDS and genes regions from first chromosome
fet_f_i = fet_f [ fet_f$seqname == chr_ref, ]
cds_f_i = cds_f [ cds_f$seqname == chr_ref, ]
gen_f_i = gen_f [ gen_f$seqname == chr_ref, ]

# rename chromosome to match names in msa
fet_f_i$seqname = sps_ref
cds_f_i$seqname = sps_ref
gen_f_i$seqname = sps_ref

# obtain 4-fold degenerated positions, which will be used to initialise the non-conserved model
message(sprintf("phast | train based on %s | get four-fold degenerated sites...", sps_ref))
msa_4fd = rphast::get4d.msa(msa, cds_f_i)

# obtain neutral phyloFit model, based on the four-fold degenerate positions
# REV is the most general model for nucleotide substitution, subject to the time-reversibility constraint. It has four frequencies and five rate parameters.
message(sprintf("phast | train based on %s | init model with phylofit...", sps_ref))
pha_mod_neutral_pf = rphast::phyloFit(msa_4fd, tree = tre, subst.mod = "REV")

# use this initial model to optimise another one with phastcons
# estimate transitions now, so that we can use this for smaller scaffolds, too
message(sprintf("phast | train based on %s | optimise model with expectation maximisation...", sps_ref))
pha_mod_neutral_pc = rphast::phastCons(msa = msa, mod = pha_mod_neutral_pf, estimate.trees = TRUE, quiet = FALSE, viterbi = FALSE, estimate.transitions = TRUE)

message(sprintf("phast | train based on %s | save model...", sps_ref))
saveRDS(pha_mod_neutral_pc, sprintf("%s/model.phastcons.ref%s.rds", out_fn, sps_ref))

