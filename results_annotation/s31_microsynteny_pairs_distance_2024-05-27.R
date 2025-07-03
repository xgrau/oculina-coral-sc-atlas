# libraries
suppressMessages(library("rtracklayer"))
suppressMessages(library("stringr"))
suppressMessages(library("igraph"))
suppressMessages(source("../scripts/helper.R"))
graphics.off()

#### Definitions ####

ogp_fn = "../data/orthology_Anthozoa_plus/orthogroup_pairs.transcripts.csv" # orthology pairs
gff_fo = "../data/reference/" # where can I find gff files?
sps_lis = c("Ocupat","Ocuarb","Spis","Amil","Nvec","Xesp") # species list
out_fn = "results_microsynteny_plus"
dir.create(out_fn, showWarnings = FALSE)

# parsing options
max_dis = 3 # how far away should we go when looking for pairs?


#### Load orthology ####

# orthology
ogp = read.table(ogp_fn, header = FALSE, sep="\t", col.names = c("gene1","gene2"), stringsAsFactors = FALSE)
# add species info
ogp$species1 = stringr::str_split(ogp$gene1, pattern = "_", simplify = TRUE)[,1]
ogp$species2 = stringr::str_split(ogp$gene2, pattern = "_", simplify = TRUE)[,1]
# ogp$species1 [ ogp$species1 == "Ocupat" ] = "Ocubis"
# ogp$species2 [ ogp$species2 == "Ocupat" ] = "Ocubis"

min_ogs_per_chr = 100

# loop ref sps
for (spi in sps_lis) {
    
    dat = data.frame()
    
    # loop other sps
    for (spj in sps_lis) {
        
        if (spi != spj) {
            
            # load og pairs
            ogp_i = ogp [ ogp$species1 %in% c(spi,spj) & ogp$species2 %in% c(spi,spj), ]
            message(sprintf("# compare chrs %s-%s", spi, spj))

            # load synteny pairs
            dpa = read.table(file = sprintf("%s/pairs_%s-%s.csv", out_fn, spi, spj), sep="\t", header = TRUE)
            dpa = dpa [ dpa$gene_distance_i <= max_dis ,]
            dpa$seqid_i = factor(dpa$seqid_i, levels = stringr::str_sort(unique(dpa$seqid_i), numeric = TRUE))
            dpa$seqid_j = factor(dpa$seqid_j, levels = stringr::str_sort(unique(dpa$seqid_j), numeric = TRUE))
            
            # which genes from sps i are in the orthology table?
            spi_genes_in_ogp = unique(c(ogp_i$gene1 [ ogp_i$species1 == spi ], ogp_i$gene2 [ ogp_i$species2 == spi ]))
            spj_genes_in_ogp = unique(c(ogp_i$gene1 [ ogp_i$species1 == spj ], ogp_i$gene2 [ ogp_i$species2 == spj ]))
            
            # which genes from sps i are in collinear pairs with sps j?
            spi_genes_in_dpa = unique(c(dpa$transcript_a_i, dpa$transcript_b_i))
            
            # store
            dat = rbind(dat, data.frame(spi = spi, spj = spj, num_ogs_shared = length(spi_genes_in_ogp), num_ogs_in_synteny = length(spi_genes_in_dpa)))
            
            # min ogs per chr
            if (spi == "Spis") {
                min_ogs_per_chr_i = 5
            } else {
                min_ogs_per_chr_i = min_ogs_per_chr
            }
            if (spj == "Spis") {
                min_ogs_per_chr_j = 5
            } else {
                min_ogs_per_chr_j = min_ogs_per_chr
            }
            
            # crosstabulate chromosomes
            xtb = table(dpa$seqid_i, dpa$seqid_j)
            xtb = xtb [ rowSums(xtb) >= min_ogs_per_chr_i , ]
            xtb = xtb [ , colSums(xtb) >= min_ogs_per_chr_j ]
            
            xtf = xtb / rowSums(xtb)
            xtf [ is.na(xtf) ] = 0
            order_cols = names(sort(apply(xtf,1,function(x) which.max(x))))
            xtf = xtf[order_cols,]
            top_chj_per_chi = colnames(xtf) [ apply(xtf, 1, which.max) ]
            
            chi_colors = viridisLite::turbo(n = ncol(xtf), begin = 0.05, end = 0.95)
            names(chi_colors) = colnames(xtf)
            pp = plot_complex_heatmap(
                xtf,
                color_mat = c("gray95", "skyblue", "dodgerblue3", "midnightblue"),
                cluster_row = FALSE,
                cluster_col = FALSE,
                use_raster = FALSE,
                categories_row = top_chj_per_chi,
                colors_row = chi_colors,
                categories_col = colnames(xtf),
                colors_col = chi_colors,
                cell_border = gpar(col = "white", lwd = 1, lty = 1),
                heatmap_border = gpar(col = "black", lwd = 1, lty = 1)
            )
            pdf(sprintf("%s/pairs_%s-%s.pdf", out_fn, spi, spj), height = (dim(xtf)[1] / 10) + 5, width = (dim(xtf)[2] / 10) + 5)
            print(pp)
	        dev.off()
            
        }
    }
    
    
    
    # plots
    dat$fraction_ogs_in_synteny = dat$num_ogs_in_synteny / dat$num_ogs_shared
    
    pdf(sprintf("%s/shared_pairs_%s.pdf", out_fn, spi), width = 3, height = 4)
    b=barplot(
        dat$fraction_ogs_in_synteny, 
        col = "lightblue2", 
        names = dat$spj, 
        ylim = c(0,1),
        las = 2, 
        main = sprintf("fraction of %s orthologs in synteny at d=%i", spi, max_dis),
        cex.main = 0.6, font.main = 1)
    text(b,0, sprintf("%.3f", dat$fraction_ogs_in_synteny), srt = 90, adj = 0)
    b=barplot(
        dat$num_ogs_in_synteny, 
        col = "lightblue2", 
        names = dat$spj, las = 2, 
        main = sprintf("num of %s orthologs in synteny at d=%i", spi, max_dis),
        cex.main = 0.6, font.main = 1)
    text(b,0, sprintf("%i", dat$num_ogs_in_synteny), srt = 90, adj = 0)
    b=barplot(
        dat$num_ogs_shared, 
        col = "lightblue2", 
        names = dat$spj, las = 2, 
        main = sprintf("num of %s orthologs shared", spi),
        cex.main = 0.6, font.main = 1)
    text(b,0, sprintf("%i", dat$num_ogs_shared), srt = 90, adj = 0)
    dev.off()
    
}