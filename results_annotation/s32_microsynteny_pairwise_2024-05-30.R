# libraries
suppressMessages(library("rtracklayer"))
suppressMessages(library("stringr"))
suppressMessages(library("igraph"))
suppressMessages(source("../scripts/helper.R"))
suppressMessages(source("../scripts/Cross_species_functions.R"))
graphics.off()

### Definitions ###

ogp_fn = "../data/orthology_Anthozoa_plus/orthogroup_pairs.transcripts.csv" # orthology pairs
gff_fo = "../data/reference/" # where can I find gff files?
# sps_lis = c("Ocupat","Amil","Nvec","Metsen","Dialin") # species list
sps_lis = c("Ocupat","Ocuarb","Amil","Nvec","Xesp") # species list
out_fn = "results_microsynteny_plus"
dir.create(out_fn, showWarnings = FALSE)


### Loop over sps ###

for (spi in sps_lis) {
    
    # loop other sps
    for (spj in sps_lis) {
        
        if (spi != spj) {
            
            # load GFF
            message(sprintf("synteny | %s - %s load...", spi, spj))
            gi = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo, spi))
            gj = rtracklayer::readGFF(sprintf("%s/%s_long.annot.gtf", gff_fo, spj))
            gig = gi[gi$type == "transcript",]
            gjg = gj[gj$type == "transcript",]

            # drop chromosomes with few genes
            message(sprintf("synteny | %s - %s ignore small chromosomes...", spi, spj))
            gig_x = table(gig$seqid)
            gjg_x = table(gjg$seqid)
            gig_x_keep = names(gig_x) [ gig_x > 500 ]
            gjg_x_keep = names(gjg_x) [ gjg_x > 500 ]
            gig = gig [ gig$seqid %in% gig_x_keep, ]
            gjg = gjg [ gjg$seqid %in% gjg_x_keep, ]
            # sort
            gig$seqid = factor(gig$seqid, levels = stringr::str_sort(unique(gig$seqid), numeric = TRUE))
            gjg$seqid = factor(gjg$seqid, levels = stringr::str_sort(unique(gjg$seqid), numeric = TRUE))
            gig = gig [ order(gig$seqid, gig$start), ]
            gjg = gjg [ order(gjg$seqid, gjg$start), ]
            
            # chr sizes in spj
            gjx = read.table(sprintf("%s/%s_gDNA.fasta.fai", gff_fo, spj))
            gjx_v = dic_from_vecs(gjx[,1], gjx[,2])

            # load og pairs
            message(sprintf("synteny | %s - %s get OG pairs...", spi, spj))
            og_pairs_ij = clean_og_pairs(
                og_pairs_fn = ogp_fn, 
                sp1 = spi, 
                sp2 = spj, 
                t2g = FALSE,
                header = TRUE
            )
            
            # retain o2o only
            message(sprintf("synteny | %s - %s clean OG pairs...", spi, spj))
            list_o2o_sp1 = names(which(table(og_pairs_ij[,1]) == 1))
	        list_o2o_sp2 = names(which(table(og_pairs_ij[,2]) == 1))
	        bool_o2o = og_pairs_ij[,1] %in% list_o2o_sp1 & og_pairs_ij[,2] %in% list_o2o_sp2
	        og_pairs_o2o = og_pairs_ij [ bool_o2o, ]

            # assign places to OG pairs
            message(sprintf("synteny | %s - %s place OG pairs...", spi, spj))
            ogc = og_pairs_o2o
            ogc$pair = sprintf("cp%05d", 1:nrow(ogc))
            ogc = merge(ogc, gig[,c("seqid","start","end","transcript_id")], all.x = FALSE, all.y = FALSE, by.x = "sp1", by.y = "transcript_id")
            ogc = merge(ogc, gjg[,c("seqid","start","end","transcript_id")], all.x = FALSE, all.y = FALSE, by.x = "sp2", by.y = "transcript_id", suffixes = c("_i","_j"))
            ogc = ogc [ order(ogc$seqid_i, ogc$start_i), ]
            
            # colors for spi chromosomes
            message(sprintf("synteny | %s - %s plot...", spi, spj))
            cols_per_chr_i = colorspace::darken(viridisLite::turbo(nlevels(ogc$seqid_i), begin = 0.1, end = 0.9), 0.1)
            names(cols_per_chr_i) = levels(ogc$seqid_i)
            
            # plot painted chromosomes
            pdf(sprintf("%s/synteny.%s-%s.pairs.pdf", out_fn, spi, spj), width = 16, height = 8)
            max_height = ceiling(max(gjx_v)/1e6+2)
            # empty plot
            plot(NA, NA, xlim = c(0,30), ylim = c(0,max_height) , las = 1, xlab = sprintf("Chromosomes %s", spj), ylab = "Mb")
            # decorations
            text(x = 1:nlevels(ogc$seqid_j), y = rep(max_height-1, nlevels(ogc$seqid_j)), levels(ogc$seqid_j), srt = 90, cex = 0.7, col = "gray5")
            title(main = sprintf("%s chroms painted according to orthologs in %s", spj, spi), cex.main = 1)
            title(sub = sprintf("n=%i o2o pairs", nrow(ogc)))
            legend("topright", names(cols_per_chr_i), fill = cols_per_chr_i, bty = "n", cex = 0.7, title = sprintf("chr %s", spi))
            # loop over chrs
            for (nn in 1:nlevels(ogc$seqid_j)) {
                
                chj = levels(ogc$seqid_j)[nn]
                ogc_i = ogc [ ogc$seqid_j == chj, ]
                rect(xleft = nn - 0.25, xright = nn + 0.25, ybottom = 0, ytop = gjx_v[chj]/1e6, col = "gray95")
                rect(
                    xleft = rep(nn - 0.25, nrow(ogc_i)),
                    xright = rep(nn + 0.25, nrow(ogc_i)),
                    ybottom = ogc_i$start_j/1e6, 
                    ytop = ogc_i$end_j/1e6,
                    col = cols_per_chr_i[as.character(ogc_i$seqid_i)],
                    border = cols_per_chr_i[as.character(ogc_i$seqid_i)],
                    lwd = 1)
            }
            
            # plot stacked barplots
            layout(matrix(1:16, nrow = 4, byrow = TRUE))
            n_genes_per_block = 25
            for (nn in 1:nlevels(ogc$seqid_j)) {
                
                chj = levels(ogc$seqid_j)[nn]
                ogc_i = ogc [ ogc$seqid_j == chj, ]
                ogc_i = ogc_i [ order(ogc_i$start_j), ]
                ogc_i$order_j = 1:nrow(ogc_i)
                ogc_i$region_j = as.integer(ogc_i$order_j / n_genes_per_block)
                ogc_i_xt = table(ogc_i$region_j, ogc_i$seqid_i)
                ogc_i_xtf = ogc_i_xt / rowSums(ogc_i_xt)
                barplot(t(ogc_i_xtf), col = cols_per_chr_i, space = 0, xlim = c(0,60), las = 2)
                title(main = chj, cex.main = 1, font.main = 1)
                if (nn == 1) {
                    title(sub = sprintf("%i genes per block", n_genes_per_block), cex.sub = 0.7)
                }

            }
            dev.off()
            
            # plot composition per chromosome
            pdf(sprintf("%s/synteny.%s-%s.pies.pdf", out_fn, spi, spj), width = 4, height = 4)
            for (nn in 1:nlevels(ogc$seqid_j)) {
                
                chj = levels(ogc$seqid_j)[nn]
                ogc_i = ogc [ ogc$seqid_j == chj, ]
                ogc_x = table(ogc_i$seqid_i)
                ogc_xf = ogc_x / sum(ogc_x)
                ogc_lab = sprintf("%s|%i|%.3f", names(ogc_x), ogc_x, ogc_xf)
                pie(
                    ogc_x[ogc_x>0],
                    col = cols_per_chr_i[ogc_x>0],
                    labels = ogc_lab[ogc_x>0],
                    main = sprintf("%s:%s\npainted after %s", chj, spj, spi),
                    sub = sprintf("n=%i o2o pairs", sum(ogc_x))
                )

            }
            dev.off()
            
            
        }
    }
}

message("all done!")

